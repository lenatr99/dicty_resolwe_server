#!/usr/bin/env python3

import os
import re
import sys
import time
from typing import Dict, List, Optional, Tuple
import gzip

import numpy as np
import scanpy as sc

import resdk


os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")


def connect() -> resdk.Resolwe:
    res = resdk.Resolwe(url="http://localhost:8000", username="admin", password="admin")
    res.login(username="admin", password="admin")
    return res


def pick_target_collection(res: resdk.Resolwe) -> Optional[resdk.resources.collection.Collection]:
    candidates = list(res.collection.filter(name="Single cell collection"))
    if not candidates:
        try:
            # Last-resort fallback: scan and pick most recent matching name/tag
            all_cols = list(res.collection.filter())
        except Exception:
            all_cols = []
        filtered = []
        for c in all_cols:
            name_match = (c.name or "").lower().strip() == "single cell collection"
            tag_match = any((t or "").lower().strip() == "community:bcm" for t in (c.tags or []))
            if name_match or tag_match:
                filtered.append(c)
        candidates = filtered
    if not candidates:
        return None
    candidates.sort(key=lambda c: getattr(c, "created", None) or 0, reverse=True)
    return candidates[0]


# Compatibility layer for Resolwe versions with entity vs sample
def compatibility_layer(res: resdk.Resolwe):
    # Force using 'entity' partitions. Resource falls back to sample only for parent creation if entity API missing.
    has_entity = hasattr(res, "entity")
    resource = res.entity if has_entity else res.sample
    part_key = "entity"
    return has_entity, resource, part_key


def ensure_in_collection(obj, collection):
    has_entity, _, _ = compatibility_layer(obj.resolwe)
    if has_entity:
        current = [c.id for c in (obj.collections or [])]
        if collection.id not in current:
            obj.collections = (obj.collections or []) + [collection]
            obj.save()
    else:
        if getattr(obj, "collection", None):
            if obj.collection.id != collection.id:
                obj.collection = collection
                obj.save()
        else:
            obj.collection = collection
            obj.save()


def get_or_make_entity_for_data(d, label: str, rep_idx: Optional[int], collection):
    has_entity, resource, _ = compatibility_layer(d.resolwe)
    parent = getattr(d, "entity", None) or getattr(d, "sample", None)
    if parent:
        ensure_in_collection(parent, collection)
        return parent
    if rep_idx is None:
        created = resource.create(name=label, data=[d], collections=[collection])
    else:
        created = resource.create(name=f"{label}_r{rep_idx}", data=[d], collections=[collection])
    return created


def collect_data_by_name(
    res: resdk.Resolwe, collection
) -> Tuple[Dict[str, List[resdk.resources.data.Data]], List[Tuple[str, resdk.resources.data.Data]], List[resdk.resources.data.Data], Dict[str, resdk.resources.data.Data]]:
    print("Scanning Data objects in collection…", flush=True)
    data_qs = list(res.data.filter(collection=collection.id))
    print(f"Found {len(data_qs)} Data objects.", flush=True)

    time_bucket: Dict[str, List[resdk.resources.data.Data]] = {}
    umap_list: List[Tuple[str, resdk.resources.data.Data]] = []
    time_var_list: List[resdk.resources.data.Data] = []
    gene_map: Dict[str, resdk.resources.data.Data] = {}

    time_re = re.compile(r"^Expression\s+(\d{2}hr)\b", re.IGNORECASE)

    for idx, d in enumerate(data_qs, start=1):
        name = (d.name or "").strip()
        m = time_re.match(name)
        if m:
            t_raw = m.group(1)  # e.g. 00hr
            label = f"hr{t_raw[:2]}"
            time_bucket.setdefault(label, []).append(d)
            if idx % 100 == 0:
                print(f"  Processed {idx}/{len(data_qs)}…", flush=True)
            continue

        if name.startswith("X_uce_umap_"):
            umap_list.append((name, d))
            continue

        if name.lower() == "single cell time":
            time_var_list.append(d)
            continue

        if name.startswith("Single cell ") and name.lower() != "single cell time":
            gene = name.replace("Single cell ", "", 1).strip()
            # Prefer latest by id if duplicates
            existing = gene_map.get(gene)
            if existing is None or getattr(d, "id", 0) > getattr(existing, "id", 0):
                gene_map[gene] = d
            continue

    print(
        "Buckets: "
        + ", ".join(
            [f"{k}:{len(v)}" for k, v in sorted(time_bucket.items(), key=lambda x: int(x[0][2:]))]
        )
        or "no time buckets",
        flush=True,
    )
    print(f"UMAP axes: {len(umap_list)}", flush=True)
    print(f"Time variable: {len(time_var_list)}", flush=True)
    print(f"Single-gene items: {len(gene_map)}", flush=True)
    return time_bucket, umap_list, time_var_list, gene_map


def write_gzip_table(path: str, header: str, rows: List[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt", compresslevel=1, encoding="utf-8") as f:
        f.write(header + "\n")
        for line in rows:
            f.write(line + "\n")


def upload_expression(res: resdk.Resolwe, collection, file_path: str, exp_name: str) -> Optional[resdk.resources.data.Data]:
    inputs = {
        "exp": file_path,
        "exp_type": "10x_UMI",
        "exp_name": exp_name,
        "source": "DICTYBASE",
        "species": "Dictyostelium discoideum",
        "build": "dd-05-2009",
    }
    try:
        d = res.run(slug="upload-expression", input=inputs, collection=collection)
    except Exception as e:
        print(f"Failed to start upload for {exp_name}: {e}")
        return None
    while d.status not in ["OK", "ER"]:
        time.sleep(2)
        d.update()
    if d.status == "ER":
        try:
            print(f"Upload error for {exp_name}: {d.process_error}")
        except Exception:
            print(f"Upload error for {exp_name}")
        return None
    print(f"Uploaded {exp_name} (Data ID: {d.id})", flush=True)
    return d


def ensure_umap_and_time(res: resdk.Resolwe, collection, existing_umap: List[Tuple[str, resdk.resources.data.Data]], existing_time_vars: List[resdk.resources.data.Data]) -> Tuple[List[Tuple[str, resdk.resources.data.Data]], List[resdk.resources.data.Data]]:
    need_umap = len(existing_umap) < 2
    need_time = len(existing_time_vars) < 1
    if not (need_umap or need_time):
        return existing_umap, existing_time_vars

    print("UMAP/time missing — generating from AX4.h5ad and uploading…", flush=True)
    adata_path = os.path.join(os.path.dirname(__file__), "AX4.h5ad")
    adata = sc.read_h5ad(adata_path)

    upload_dir = os.path.join(os.path.dirname(__file__), "..", "data", "upload")
    upload_dir = os.path.abspath(upload_dir)

    # UMAP axes
    if need_umap:
        if "X_uce_umap" not in adata.obsm.keys():
            raise RuntimeError("'X_uce_umap' not found in AnnData.obsm")
        for i in range(2):
            values = np.asarray(adata.obsm["X_uce_umap"])[:, i].ravel()
            rows = [f"{adata.obs_names[idx]}\t{values[idx]}" for idx in range(values.shape[0])]
            out_path = os.path.join(upload_dir, f"temp_X_uce_umap_{i}.tab.gz")
            write_gzip_table(out_path, header="Gene\tExpression", rows=rows)
            d = upload_expression(res, collection, out_path, f"X_uce_umap_{i}")
            if d is not None:
                existing_umap.append((f"X_uce_umap_{i}", d))

    # Time variable (first two chars of time as int, like test3)
    if need_time:
        if "time" not in adata.obs.columns:
            raise RuntimeError("'time' not found in AnnData.obs")
        time_numeric = adata.obs["time"].astype(str).str[:2].astype(int).to_numpy()
        rows = [f"{adata.obs_names[idx]}\t{time_numeric[idx]}" for idx in range(time_numeric.shape[0])]
        out_path = os.path.join(upload_dir, "temp_time.tab.gz")
        write_gzip_table(out_path, header="Gene\tExpression", rows=rows)
        d = upload_expression(res, collection, out_path, "Single cell time")
        if d is not None:
            existing_time_vars.append(d)

    return existing_umap, existing_time_vars


def build_time_partitions(res: resdk.Resolwe, collection, time_bucket, umap_list, time_var_list):
    _, _, part_key = compatibility_layer(res)

    labels_sorted = sorted(time_bucket.keys(), key=lambda k: int(k[2:]))
    label_to_pos = {lbl: i + 1 for i, lbl in enumerate(labels_sorted)}

    partitions = []
    for lbl in labels_sorted:
        items = time_bucket[lbl]
        total = len(items)
        print(f"Building time partitions for {lbl} ({total})…", flush=True)
        for i, d in enumerate(items, start=1):
            parent = get_or_make_entity_for_data(d, lbl, i, collection)
            partitions.append({part_key: parent.id, "label": lbl, "position": label_to_pos[lbl]})
            if i % 50 == 0 or i == total:
                print(f"  {lbl}: {i}/{total}", flush=True)

    umap_list_sorted = sorted(umap_list, key=lambda x: x[0])
    base_pos = len(labels_sorted)
    print(f"Adding UMAP axes ({len(umap_list_sorted)})…", flush=True)
    for i, (name, d) in enumerate(umap_list_sorted):
        parent = get_or_make_entity_for_data(d, name, None, collection)
        partitions.append({part_key: parent.id, "label": name, "position": base_pos + i + 1})

    print(f"Adding time variable items ({len(time_var_list)})…", flush=True)
    for i, d in enumerate(time_var_list):
        parent = get_or_make_entity_for_data(d, "time", i, collection)
        partitions.append({part_key: parent.id, "label": "time", "position": base_pos + len(umap_list_sorted) + i + 1})

    return partitions


def build_gene_partitions(res: resdk.Resolwe, collection, gene_map):
    _, _, part_key = compatibility_layer(res)
    genes_sorted = sorted(gene_map.keys())
    partitions = []
    total = len(genes_sorted)
    print(f"Building single-gene partitions ({total})…", flush=True)
    for pos, gene in enumerate(genes_sorted, start=1):
        d = gene_map[gene]
        parent = get_or_make_entity_for_data(d, gene, 1, collection)
        partitions.append({part_key: parent.id, "label": gene, "position": pos})
        if pos % 100 == 0 or pos == total:
            print(f"  genes: {pos}/{total}", flush=True)
    return partitions


def main():
    print("Connecting to Resolwe…", flush=True)
    res = connect()

    collection = pick_target_collection(res)
    if not collection:
        print("No target collection found. Ensure your earlier uploads used 'Single cell collection'.")
        sys.exit(1)

    print(f"Using collection: {collection.name} (ID: {collection.id})", flush=True)

    time_bucket, umap_list, time_var_list, gene_map = collect_data_by_name(res, collection)

    if not (time_bucket or umap_list or time_var_list or gene_map):
        print("No matching Data objects found in the target collection.")
        sys.exit(2)

    # Ensure required UMAP/time exist (create by uploading from AX4.h5ad if missing)
    umap_list, time_var_list = ensure_umap_and_time(res, collection, umap_list, time_var_list)

    time_partitions = build_time_partitions(res, collection, time_bucket, umap_list, time_var_list)
    gene_partitions = build_gene_partitions(res, collection, gene_map)

    descriptor = {
        "growth": "HL5",
        "strain": "AX4",
        "details": "AX4",
        "project": "Single cell",
        "citation": {"url": "TO-DO", "name": "Katoh-Kurasawa M et. al."},
        "treatment": "Filter Development",
    }

    ds_qs = res.descriptor_schema.filter(slug="dicty-time-series")
    descriptor_schema = ds_qs[0] if ds_qs else None

    ts_slug = f"relation-{collection.id}-rebuilt"
    sc_slug = f"relation-{collection.id}-rebuilt-sc"

    # Create/overwrite relations with the chosen slugs if they already exist
    try:
        existing = list(res.relation.filter(collection=collection.id, slug=ts_slug))
        for r in existing:
            r.delete()
    except Exception:
        pass
    try:
        existing = list(res.relation.filter(collection=collection.id, slug=sc_slug))
        for r in existing:
            r.delete()
    except Exception:
        pass

    relation = None
    relation2 = None

    if time_partitions:
        print("Creating time-series relation…", flush=True)
        relation = res.relation.create(
            type="series",
            category="Time series",
            collection=collection,
            slug=ts_slug,
            partitions=time_partitions,
        )
        relation.descriptor = descriptor
        if descriptor_schema:
            relation.descriptor_schema = descriptor_schema
        relation.save()
        print(
            f"Created time-series relation '{relation.slug}' (ID: {relation.id}) with {len(time_partitions)} partitions.",
            flush=True,
        )
    else:
        print("No time-series partitions found; skipping time-series relation.", flush=True)

    if gene_partitions:
        print("Creating single-cell relation…", flush=True)
        relation2 = res.relation.create(
            type="series",
            category="Single cell series",
            collection=collection,
            slug=sc_slug,
            partitions=gene_partitions,
        )
        relation2.descriptor = descriptor
        if descriptor_schema:
            relation2.descriptor_schema = descriptor_schema
        relation2.save()
        print(
            f"Created single-cell relation '{relation2.slug}' (ID: {relation2.id}) with {len(gene_partitions)} partitions.",
            flush=True,
        )
    else:
        print("No single-gene partitions found; skipping single-cell relation.", flush=True)

    # Final summary
    summary_parts = []
    if relation is not None:
        summary_parts.append(
            f"time-series '{relation.slug}' ({len(time_partitions)} partitions)"
        )
    if relation2 is not None:
        summary_parts.append(
            f"single-cell '{relation2.slug}' ({len(gene_partitions)} partitions)"
        )
    if summary_parts:
        print("Successfully created: " + "; ".join(summary_parts), flush=True)
    else:
        print("No relations were created (no partitions found).", flush=True)


if __name__ == "__main__":
    main()


