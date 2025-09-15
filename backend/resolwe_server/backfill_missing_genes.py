#!/usr/bin/env python3

import os
import sys
import time
import gzip
from typing import Dict, List, Optional, Tuple
import concurrent.futures

import numpy as np
import scanpy as sc
import resdk


os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")


def connect() -> resdk.Resolwe:
    import os
    admin_user = os.environ.get('ADMIN_USERNAME', 'admin')
    admin_pass = os.environ.get('ADMIN_PASSWORD', 'admin')
    url = os.environ.get('RESOLWE_URL', 'http://localhost:8000')
    
    res = resdk.Resolwe(url=url, username=admin_user, password=admin_pass)
    res.login(username=admin_user, password=admin_pass)
    return res


def pick_target_collection(res: resdk.Resolwe) -> Optional[resdk.resources.collection.Collection]:
    candidates = list(res.collection.filter(name="Single cell collection"))
    if not candidates:
        try:
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


def compatibility_layer(res: resdk.Resolwe):
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
    print(f"Uploaded {exp_name} (Data ID: {d.id})")
    return d


def find_single_cell_relation(res: resdk.Resolwe, collection):
    # Prefer rebuilt slug if present
    preferred_slug = f"relation-{collection.id}-sc"
    rels = list(res.relation.filter(collection=collection.id))
    target = None
    for r in rels:
        if getattr(r, "slug", None) == preferred_slug:
            target = r
            break
    if target is None:
        # Fallback: pick latest with category 'Single cell series'
        rels_sc = [r for r in rels if getattr(r, "category", "") == "Single cell series"]
        if not rels_sc:
            return None
        rels_sc.sort(key=lambda r: getattr(r, "id", 0), reverse=True)
        target = rels_sc[0]
    # Ensure partitions are loaded
    target.update()
    return target


def build_gene_map(res: resdk.Resolwe, collection) -> Dict[str, resdk.resources.data.Data]:
    data_qs = list(res.data.filter(collection=collection.id))
    gene_map: Dict[str, resdk.resources.data.Data] = {}
    for d in data_qs:
        name = (getattr(d, "name", "") or "").strip()
        if name.startswith("Single cell ") and name.lower() != "single cell time":
            gene = name.replace("Single cell ", "", 1).strip()
            existing = gene_map.get(gene)
            if existing is None or getattr(d, "id", 0) > getattr(existing, "id", 0):
                gene_map[gene] = d
    return gene_map


def main():
    print("Connecting to Resolwe…", flush=True)
    res = connect()

    collection = pick_target_collection(res)
    if not collection:
        print("No target collection found.")
        sys.exit(1)
    print(f"Using collection: {collection.name} (ID: {collection.id})", flush=True)

    relation = find_single_cell_relation(res, collection)
    if relation is None:
        print("No existing single-cell relation found to update.")
        sys.exit(2)
    print(f"Updating relation: {relation.slug} (ID: {relation.id})", flush=True)

    # Existing labels in relation
    existing_labels = [p.get("label") for p in (relation.partitions or []) if p.get("label")]
    existing_set = set(existing_labels)
    print(f"Existing gene partitions: {len(existing_set)}", flush=True)

    # Load data and determine missing genes
    adata_path = os.path.join(os.path.dirname(__file__), "AX4.h5ad")
    print("Loading AnnData…", flush=True)
    adata = sc.read_h5ad(adata_path)
    genes = list(adata.var["gene_ids"])  # type: ignore[index]
    gene_index = {g: i for i, g in enumerate(genes)}
    missing_genes = [g for g in genes if g not in existing_set]
    print(f"Missing genes to backfill: {len(missing_genes)}", flush=True)
    if not missing_genes:
        print("No missing genes found; nothing to do.")
        return

    upload_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "upload"))

    # Discover existing Data objects to reuse when available
    gene_to_data = build_gene_map(res, collection)

    # Prepare uploads for those missing genes that don't already have a Data object
    to_upload: List[str] = [g for g in missing_genes if g not in gene_to_data]
    print(f"Genes needing upload: {len(to_upload)}", flush=True)

    # Create temp files for uploads (sequential I/O)
    total_up = len(to_upload)
    for idx, gene in enumerate(to_upload, start=1):
        gi = gene_index.get(gene)
        if gi is None:
            print(f"Gene not found in AnnData (skipping): {gene}")
            continue
        X = adata.X
        try:
            col = X[:, gi]
            values = col.toarray().ravel() if hasattr(col, "toarray") else np.asarray(col).ravel()
        except Exception:
            values = np.asarray(X[:, gi]).ravel()
        rows = [f"{adata.obs_names[i]}\t{values[i]}" for i in range(values.shape[0])]
        out_path = os.path.join(upload_dir, f"temp_{gene}.tab.gz")
        write_gzip_table(out_path, header="Gene\tExpression", rows=rows)
        if idx % 100 == 0 or idx == total_up:
            print(f"Prepared files: {idx}/{total_up}", flush=True)

    # Parallel uploads
    WORKERS = 10
    upload_tasks: List[Tuple[str, str, int, int]] = []
    for idx, gene in enumerate(to_upload, start=1):
        out_path = os.path.join(upload_dir, f"temp_{gene}.tab.gz")
        exp_name = f"Single cell {gene}"
        upload_tasks.append((out_path, exp_name, idx, total_up))

    def upload_task_runner(task: Tuple[str, str, int, int]):
        file_path, exp_name, ind, total = task
        print(f"⏳ Uploading {ind}/{total}: {exp_name}", flush=True)
        return exp_name, upload_expression(res, collection, file_path, exp_name)

    uploaded_map: Dict[str, resdk.resources.data.Data] = {}
    if upload_tasks:
        with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as executor:
            futures = [executor.submit(upload_task_runner, t) for t in upload_tasks]
            for f in concurrent.futures.as_completed(futures):
                try:
                    exp_name, data_obj = f.result()
                    if data_obj is not None:
                        uploaded_map[exp_name] = data_obj
                except Exception as exc:
                    print(f"Upload generated exception: {exc}")

    # Build new partitions to append (combine reused and uploaded)
    _, _, part_key = compatibility_layer(res)
    new_partitions: List[Dict] = []

    total = len(missing_genes)
    for idx, gene in enumerate(missing_genes, start=1):
        d = gene_to_data.get(gene) or uploaded_map.get(f"Single cell {gene}")
        if d is None:
            # Could not upload or not found
            continue
        parent = get_or_make_entity_for_data(d, gene, 1, collection)
        new_partitions.append({part_key: parent.id, "label": gene})

        if idx % 50 == 0 or idx == total:
            print(f"Processed {idx}/{total} missing genes…", flush=True)

    if not new_partitions:
        print("No new partitions to add.")
        return

    # Compute positions: continue after current max
    current_parts = relation.partitions or []
    current_len = len(current_parts)
    for offset, part in enumerate(new_partitions, start=1):
        part["position"] = current_len + offset

    # Append and save
    relation.partitions = current_parts + new_partitions
    relation.save()
    print(f"Appended {len(new_partitions)} gene partitions. New total: {len(relation.partitions)}", flush=True)


if __name__ == "__main__":
    main()


