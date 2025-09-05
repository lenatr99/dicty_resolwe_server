#!/usr/bin/env python3

import os
import sys
import time
import gzip
from typing import Dict, List, Optional, Tuple
import concurrent.futures

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import resdk


# Ensure Django settings are available for Resolwe SDK
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")


def connect() -> resdk.Resolwe:
    res = resdk.Resolwe(url="http://localhost:8000", username="admin", password="admin")
    res.login(username="admin", password="admin")
    return res


def pick_target_collection(res: resdk.Resolwe, name: str = "Single cell collection") -> Optional[resdk.resources.collection.Collection]:
    candidates = list(res.collection.filter(name=name))
    if not candidates:
        try:
            all_cols = list(res.collection.filter())
        except Exception:
            all_cols = []
        filtered = []
        for c in all_cols:
            name_match = (c.name or "").lower().strip() == name.lower().strip()
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


def find_existing_time_relation(res: resdk.Resolwe, collection):
    rels = list(res.relation.filter(collection=collection.id))
    # Prefer category containing "time" (case-insensitive)
    rels_time = [r for r in rels if "time" in (getattr(r, "category", "") or "").lower()]
    rels_time.sort(key=lambda r: getattr(r, "id", 0), reverse=True)
    return rels_time[0] if rels_time else None


def delete_old_replicate_data(res: resdk.Resolwe, collection, keep_data_ids: List[int]) -> None:
    data_qs = list(res.data.filter(collection=collection.id))
    keep_set = set(keep_data_ids)
    to_delete = []
    for d in data_qs:
        name = (getattr(d, "name", "") or "").strip()
        # Keep only our new averaged timepoint datasets labeled as "Expression hrXX"
        if d.id in keep_set:
            continue
        if name.startswith("Expression hr"):  # old averaged or per-replicate timepoint uploads
            to_delete.append(d)
        elif name.startswith("Expression ") and any(tag in name for tag in ["00hr", "04hr", "08hr", "12hr", "16hr", "20hr"]):
            to_delete.append(d)
    # Perform deletions (forced, non-interactive)
    def _force_delete(dobj):
        try:
            print(f"Deleting old Data ID {dobj.id} (name: {dobj.name})")
            dobj.delete(force=True)
        except Exception as e:
            print(f"Failed to delete Data ID {dobj.id}: {e}")
    if to_delete:
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
            list(ex.map(_force_delete, to_delete))


def main():
    print("Connecting to Resolwe…", flush=True)
    res = connect()

    collection = pick_target_collection(res)
    if not collection:
        print("No target collection found.")
        sys.exit(1)
    print(f"Using collection: {collection.name} (ID: {collection.id})", flush=True)

    # Load single-cell matrix
    adata_path = os.path.join(os.path.dirname(__file__), "AX4.h5ad")
    print(f"Loading AnnData from {adata_path}…", flush=True)
    adata = sc.read_h5ad(adata_path)

    # Normalize time labels to hrXX strings
    if "time" not in adata.obs.columns:
        print("AnnData missing 'time' in obs.")
        sys.exit(2)
    # Convert possible formats like "04hr" or "4hr" to hr04 etc.
    obs_time = adata.obs["time"].astype(str)
    obs_time = obs_time.str.extract(r"(\d+)", expand=False).fillna("0")
    obs_time = obs_time.str.zfill(2)
    time_labels = obs_time.apply(lambda s: f"hr{s}")

    # Compute per-timepoint mean for each gene over all cells
    X = adata.layers.get(None, adata.X)
    gene_ids = list(adata.var["gene_ids"])  # type: ignore[index]
    unique_times = sorted(time_labels.unique(), key=lambda t: int(t.replace("hr", "")))

    upload_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "upload"))
    os.makedirs(upload_dir, exist_ok=True)

    tasks: List[Tuple[str, str]] = []  # (file_path, exp_name)
    for lbl in unique_times:
        mask = time_labels == lbl
        if mask.sum() == 0:
            continue
        sub = X[mask.values, :]
        if sp.issparse(sub):
            mean_vals = np.asarray(sub.mean(axis=0)).ravel()
        else:
            mean_vals = np.asarray(sub).mean(axis=0).ravel()
        df = pd.DataFrame({"Gene": gene_ids, "Expression": mean_vals})
        file_path = os.path.join(upload_dir, f"AX4_{lbl}.tab.gz")
        df.to_csv(file_path, sep="\t", index=False, compression="gzip")
        exp_name = f"Expression {lbl}"
        tasks.append((file_path, exp_name))

    # Upload averaged expressions per timepoint
    print(f"Prepared {len(tasks)} averaged timepoint files. Uploading…", flush=True)
    uploaded: Dict[str, resdk.resources.data.Data] = {}

    def run_upload(task: Tuple[str, str]):
        file_path, exp_name = task
        d = upload_expression(res, collection, file_path, exp_name)
        return exp_name, d

    if tasks:
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as ex:
            futs = [ex.submit(run_upload, t) for t in tasks]
            for f in concurrent.futures.as_completed(futs):
                try:
                    exp_name, d = f.result()
                    if d is not None:
                        uploaded[exp_name] = d
                except Exception as e:
                    print(f"Upload failed: {e}")

    # Build time-series relation with a single partition per time label
    if not uploaded:
        print("No averaged datasets uploaded; aborting relation rebuild.")
        sys.exit(3)

    has_entity, _, part_key = compatibility_layer(res)
    partitions = []
    sorted_labels = sorted(uploaded.keys(), key=lambda nm: int(nm.replace("Expression ", "").replace("hr", "")))
    for position, exp_name in enumerate(sorted_labels, start=1):
        d = uploaded[exp_name]
        lbl = exp_name.replace("Expression ", "")
        parent = get_or_make_entity_for_data(d, lbl, None, collection)
        partitions.append({part_key: parent.id, "label": lbl, "position": position})

    # Delete existing time-series relations and create a new clean one
    existing_rel = find_existing_time_relation(res, collection)
    if existing_rel is not None:
        try:
            print(f"Deleting existing time-series relation {existing_rel.slug} (ID {existing_rel.id})")
            existing_rel.delete(force=True)
        except Exception as e:
            print(f"Failed to delete existing relation: {e}")

    # Create relation without '-avg' suffix; if slug collides, retry letting server generate slug
    try:
        relation = res.relation.create(
            type="series",
            category="Time series",
            collection=collection,
            slug=f"relation-{collection.id}",
            partitions=partitions,
        )
    except Exception as e:
        print(f"Slug 'relation-{collection.id}' may exist, retrying without explicit slug: {e}")
        relation = res.relation.create(
            type="series",
            category="Time series",
            collection=collection,
            partitions=partitions,
        )
    relation.save()
    print(f"Created relation '{relation.slug}' with {len(partitions)} partitions.")

    # Delete old replicate data objects
    keep_ids = [d.id for d in uploaded.values()]
    delete_old_replicate_data(res, collection, keep_ids)


if __name__ == "__main__":
    main()


