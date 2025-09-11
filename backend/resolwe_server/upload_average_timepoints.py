#!/usr/bin/env python3

import os
import sys
import time
import gzip
from typing import Dict, List, Optional, Tuple
import concurrent.futures
from datetime import datetime

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


def upload_expression(res: resdk.Resolwe, collection, file_path: str, exp_name: str, timestamp: str) -> Optional[resdk.resources.data.Data]:
    # Add timestamp to make names unique
    unique_exp_name = f"{exp_name} {timestamp}"
    inputs = {
        "exp": file_path,
        "exp_type": "10x_UMI",
        "exp_name": unique_exp_name,
        "source": "DICTYBASE",
        "species": "Dictyostelium discoideum",
        "build": "dd-05-2009",
    }
    try:
        d = res.run(slug="upload-expression", input=inputs, collection=collection)
    except Exception as e:
        error_msg = str(e)
        if "duplicate key" in error_msg.lower():
            print(f"Skipping {exp_name} - already exists (duplicate key error)")
            return None
        else:
            print(f"Failed to start upload for {exp_name}: {e}")
            return None
    
    # Wait for completion
    while d.status not in ["OK", "ER"]:
        time.sleep(2)
        d.update()
        
    if d.status == "ER":
        try:
            print(f"Upload error for {exp_name}: {d.process_error}")
        except Exception:
            print(f"Upload error for {exp_name}")
        return None
    print(f"Uploaded {unique_exp_name} (Data ID: {d.id})")
    return d


def main():
    print("Connecting to Resolwe…", flush=True)
    res = connect()

    collection = pick_target_collection(res)
    if not collection:
        print("No target collection found.")
        sys.exit(1)
    print(f"Using collection: {collection.name} (ID: {collection.id})", flush=True)

    # Generate timestamp for unique naming
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"Using timestamp: {timestamp}", flush=True)

    # Create a new collection for averaged data to avoid relation category constraint on the original collection
    new_collection_name = f"{collection.name} Averaged {timestamp}"
    print(f"Creating new collection for averaged data: {new_collection_name}", flush=True)
    try:
        new_collection = res.collection.create(name=new_collection_name)
    except Exception as e:
        print(f"Failed to create '{new_collection_name}': {e}. Retrying with fallback name.")
        fallback_name = f"{new_collection_name} Copy"
        new_collection = res.collection.create(name=fallback_name)
    print(f"Using averaged collection: {new_collection.name} (ID: {new_collection.id})", flush=True)

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
        file_path = os.path.join(upload_dir, f"AX4_{lbl}_avg_{timestamp}.tab.gz")
        df.to_csv(file_path, sep="\t", index=False, compression="gzip")
        exp_name = f"Expression {lbl} Average"
        tasks.append((file_path, exp_name))

    # Upload averaged expressions per timepoint
    print(f"Prepared {len(tasks)} averaged timepoint files. Uploading…", flush=True)
    uploaded: Dict[str, resdk.resources.data.Data] = {}

    def run_upload(task: Tuple[str, str]):
        file_path, exp_name = task
        d = upload_expression(res, new_collection, file_path, exp_name, timestamp)
        return exp_name, d

    if tasks:
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:  # Reduced workers to avoid overwhelming
            futs = [ex.submit(run_upload, t) for t in tasks]
            completed_count = 0
            for f in concurrent.futures.as_completed(futs):
                try:
                    exp_name, d = f.result()
                    completed_count += 1
                    print(f"Progress: {completed_count}/{len(tasks)} uploads processed", flush=True)
                    if d is not None:
                        uploaded[exp_name] = d
                except Exception as e:
                    print(f"Upload failed: {e}")
                    completed_count += 1

    # Build time-series relation with a single partition per time label
    if not uploaded:
        print("No averaged datasets uploaded; cannot create relation.")
        sys.exit(3)

    has_entity, _, part_key = compatibility_layer(res)
    partitions = []
    sorted_labels = sorted(uploaded.keys(), key=lambda nm: int(nm.replace("Expression ", "").replace(" Average", "").replace("hr", "")))
    for position, exp_name in enumerate(sorted_labels, start=1):
        d = uploaded[exp_name]
        lbl = exp_name.replace("Expression ", "").replace(" Average", "")
        parent = get_or_make_entity_for_data(d, lbl, None, new_collection)
        partitions.append({part_key: parent.id, "label": lbl, "position": position})

    # Create new time-series relation with unique slug
    try:
        relation = res.relation.create(
            type="series",
            category="Time series",
            collection=new_collection,
            slug=f"avg-relation-{new_collection.id}-{timestamp}",
            partitions=partitions,
        )
    except Exception as e:
        print(f"Slug 'avg-relation-{new_collection.id}-{timestamp}' may exist, retrying without explicit slug: {e}")
        relation = res.relation.create(
            type="series",
            category="Time series",
            collection=new_collection,
            partitions=partitions,
        )
    relation.save()
    print(f"Created averaged time-series relation '{relation.slug}' with {len(partitions)} partitions.")
    print(f"✅ Upload complete! Successfully uploaded {len(uploaded)} out of {len(tasks)} timepoints.")


if __name__ == "__main__":
    main()
