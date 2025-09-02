#!/usr/bin/env python3

import os
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import resdk
import gzip
import concurrent.futures

# Set environment variable FIRST before any imports
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")



# Connect to Resolwe and setup
res = resdk.Resolwe(url="http://localhost:8000", username="admin", password="admin")
res.login(username="admin", password="admin")

test_collection = res.collection.create(
    name="Single cell collection",
    description="Single cell collection",
    tags=["community:bcm"],
)
print(f"✅ Created collection: {test_collection.name} (ID: {test_collection.id})")

adata = sc.read_h5ad("AX4.h5ad")
WORKERS = 10
GENES_NUM = adata.shape[1]
GENES_NUM = 10

#make directory for temp_files
os.makedirs("data/upload", exist_ok=True)

# Compatibility layer
HAS_ENTITY = hasattr(res, "entity")
RESOURCE = res.entity if HAS_ENTITY else res.sample
PART_KEY = "entity" if HAS_ENTITY else "sample"


def ensure_in_collection(obj, collection):
    if HAS_ENTITY:
        if collection.id not in [c.id for c in (obj.collections or [])]:
            obj.collections = (obj.collections or []) + [collection]
            obj.save()
    else:
        if hasattr(obj, "collection") and obj.collection:
            if obj.collection.id != collection.id:
                obj.collection = collection
                obj.save()
        else:
            obj.collection = collection
            obj.save()


def get_or_make_entity_for_data(d, label, rep_idx, collection):
    
    parent = getattr(d, "entity", None) if HAS_ENTITY else getattr(d, "sample", None)
    if parent:
        ensure_in_collection(parent, collection)
        print(f"🧪 Using existing {PART_KEY} {parent.id} for Data {d.id}")
        return parent
    if rep_idx is None:
        created = RESOURCE.create(
            name=label, data=[d], collections=[collection]
        )
    else:
        created = RESOURCE.create(
            name=f"{label}_r{rep_idx}", data=[d], collections=[collection]
        )
    print(f"🧪 Created {PART_KEY} {created.id} for Data {d.id}")
    return created


def upload_expression(args):
    """Unified upload function for both bulk and single gene data"""
    file_path, exp_name, ind, total = args
    start_time = time.time()

    inputs = {
        "exp": file_path,
        "exp_type": "10x_UMI",
        "exp_name": exp_name,
        "source": "DICTYBASE",
        "species": "Dictyostelium discoideum",
        "build": "dd-05-2009",
    }

    ref_seq = res.run(
        slug="upload-expression", input=inputs, collection=test_collection
    )
    print(f"⏳ Processing {ind+1}/{total}: {exp_name}")

    while ref_seq.status not in ["OK", "ER"]:
        time.sleep(2)
        ref_seq.update()

    if ref_seq.status == "ER":
        print(f"❌ Upload {ind+1} failed: {ref_seq.process_error}")
        return None

    elapsed = time.time() - start_time
    print(f"✅ Upload {ind+1}/{total} completed in {elapsed:.2f}s")
    return ref_seq, exp_name


# 1. Generate bulk expression files for time points and markers
bulk_tasks = []
times = ["00hr", "04hr", "08hr", "12hr", "16hr", "20hr"]
reps = adata.obs["marker"].unique()
time_exps = {t: [] for t in times}

for t in times:
    A = adata[adata.obs["time"] == t]
    for r in reps:
        B = A[A.obs["marker"] == r]
        if B.n_obs == 0:
            continue
        X = B.layers.get(None, B.X)
        m = (
            np.asarray(X.mean(axis=0)).ravel()
            if sp.issparse(X)
            else np.asarray(X.mean(axis=0))
        )
        out = f"data/upload/AX4_{t}_{r}.tab.gz"
        pd.DataFrame({"Gene": B.var["gene_ids"], "Expression": m}).to_csv(
            out, sep="\t", index=False, compression="gzip"
        )
        bulk_tasks.append(
            (out, f"Expression {t} {r}", len(bulk_tasks), len(times) * len(reps))
        )

# 2. Generate single gene files
adata_subset = adata[:, :GENES_NUM].copy()
gene_cache = {gene_id: idx for idx, gene_id in enumerate(adata_subset.var["gene_ids"])}

for ind, gene in enumerate(adata_subset.var["gene_ids"]):
    i = gene_cache[gene]
    X = adata_subset.X

    if sp.isspmatrix_csr(X):
        values = X[:, i].toarray().ravel()
    elif sp.isspmatrix_csc(X):
        values = X[:, i].toarray().ravel()
    else:
        values = np.asarray(X[:, i]).ravel()

    lines = ["Gene\tExpression"]
    lines.extend(
        f"{adata_subset.obs_names[idx]}\t{value}" for idx, value in enumerate(values)
    )

    out_path = f"data/upload/temp_{gene}.tab.gz"
    with gzip.open(out_path, "wt", compresslevel=1, encoding="utf-8") as f:
        f.write("\n".join(lines))

    bulk_tasks.append(
        (out_path, f"Single cell {gene}", len(bulk_tasks), len(bulk_tasks) + GENES_NUM)
    )

for i in range(2):
    values = adata_subset.obsm["X_uce_umap"][:, i]
    print(len(values))

    lines = ["Gene\tExpression"]
    lines.extend(
        f"{adata_subset.obs_names[idx]}\t{value}" for idx, value in enumerate(values)
    )

    out_path = f"data/upload/temp_X_uce_umap_{i}.tab.gz"
    with gzip.open(out_path, "wt", compresslevel=1, encoding="utf-8") as f:
        f.write("\n".join(lines))

    bulk_tasks.append(
        (out_path, f"X_uce_umap_{i}", len(bulk_tasks), len(bulk_tasks) + GENES_NUM)
    )


#make it to times by taking first 2 characters
adata_subset.obs["time"] = adata_subset.obs["time"].str[:2].astype(int)
values = adata_subset.obs["time"]
print(len(values))

lines = ["Gene\tExpression"]
lines.extend(
    f"{adata_subset.obs_names[idx]}\t{value}" for idx, value in enumerate(values)
)

out_path = f"data/upload/temp_time.tab.gz"
with gzip.open(out_path, "wt", compresslevel=1, encoding="utf-8") as f:
    f.write("\n".join(lines))

bulk_tasks.append(
    (out_path, f"Single cell time", len(bulk_tasks), len(bulk_tasks) + GENES_NUM)
)

# 3. Process all uploads in parallel batches
all_results = {}
batch_size = WORKERS
total_batches = (len(bulk_tasks) + batch_size - 1) // batch_size

for batch_num in range(total_batches):
    batch_start = batch_num * batch_size
    batch_end = min(batch_start + batch_size, len(bulk_tasks))
    current_batch = bulk_tasks[batch_start:batch_end]

    print(
        f"📦 Processing batch {batch_num + 1}/{total_batches} with {len(current_batch)} uploads"
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as executor:
        futures = [executor.submit(upload_expression, task) for task in current_batch]

        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                if result:
                    ref_seq, name = result
                    all_results[name] = ref_seq

                    # Categorize results
                    if "Single cell" in name:
                        gene = name.replace("Single cell ", "")
                        # Will be used for single gene partitions
                    else:
                        # Extract time from name like "Expression 00hr mCerulean"
                        parts = name.split()
                        if len(parts) >= 2:
                            time_key = parts[1]  # "00hr"
                            time_exps[time_key].append(ref_seq)
            except Exception as exc:
                print(f"❌ Upload generated exception: {exc}")

print(f"✅ Successfully uploaded {len(all_results)} files")

# 4. Create relations and partitions
time_partitions = []

# Time series partitions
samples_by_time = {f"hr{t[:2]}": time_exps[t] for t in times}
sorted_time_labels = sorted(samples_by_time.keys(), key=lambda x: int(x[2:]))
time_pos = {lbl: i + 1 for i, lbl in enumerate(sorted_time_labels)}

for lbl in sorted_time_labels:
    for i, d in enumerate(samples_by_time[lbl], start=1):
        parent = get_or_make_entity_for_data(d, lbl, i, test_collection)
        time_partitions.append(
            {"entity": parent.id, "label": lbl, "position": time_pos[lbl]}
        )

# Make entities for X_uce_umap
umap_results = {k: v for k, v in all_results.items() if "X_uce_umap" in k}
for i, (name, data_obj) in enumerate(umap_results.items()):
    parent = get_or_make_entity_for_data(data_obj, name, None, test_collection)
    time_partitions.append(
        {"entity": parent.id, "label": name, "position": len(sorted_time_labels) + i + 1}
    )

umap_results = {k: v for k, v in all_results.items() if "time" in k}
for i, (name, data_obj) in enumerate(umap_results.items()):
    parent = get_or_make_entity_for_data(data_obj, "time", i, test_collection)
    time_partitions.append(
        {"entity": parent.id, "label": f"time", "position": len(sorted_time_labels) + i + 1}
    )

gene_partitions = []
# Single gene partitions
gene_results = {
    k.replace("Single cell ", ""): v
    for k, v in all_results.items()
    if "Single cell" in k
}
sorted_genes = sorted(gene_results.keys())
gene_pos = {gene: i + 1 for i, gene in enumerate(sorted_genes)}

for gene, data_obj in gene_results.items():
    parent = get_or_make_entity_for_data(data_obj, gene, 1, test_collection)
    gene_partitions.append({"entity": parent.id, "label": gene, "position": gene_pos[gene]})

# 5. Create relation
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


ID_NUMBER = 5008
relation = res.relation.create(
    type="series",
    category="Time series",
    collection=test_collection,
    slug=f"relation-{ID_NUMBER}",
    partitions=time_partitions,
)

relation.descriptor = descriptor
if descriptor_schema:
    relation.descriptor_schema = descriptor_schema
relation.save()

ds_qs = res.descriptor_schema.filter(slug="dicty-time-series")
descriptor_schema = ds_qs[0] if ds_qs else None

relation2 = res.relation.create(
    type="series",
    category="Single cell series",
    collection=test_collection,
    slug=f"relation-{ID_NUMBER}_sc",
    partitions=gene_partitions,
)

relation2.descriptor = descriptor
if descriptor_schema:
    relation2.descriptor_schema = descriptor_schema
relation2.save()

print(
    f"✅ Created relation '{relation.slug}' (ID: {relation.id}) with {len(time_partitions)} time partitions and relation '{relation2.slug}' (ID: {relation2.id}) with {len(gene_partitions)} gene partitions."
)
