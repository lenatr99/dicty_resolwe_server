#!/usr/bin/env python3

import os
import sys
import logging

# Set environment variable FIRST before any imports
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")

# Now import the other modules
import resdk
import time
import asyncio
import django
# from resolwe_bio.tools import expression2storage

# Connect to Resolwe and login as admin
res = resdk.Resolwe(url="http://localhost:8000", username="admin", password="admin")
res.login(username="admin", password="admin")

test_collection = res.collection.create(name="Milestones collection", description="Milestones collection", tags=["community:bcm"])
print(
    f"✅ Successfully created collection: {test_collection.name} (ID: {test_collection.id})"
)

times = ["00hr", "04hr", "08hr", "12hr", "16hr", "20hr"]
reps = ["r1", "r2", "r3"]

exps = {}

for time_point in times:
    if time_point not in exps:
        exps[time_point] = []
    for rep in reps:
        inputs = {
            "rc": f"mutual_supressors_data/AX4_{time_point}_{rep}_rc.tab.gz",
            "exp": f"mutual_supressors_data/AX4_{time_point}_{rep}_polya.tab.gz",
            "exp_type": "polyA",
            "exp_name": "Expression",
            "source": "DICTYBASE",
            "species": "Dictyostelium discoideum",
            "build": "dd-05-2009",
        }

        ref_seq = res.run(
            slug="upload-expression",
            input=inputs,
            collection=test_collection,
        )

        # Wait until the first process is done
        print("⏳ Waiting for reference sequence upload to complete...")
        while ref_seq.status != "OK":
            if ref_seq.status == "ER":
                print(f"❌ Reference sequence process failed with status: {ref_seq.status}")
                print(
                    f"❌ Reference sequence process failed with error: {ref_seq.process_error}"
                )
                break
            print(f"📊 Current status: {ref_seq.status}")
            time.sleep(2)
            ref_seq.update()

        if ref_seq.status == "OK":
            print(f"✅ Reference sequence upload completed successfully!")
            exps[time_point].append(ref_seq)
        

# 0) CONFIG: list your Data objects per timepoint (replace with your actual Data objs)
# Example with one file you already uploaded, expand as you add uploads:
samples_by_time = {
    "hr00": exps["00hr"],             # e.g. 3 replicates -> [d00_r1, d00_r2, d00_r3]
    "hr04": exps["04hr"],
    "hr08": exps["08hr"],
    "hr12": exps["12hr"],
    "hr16": exps["16hr"],
    "hr20": exps["20hr"],
}

# 1) Compatibility layer: 'entity' vs 'sample'
HAS_ENTITY = hasattr(res, "entity")
RESOURCE = res.entity if HAS_ENTITY else res.sample
PART_KEY = "entity" if HAS_ENTITY else "sample"

def ensure_in_collection(obj, collection):
    # Handle both Entity (.collections) and Sample (.collection) objects
    if HAS_ENTITY:
        # Entity objects have .collections (plural)
        if collection.id not in [c.id for c in (obj.collections or [])]:
            obj.collections = (obj.collections or []) + [collection]
            obj.save()
    else:
        # Sample objects have .collection (singular)
        if hasattr(obj, 'collection') and obj.collection:
            if obj.collection.id != collection.id:
                obj.collection = collection
                obj.save()
        else:
            obj.collection = collection
            obj.save()

def get_or_make_entity_for_data(d, label, rep_idx, collection):
    # Prefer a preexisting entity/sample on the Data (processes often create it)
    parent = getattr(d, "entity", None) if HAS_ENTITY else getattr(d, "sample", None)
    if parent:
        ensure_in_collection(parent, collection)
        print(f"🧪 Using existing {PART_KEY} {parent.id} for Data {d.id}")
        return parent
    # Otherwise create a new one and attach the Data
    created = RESOURCE.create(
        name=f"{label}_r{rep_idx}",
        data=[d],
        collections=[collection],
    )
    print(f"🧪 Created {PART_KEY} {created.id} for Data {d.id}")
    return created

# 2) Sort labels and assign positions in time order
def label_to_hour(lbl):
    return int(lbl.lower().replace("hr", ""))

sorted_labels = sorted(samples_by_time.keys(), key=label_to_hour)
pos_by_label = {lbl: i + 1 for i, lbl in enumerate(sorted_labels)}

# 3) Build partitions (like your “AX4_2023” example)
partitions = []
for lbl in sorted_labels:
    for i, d in enumerate(samples_by_time[lbl], start=1):
        parent = get_or_make_entity_for_data(d, lbl, i, test_collection)
        partitions.append({
            "entity": parent.id,   # Relations API always expects 'entity' field
            "label": lbl,
            "position": pos_by_label[lbl],
        })

# 4) Optional descriptor/schema to mirror your example
descriptor = {
    "growth": "HL5",
    "strain": "AX4",
    "details": "AX4",
    "project": "12. Mutual suppression",
    "citation": {"url": "TO-DO", "name": "Katoh-Kurasawa M et. al."},
    "treatment": "Filter Development",
}
ds_qs = res.descriptor_schema.filter(slug="dicty-time-series")
descriptor_schema = ds_qs[0] if ds_qs else None

# 5) Create the relation (note lowercase category)
relation = res.relation.create(
    type="series",
    category="time series",
    collection=test_collection,
    slug="relation-60",  # change if it collides
    partitions=partitions,
)

relation.descriptor = descriptor
if descriptor_schema:
    relation.descriptor_schema = descriptor_schema
relation.save()

print(f"✅ Created relation '{relation.slug}' (ID: {relation.id}) with {len(partitions)} partitions.")
