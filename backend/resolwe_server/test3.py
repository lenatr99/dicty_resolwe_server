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

inputs = {
    "h5ad_file": "AX4_2.h5ad.gz",
    "sc_name": "Test",
    "species": "Dictyostelium discoideum",
}

ref_seq = res.run(
    slug="upload-single-cell",
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
        