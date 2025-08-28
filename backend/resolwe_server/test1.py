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

# Connect to Resolwe and login as admin
res = resdk.Resolwe(url="http://localhost:8000", username="admin", password="admin")
res.login(username="admin", password="admin")

test_collection = res.collection.create(name="Test collection")
print(
    f"✅ Successfully created collection: {test_collection.name} (ID: {test_collection.id})"
)

inputs = {
    "src": "./genome.fasta.gz",
    "species": "Dictyostelium discoideum",
    "build": "dd-05-2009",
}

ref_seq = res.run(
    slug="upload-fasta-nucl",
    input=inputs,
    collection=test_collection,
)

print(f"✅ Successfully created reference sequence: {ref_seq.id}")
print(f"📊 Initial status: {ref_seq.status}")

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

    data = {
        "ref_seq": ref_seq.id,
    }

    hisat2_index = res.run(
        slug="hisat2-index",
        input=data,
        collection=test_collection,
    )

    print(f"✅ Successfully created hisat2 index process: {hisat2_index.id}")
    print(f"📊 Initial status: {hisat2_index.status}")

    # Wait until the second process is done
    print("⏳ Waiting for hisat2 index to complete...")
    while hisat2_index.status != "OK":
        if hisat2_index.status == "ER":
            print(f"❌ Hisat2 index process failed with status: {hisat2_index.status}")
            break
        print(f"📊 Current status: {hisat2_index.status}")
        time.sleep(2)
        hisat2_index.update()

    if hisat2_index.status == "OK":
        print(f"✅ Successfully created hisat2 index: {hisat2_index.output}")
    else:
        print(f"❌ Hisat2 index failed with status: {hisat2_index.status}")
        print(f"❌ Hisat2 index failed with error: {hisat2_index.process_error}")
else:
    print(
        f"❌ Cannot proceed to hisat2 index - reference sequence failed with status: {ref_seq.status}"
    )


reads = res.run(
    slug="upload-fastq-single",
    input={
        "src": "./reads.fastq.gz",
    },
    collection=test_collection,
)

print(f"✅ Successfully created reads upload process: {reads.id}")
print(f"📊 Initial status: {reads.status}")

# Wait until the third process is done
print("⏳ Waiting for reads upload to complete...")
while reads.status != "OK":
    if reads.status == "ER":
        print(f"❌ Reads upload process failed with status: {reads.status}")
        break
    print(f"📊 Current status: {reads.status}")
    time.sleep(2)
    reads.update()

if reads.status == "OK":
    print(f"✅ Reads upload completed successfully!")


inputs = {
    "src": "annotation_dicty.gff.gz",
    "source": "DICTYBASE",
    "species": "Dictyostelium discoideum",
    "build": "dd-05-2009",
}

annotation = res.run(
    slug="upload-gff3",
    input=inputs,
    collection=test_collection,
)
print(f"✅ Successfully created annotation process: {annotation.id}")


while annotation.status != "OK":
    if annotation.status == "ER":
        print(f"❌ Annotation process failed with status: {annotation.status}")
        break
    print(f"📊 Current status: {annotation.status}")
    time.sleep(2)
    annotation.update()

if annotation.status == "OK":
    print(f"✅ Successfully created annotation: {annotation.output}")
else:
    print(f"❌ Annotation failed with status: {annotation.status}")
    print(f"❌ Annotation failed with error: {annotation.process_error}")


aligned_reads = res.run(
    slug="alignment-hisat2",
    input={"genome": hisat2_index, "reads": reads},
    collection=test_collection,
)

while aligned_reads.status != "OK":
    if aligned_reads.status == "ER":
        print(f"❌ Alignment process failed with status: {aligned_reads.status}, error: {aligned_reads.process_error}")
        break
    print(f"📊 Current status: {aligned_reads.status}")
    time.sleep(2)
    aligned_reads.update()

if aligned_reads.status == "OK":
    print(f"✅ Successfully created alignment process: {aligned_reads.id}")

mappa = res.run(
    slug="upload-mappability",
    input={"src": "purpureum_mappability_50.tab.gz"},
    collection=test_collection,
)
print(f"✅ Successfully created mappability process: {mappa.id}")

while mappa.status != "OK":
    if mappa.status == "ER":
        print(f"❌ Mappability process failed with status: {mappa.status}")
        break
    print(f"📊 Current status: {mappa.status}")
    time.sleep(2)
    mappa.update()

if mappa.status == "OK":
    print(f"✅ Successfully created mappability: {mappa.output}")
else:
    print(f"❌ Mappability failed with status: {mappa.status}")
    print(f"❌ Mappability failed with error: {mappa.process_error}")



inputs = {
    "alignment": aligned_reads.id,
    "gff": annotation.id,
    "mappable": mappa.id,
}
expression = res.run(
    slug="expression-dicty",
    input=inputs,
    collection=test_collection,
)
print(f"✅ Successfully created expression process: {expression.id}")
while expression.status != "OK":
    if expression.status == "ER":
        print(f"❌ Expression process failed with status: {expression.status}")
        break
    print(f"📊 Current status: {expression.status}")
    time.sleep(2)
    expression.update()

if expression.status == "OK":
    print(f"✅ Successfully created expression: {expression.output}")