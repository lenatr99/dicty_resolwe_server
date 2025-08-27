#!/usr/bin/env python3
"""
Complete annotation setup script based on resolwe-bio test patterns.
This creates all the annotation fields that the bioinformatics processes expect.
"""

import os
import sys
import django

# Set up Django environment
sys.path.insert(0, '/Users/lenatrnovec/backend')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "resolwe_server.settings")
django.setup()

from django.apps import apps
from resolwe.flow.models.annotations import AnnotationField, AnnotationGroup

def setup_complete_annotations():
    """Set up complete annotation schema like in BioProcessTestCase."""
    
    print("Setting up complete annotation schema based on resolwe-bio test patterns...")
    
    try:
        # Check what already exists
        print("🔍 Checking existing annotation schema...")
        existing_groups = AnnotationGroup.objects.all()
        existing_fields = AnnotationField.objects.all()
        
        print(f"  Found {existing_groups.count()} groups and {existing_fields.count()} fields")
        for group in existing_groups:
            print(f"    Group: {group.name}")
        for field in existing_fields:
            print(f"    Field: {field.group.name}.{field.name} (type: {field.type})")
        
        # Don't delete existing - just add missing ones
        
        # Create or get general annotation group
        general_group, created = AnnotationGroup.objects.get_or_create(
            name="general",
            defaults={"sort_order": 1, "label": "General Information"}
        )
        if created:
            print("✅ Created general annotation group")
        else:
            print("🔍 General annotation group already exists")
        
        species_field, created = AnnotationField.objects.get_or_create(
            name="species",
            group=general_group,
            defaults={
                "sort_order": 1,
                "type": "STRING",  # Using the constraint-approved type
                "label": "Species",
                "description": "Species information",
                "required": False,
                "vocabulary": {
                    "Caenorhabditis elegans": "Caenorhabditis elegans",
                    "Cricetulus griseus": "Cricetulus griseus",
                    "Dictyostelium discoideum": "Dictyostelium discoideum",
                    "Dictyostelium purpureum": "Dictyostelium purpureum",
                    "Drosophila melanogaster": "Drosophila melanogaster",
                    "Homo sapiens": "Homo sapiens",
                    "Macaca mulatta": "Macaca mulatta",
                    "Mus musculus": "Mus musculus",
                    "Rattus norvegicus": "Rattus norvegicus",
                    "other": "Other",
                },
            }
        )
        if created:
            print("✅ Created general.species field")
        else:
            print("🔍 General.species field already exists")
        
        # Add other general fields that might be useful
        print("✅ Creating additional general fields...")
        AnnotationField.objects.create(
            name="description",
            sort_order=2,
            group=general_group,
            type="STRING",
            label="Description",
            description="General description",
            required=False,
        )
        
        AnnotationField.objects.create(
            name="biosample_source",
            sort_order=3,
            group=general_group,
            type="STRING",
            label="Biosample Source",
            description="Source of the biosample",
            required=False,
        )
        
        AnnotationField.objects.create(
            name="biosample_treatment",
            sort_order=4,
            group=general_group,
            type="STRING",
            label="Biosample Treatment",
            description="Treatment applied to the biosample",
            required=False,
        )
        
        # Create biospecimen_information group and fields
        print("✅ Creating biospecimen_information annotation group...")
        biospecimen_group = AnnotationGroup.objects.create(
            name="biospecimen_information", 
            sort_order=2,
            label="Biospecimen Information"
        )
        
        AnnotationField.objects.create(
            name="organ",
            sort_order=1,
            group=biospecimen_group,
            type="STRING",
            label="Organ",
            description="Organ information",
            required=False,
        )
        
        # Create QC group and fields
        print("✅ Creating qc annotation group...")
        qc_group = AnnotationGroup.objects.create(name="qc", sort_order=3, label="Quality Control")
        
        AnnotationField.objects.create(
            name="status",
            sort_order=1,
            group=qc_group,
            type="STRING",
            label="QC Status",
            description="Quality control status",
            required=False,
        )
        
        AnnotationField.objects.create(
            name="message",
            sort_order=2,
            group=qc_group,
            type="STRING",
            label="QC Message",
            description="Quality control message",
            required=False,
        )
        
        print("\n🎉 Complete annotation schema setup successful!")
        
        # Verify what was created
        print("📋 Summary of created annotation schema:")
        for group in AnnotationGroup.objects.all().order_by('sort_order'):
            print(f"  📁 Group: {group.name} ({group.label})")
            for field in AnnotationField.objects.filter(group=group).order_by('sort_order'):
                print(f"    📄 Field: {group.name}.{field.name} ({field.label}) - type: {field.type}")
                if field.vocabulary:
                    print(f"        🔤 Vocabulary: {list(field.vocabulary.keys())}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error setting up annotation schema: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    try:
        if setup_complete_annotations():
            print("\n✅ Annotation schema is now ready for bioinformatics processes!")
            print("The HISAT2 alignment process should now work correctly.")
        else:
            print("\n❌ Failed to set up annotation schema.")
            sys.exit(1)
    except Exception as e:
        print(f"❌ Fatal error: {e}")
        sys.exit(1)
