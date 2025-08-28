#!/usr/bin/env python3
"""
Script to rename files from long GSM format to clean format.

Original format: GSM8583751_AX4_00hr_r1_1_S1_L001_R1_001_preprocessed_mapped_expression_rc.tab.gz
Clean format: AX4_00hr_r1_rc.tab.gz

Original format: GSM8583751_AX4_00hr_r1_1_S1_L001_R1_001_preprocessed_mapped_expression_rpkum_polya.tab.gz
Clean format: AX4_00hr_r1_polya.tab.gz
"""

import os
import re
import sys

def clean_filename(filename):
    """
    Convert long GSM filename to clean format.
    
    Args:
        filename (str): Original filename
        
    Returns:
        str: Clean filename or None if pattern doesn't match
    """
    # Pattern to match the original filename format
    pattern = r'GSM\d+_(AX4_\d+hr)_r(\d+)_\d+_S\d+_L001_R1_001_preprocessed_mapped_expression_(rc|rpkum_polya)\.tab\.gz'
    
    match = re.match(pattern, filename)
    if not match:
        return None
    
    time_strain = match.group(1)  # e.g., "AX4_00hr"
    replicate = match.group(2)    # e.g., "1", "3", "4"
    exp_type = match.group(3)     # "rc" or "rpkum_polya"
    
    # Convert rpkum_polya to polya
    if exp_type == "rpkum_polya":
        exp_type = "polya"
    
    # Create clean filename
    clean_name = f"{time_strain}_r{replicate}_{exp_type}.tab.gz"
    
    return clean_name

def rename_files_in_directory():
    """
    Rename all matching files in the current directory.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Get all .gz files in the directory
    files = [f for f in os.listdir(script_dir) if f.endswith('.tab.gz')]
    
    renamed_count = 0
    skipped_count = 0
    
    print(f"Found {len(files)} .tab.gz files to process...")
    print()
    
    for filename in files:
        clean_name = clean_filename(filename)
        
        if clean_name is None:
            print(f"SKIP: {filename} (doesn't match expected pattern)")
            skipped_count += 1
            continue
        
        if clean_name == filename:
            print(f"SKIP: {filename} (already has clean name)")
            skipped_count += 1
            continue
        
        old_path = os.path.join(script_dir, filename)
        new_path = os.path.join(script_dir, clean_name)
        
        # Check if target file already exists
        if os.path.exists(new_path):
            print(f"ERROR: {filename} -> {clean_name} (target already exists)")
            skipped_count += 1
            continue
        
        try:
            os.rename(old_path, new_path)
            print(f"RENAMED: {filename} -> {clean_name}")
            renamed_count += 1
        except OSError as e:
            print(f"ERROR: Failed to rename {filename}: {e}")
            skipped_count += 1
    
    print()
    print(f"Summary: {renamed_count} files renamed, {skipped_count} files skipped")

def preview_renames():
    """
    Preview what renames would be performed without actually renaming.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Get all .gz files in the directory
    files = [f for f in os.listdir(script_dir) if f.endswith('.tab.gz')]
    
    print(f"Found {len(files)} .tab.gz files to process...")
    print()
    print("Preview of renames (no files will be changed):")
    print("-" * 80)
    
    for filename in files:
        clean_name = clean_filename(filename)
        
        if clean_name is None:
            print(f"SKIP: {filename} (doesn't match expected pattern)")
            continue
        
        if clean_name == filename:
            print(f"SKIP: {filename} (already has clean name)")
            continue
        
        print(f"{filename}")
        print(f"  -> {clean_name}")
        print()

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--preview":
        preview_renames()
    else:
        print("File Renaming Script")
        print("=" * 50)
        print()
        
        response = input("Do you want to preview the renames first? (y/n): ").lower().strip()
        if response in ['y', 'yes']:
            preview_renames()
            print()
            response = input("Proceed with renaming? (y/n): ").lower().strip()
            if response not in ['y', 'yes']:
                print("Renaming cancelled.")
                sys.exit(0)
        
        rename_files_in_directory()
