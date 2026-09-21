#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
from warnings import warn

def extract_attribute(attribute_string, attribute_name):
    """
    Extract the value for attribute_name safely, using regex anchor to avoid
    matching substrings like copy_num_ID= when searching for ID=.
    """
    match = re.search(r"(?:^|;)\s*" + re.escape(attribute_name) + r"=([^;\s]+)", attribute_string)
    return match.group(1) if match else None

def add_mRNA_feature(input_gff_path, output_gff_path="temp.gff"):
    """
    This function takes a GFF file as input and adds mRNA features where they are missing.
    It supports multi-isoform CDS records under a single gene by grouping CDS records by
    transcript/protein identifier (ID, protein_id, transcript_id) to avoid merging distinct
    isoforms into a single mRNA (which causes memory corruption in gffread).
    It also preserves all original attributes of CDS features while adjusting Parentage.
    """
    with open(input_gff_path, 'r') as file:
        lines = file.readlines()

    # Pass 1: Collect pre-existing mRNA / transcript IDs
    existing_mrna_ids = set()
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        feature_type = parts[2].strip()
        attributes = parts[8].strip()

        if feature_type in ["mRNA", "transcript"]:
            mRNA_ID = extract_attribute(attributes, "ID")
            if mRNA_ID:
                existing_mrna_ids.add(mRNA_ID)

    # Pass 2: Identify CDS features needing synthetic mRNA parent
    cds_to_new_parent = {}
    synthetic_mrna_defs = {}  # key: (parent_gene_id, mrna_id) -> metadata dict

    for idx, line in enumerate(lines):
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        feature_type = parts[2].strip()
        attributes = parts[8].strip()

        if feature_type == "CDS":
            parent = extract_attribute(attributes, "Parent")
            if parent and parent not in existing_mrna_ids:
                # Parent points to a gene or feature that is not an existing mRNA
                cds_id = extract_attribute(attributes, "ID")
                prot_id = extract_attribute(attributes, "protein_id")
                tx_id = extract_attribute(attributes, "transcript_id")

                sub_id = cds_id or prot_id or tx_id
                if sub_id:
                    mrna_id = f"{sub_id}_mRNA" if not sub_id.endswith("_mRNA") else sub_id
                else:
                    mrna_id = f"{parent}_mRNA"

                cds_to_new_parent[idx] = mrna_id

                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    start, end = 0, 0
                contig, source, strand = parts[0], parts[1], parts[6]

                key = (parent, mrna_id)
                if key not in synthetic_mrna_defs:
                    synthetic_mrna_defs[key] = {
                        'contig': contig, 'source': source,
                        'start': start, 'end': end,
                        'strand': strand
                    }
                else:
                    synthetic_mrna_defs[key]['start'] = min(synthetic_mrna_defs[key]['start'], start)
                    synthetic_mrna_defs[key]['end'] = max(synthetic_mrna_defs[key]['end'], end)

    # Pass 3: Rebuild GFF with synthetic mRNA features and updated CDS parents
    modified_lines = []
    written_synthetic_mrnas = set()

    for idx, line in enumerate(lines):
        if not line.strip() or line.startswith("#"):
            modified_lines.append(line)
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            modified_lines.append(line)
            continue

        feature_type = parts[2].strip()
        attributes = parts[8].strip()

        if feature_type == "gene":
            modified_lines.append(line)
            gene_id = extract_attribute(attributes, "ID")
            if gene_id:
                # Output any synthetic mRNA defined for this gene
                for (parent_gid, mrna_id), mdata in synthetic_mrna_defs.items():
                    if parent_gid == gene_id and mrna_id not in written_synthetic_mrnas:
                        mrna_line = f"{mdata['contig']}\t{mdata['source']}\tmRNA\t{mdata['start']}\t{mdata['end']}\t.\t{mdata['strand']}\t.\tID={mrna_id};Parent={gene_id}\n"
                        modified_lines.append(mrna_line)
                        written_synthetic_mrnas.add(mrna_id)

        elif feature_type == "CDS" and idx in cds_to_new_parent:
            new_parent = cds_to_new_parent[idx]
            if "Parent=" in attributes:
                new_attributes = re.sub(r"Parent=[^;\s]+", f"Parent={new_parent}", attributes)
            else:
                new_attributes = f"Parent={new_parent};" + attributes
            parts[8] = new_attributes
            modified_lines.append("\t".join(parts) + "\n")
        else:
            modified_lines.append(line)

    with open(output_gff_path, 'w') as file:
        file.writelines(modified_lines)

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add mRNA features to GFF file")
    parser.add_argument("-i", "--input_gff", required=True, help="Path to input GFF file")
    parser.add_argument("-o", "--output_gff", default="temp.gff", help="Path to output GFF file")
    args = parser.parse_args()

    add_mRNA_feature(args.input_gff, args.output_gff)
