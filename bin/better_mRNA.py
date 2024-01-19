#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def add_mRNA_feature(input_gff_path, output_gff_path = "temp.gff"):


    def extract_attribute(attribute_string, attribute_name):
        # Using regex to find the value for the specified attribute, case-insensitive
        match = re.search(f"{attribute_name}=([^;]+)", attribute_string)
        return match.group(1) if match else None


    """
    This function takes a GFF file as input and adds mRNA features where they are missing.
    It also adjusts the parentage of CDS features to match the new mRNA features.
    At the end, each gene should have a corresponding mRNA feature and each mRNA should have its CDSs as children.

    """

    with open(input_gff_path, 'r') as file:
        lines = file.readlines()


    ##############################
    # First pass to check for mRNA
    # in each gene and store in a
    # dictionary
    ##############################
        
    modified_lines = []
    gene_to_mrna = {}
    for line in lines:
        if line.strip() and not line.startswith("#"):
            parts = line.split("\t")

            # spec from https://m.ensembl.org/info/website/upload/gff3.html
            feature_type = parts[2].strip()
            attributes = parts[8].strip()

            if feature_type == "gene":

                # Hard coded because uses output of gffread -C command
                gene_id = extract_attribute(attributes, "ID")
                if gene_id:
                    gene_to_mrna[gene_id] = None  # Initially no mRNA for this gene
                else:
                    print("Gene ID not found in GFF")
                    warn(f"Gene ID not found in line: {line}")
                    pass

            elif feature_type == "mRNA":

                parent_gene_id = extract_attribute(attributes, "Parent")
                mRNA_ID = extract_attribute(attributes, "ID")
                if parent_gene_id and parent_gene_id in gene_to_mrna:
                    gene_to_mrna[parent_gene_id] = mRNA_ID 

                else:
                    pass
                    
    ##############################
    # Second pass to add mRNA and
    # adjust CDS parentage
    # Not optimal but heh
    ##############################
    for line in lines:
        if line.strip() and not line.startswith("#"):
            parts = line.split("\t")
            feature_type = parts[2].strip()
            attributes = parts[8].strip()

            if feature_type == "gene":

                gene_id = extract_attribute(attributes, "ID")

                # If mRNA doesn't exist, create it
                if gene_id and not gene_to_mrna[gene_id]:  

                    # mRNA copies the gene line, but with mRNA as the feature type and a new ID
                    mrna_line = "\t".join(parts[:2] + ["mRNA"] + parts[3:8] + [f"ID={gene_id}_mRNA;Parent={gene_id}\n"])
                    modified_lines.append(mrna_line)
                modified_lines.append(line)

            elif feature_type == "CDS":
                parent_gene_id = extract_attribute(attributes, "Parent")
                # Check if the parent gene ID exists, if not, skip or handle appropriately
                if parent_gene_id in gene_to_mrna:
                    if not gene_to_mrna[parent_gene_id]:  # If mRNA was added, adjust CDS parent
                        attributes = f"Parent={parent_gene_id}_mRNA\n"
                    parts[8] = attributes
                cds_line = "\t".join(parts)
                modified_lines.append(cds_line)
            else:
                # Handle other feature types or missing parent gene case
                modified_lines.append(line)

        else:
            modified_lines.append(line)  # Add comments and empty lines as they are

        
    with open(output_gff_path, 'w') as file:
        file.writelines(modified_lines)

    #os.remove("coding.gff")
    return 0


if __name__ == "__main__":

    import argparse
    import os 
    import re
    from warnings import warn

    parser = argparse.ArgumentParser(description="Add mRNA features to GFF file")
    parser.add_argument("-i", "--input_gff", help="Path to input GFF file")
    parser.add_argument("-o", "--output_gff", default="temp.gff", help="Path to output GFF file")
    args = parser.parse_args()

    add_mRNA_feature(args.input_gff, args.output_gff)


