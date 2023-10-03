#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
--Check Synteny--


Author : Paul Roginski
mail : plm.roginski@gmail.com
version : sept-2023


This program tests if a given object from genome A is in micro-synteny with a
given object from genome B.

In other words, it looks for genes around a given position in genome A and try
to find their ortholog around the given position in genome B.

It can perform a synteny check between any pair objects A-B. Object A or
B can be :
    - any genomic region defined with a sequence (often a chromosome), a start
    (0-based), and a stop (included),
    - any genomic features with an 'ID' attribute in the GFF3 file.


Examples of questions it can answer :

    - "Is the gene-HBB in Human in synteny with gene-HBB in Chimpanzee ?"

    - "Is transcript:Os04t0508150-01 in Oryza sativa in synteny with the region
    20605296 to 20604985 on the fourth chromosome of Oryza meridionalis ?


It takes as input :
    - the GFF3 annotation file for genome A
    - the GFF3 annotation file for genome B
    - a .tsv file with object from genome A and object from genome B
    (see --help for details)
    - a .tsv file with pairs of orthologous genes bewteen the two genomes :
        - column 1 : gene from genome A
        - column 2 : its orthologous gene from genome B


A few reminders about GFF3 (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) :
    - col 4 (start) <= col 5 (end)
    - col 9 (attributes) : 'Multiple tag=value pairs are separated by semicolons'

Orthologous genes from the tsv must be an 'ID' value in the corresponding GFF3
file, ortherwise they won't be mapped.
"""


import configargparse
import yaml
from pybedtools.bedtool import BedTool
import re
import pandas as pd
import sys
import string
import numpy as np
from dcjoin import calculate_distance
import csv
import prettytable as pt
import math
import os


# Define how to represent :
EMPTY_CHAR = " "  # No gene
ANY_CHAR = "."  # A gene that has not ortholog.
DEV_MODE = False  # THIS IS FOR THE DEVELOPPER
HOME = "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/"  # YOU CAN REMOVE IT


def micro_syteny_mapping_interpreter(mapping,
                                     ins_max,
                                     del_max,
                                     rea_max,
                                     up_min,
                                     ov_min,
                                     do_min):
    """
    Provides a list of statistics from a microsynteny mapping dictionnary.

    Counts the number of :
        - insertions (e.i. genes that "are in B and not in A"),
        - deletions (e.i. genes from A that are not found (their prtholog) in
        B),
        - rearrangments (double DNA cut and join, minimal number of movments),
        - ortho_upstream/overlap/downstream (number of genes from A that are
        found (their homolog), in a specific position with respect to the
        reference interval in B).
    """

    A_letters = mapping['A ortholog']
    B_letters = mapping['B ortholog']

    A_letters_formated = [
        ele for ele in A_letters if ele not in [
            ANY_CHAR, EMPTY_CHAR]]
    B_letters_formated = []

    nb_insertions = 0

    for ele in B_letters:

        if ele not in [ANY_CHAR, EMPTY_CHAR]:

            if ele[-3:] == '(r)':

                B_letters_formated.append(ele[:-3].lower())

            else:

                B_letters_formated.append(ele)

        elif ele == ANY_CHAR:

            nb_insertions += 1

    nb_deletions = len(A_letters_formated) - len(B_letters_formated)

    nb_rearrangements = None
    try:
        # The number of rearrangment events is computed using a customized version of
        # https://github.com/mlliou112/py-dcj, which implements the algorithm proposed
        # in Braga et al. (2010). For our purpose, it is mainly based in the reference
        # paper Yancopoulos et al. (2005).
        # Consider only the order and the orientation (straight or upside-down) of
        # the letters in the analysis vector.

        # Flanking "." represent the ends of a linear sequence.
        # "y" and "z" fix the strand orientation
        nb_rearrangements = int(calculate_distance(
            ['.y' + ''.join(A_letters_formated) + 'z.'],
            ['.y' + ''.join(B_letters_formated) + 'z.'], method="dcj"
        )
        )

    except ValueError:
        print("ValueError: Duplicated markers are not allowed.")
        nb_rearrangements = -1
        
    # Count the number of orthologs found in B for each 'side' (upstream,
    # downstream, and overlapping).

    # The index of genes in B that have their ortholog in A.
    B_ortho_index = [i for i in range(len(
        mapping['A position'])) if mapping['B ortholog'][i] not in [EMPTY_CHAR, ANY_CHAR]]

    nb_ortho_up = len(
        [i for i in B_ortho_index if 'upstream' in mapping['B position'][i]])
    nb_ortho_ov = len(
        [i for i in B_ortho_index if 'overlaps' in mapping['B position'][i]])
    nb_ortho_do = len(
        [i for i in B_ortho_index if 'downstream' in mapping['B position'][i]])

    # Interpretation
    print("nb_ins {}".format(nb_insertions))
    print("nb_del {}".format(nb_deletions))
    print("nb_rea {}".format(nb_rearrangements))
    print("nb_up {}".format(nb_ortho_up))
    print("nb_ov {}".format(nb_ortho_ov))
    print("nb_do {}".format(nb_ortho_do))

    if nb_insertions <= ins_max and nb_deletions <= del_max and nb_rearrangements <= rea_max and nb_ortho_up >= up_min and nb_ortho_ov >= ov_min and nb_ortho_do >= do_min:

        isSyntenic = True

    else:

        isSyntenic = False

    return ([nb_insertions, nb_deletions, nb_rearrangements,
            nb_ortho_up, nb_ortho_ov, nb_ortho_do, isSyntenic])


def micro_synteny_mapper(
        A_objs,
        A_obj,
        B_objs,
        B_obj,
        orthologs_df,
        A_strands,
        B_strands,
        ins_max,
        del_max,
        rea_max,
        up_min,
        ov_min,
        do_min):
    """
    Puts face-to-face :
        - a set of ordered genes from genome A
        - a set of ordered genes from genome B.

    Returns a dictionnary for the A-B pair with :
        - 'A' : the name of the object from genome A,
        - 'B' : the name of the object from genome B,
        - 'mapping' : the mapping dictionnary
        - 'stats' : the list of statistics about the mapping.

    Orders genes in 5' -> 3' sens for both genomic regions.
    Makes a transformed version of the two ordered sets where orthologous genes
    between A and B are represented with the same letter.

    Genes in B that do not have the genomic region of A are replaced with a dot
    (see ANY_CHAR).
    Genes in B that are not in the same strand as their ortholog in A are
    flanked with a '(r)' (e.g. 'D(r)' means that the ortholog of gene D from A
    is found in B on the opposite strand.).
    """

    # Each pair of orthologs will be associated with a given letter
    alphabet = list(string.ascii_uppercase)

    # Be sure that A_obj and B_obj have the same number of "gene emplacements"
    for side in ['upstream', 'overlap', 'downstream']:

        A = A_objs[A_obj][side]
        B = B_objs[B_obj][side]

        A_len = len(A)
        B_len = len(B)

        if B_len > A_len:
            A_objs[A_obj][side] = A + [EMPTY_CHAR] * (B_len - A_len)
        elif A_len > B_len:
            B_objs[B_obj][side] = B + [EMPTY_CHAR] * (A_len - B_len)

    # Reverse the upstream genes (currently ordered by increasing distance to
    # the query interval)
    # and thus list upstream -> overlap -> downstream genes in the 5' -> 3'
    # sens.
    A_objs[A_obj]['upstream'].reverse()
    A_objs[A_obj]['genes'] = A_objs[A_obj]['upstream'] + \
        A_objs[A_obj]['overlap'] + A_objs[A_obj]['downstream']

    B_objs[B_obj]['upstream'].reverse()
    B_objs[B_obj]['genes'] = B_objs[B_obj]['upstream'] + \
        B_objs[B_obj]['overlap'] + B_objs[B_obj]['downstream']

    # Add comments
    A_objs[A_obj]['comments'] = ['upstream'] * len(A_objs[A_obj]['upstream']) + ["overlaps " + A_obj] * len(
        A_objs[A_obj]['overlap']) + ['downstream'] * len(A_objs[A_obj]['downstream'])
    B_objs[B_obj]['comments'] = ['upstream'] * len(B_objs[B_obj]['upstream']) + ["overlaps " + B_obj] * len(
        B_objs[B_obj]['overlap']) + ['downstream'] * len(B_objs[B_obj]['downstream'])

    # Each pair of orthologs will be associated with a given letter

    # If there is a problem with A or B, return an empty mapping (no ortholog
    # letters).
    if "FOUND" in A_objs[A_obj]['genes'] or "FOUND" in B_objs[B_obj]['genes']:

        A_objs[A_obj]['letters'] = [EMPTY_CHAR] * len(A_objs[A_obj]['genes'])
        B_objs[B_obj]['letters'] = [EMPTY_CHAR] * len(B_objs[B_obj]['genes'])

    # If there is no problem (normal case) :
    else:

        # For A_obj, just rename genes with alphabet letters, in the 5' -> 3'
        # order
        A_objs[A_obj]['letters'] = []
        alphabet = list(string.ascii_uppercase)
        i = 0
        for gene in A_objs[A_obj]['genes']:

            if gene == EMPTY_CHAR:

                A_objs[A_obj]['letters'].append(EMPTY_CHAR)

            else:

                A_objs[A_obj]['letters'].append(alphabet[i])
                i += 1

        B_objs[B_obj]['letters'] = []
        for gene in B_objs[B_obj]['genes']:

            if gene == EMPTY_CHAR:

                B_objs[B_obj]['letters'].append(EMPTY_CHAR)

            elif gene in orthologs_df['B'].values:

                ortholog = orthologs_df['A'][np.where(
                    orthologs_df['B'] == gene)[0][0]]

                if ortholog in A_objs[A_obj]['genes']:

                    letter = A_objs[A_obj]['letters'][A_objs[A_obj][
                        'genes'].index(ortholog)]

                    if A_strands[ortholog] != B_strands[gene]:

                        B_objs[B_obj]['letters'].append(letter + "(r)")

                    else:
                        B_objs[B_obj]['letters'].append(letter)

                else:

                    B_objs[B_obj]['letters'].append(ANY_CHAR)

            else:

                B_objs[B_obj]['letters'].append(ANY_CHAR)

    # Build a mapping dictionnary
    mapping = {'A position': A_objs[A_obj]['comments'],
               'A gene name': A_objs[A_obj]['genes'],
               'A ortholog': A_objs[A_obj]['letters'],
               'B ortholog': B_objs[B_obj]['letters'],
               'B gene name': B_objs[B_obj]['genes'],
               'B position': B_objs[B_obj]['comments']
               }

    # Get a list of statistics about the mapping
    stats = micro_syteny_mapping_interpreter(mapping,
                                             ins_max,
                                             del_max,
                                             rea_max,
                                             up_min,
                                             ov_min,
                                             do_min)

    return ({'A': A_obj, 'B': B_obj, 'mapping': mapping, 'stats': stats})


def closest_genes(bed, gff_path, n):
    """
    For every interval in the BED, gets the n closest genes up and downstream.

    Mimics the 'betools clostest -a bed -b gff -k n -D'

    Returns a ditionnary with the BED's intervals as keys and as values a
    dictionnary with as keys 'upstream','overlap', and 'downstream' and as
    values the corresponding genes.
    This function takes as input :
    - a BED with intervals of interest
    - a complete GFF3 of the genome
    - 'n' as the number of genes to collect for upstream and downstream genes.
    """

    # Collect the 1000 closest intervals to each object (we cannot filter gene
    # objects directly)
    gff = BedTool(gff_path).sort()
    closest_intervals = BedTool(bed).closest(gff, k=1000, D="ref")

    # Filter and keep only genes
    clst_genes = BedTool(
        [interval for interval in closest_intervals if interval.fields[5] == "gene"])

    # Get the n non-overlapping upstream and downstream genes

    bed_genes_dic = {}

    # This iteration relies on the fact that closest_intervals are sorted by
    # increasing absolute distance
    for line in clst_genes:

        obj_int = '_'.join(line.fields[0:3])

        if obj_int not in bed_genes_dic:

            bed_genes_dic[obj_int] = {
                'upstream': [],
                'overlap': [],
                'downstream': []
            }

            read_on = True

        if read_on:

            # line.count represents the distance in nucleotides, to the object's
            # interval borders

            # If the gene is included in or overlapping the object :
            if line.count == 0:

                bed_genes_dic[obj_int]['overlap'].append(
                    re.search(r".*ID=([^;]+).*", line.fields[11]).group(1))

            # If the gene is upstream the object and we need more upstream
            # genes :
            elif len(bed_genes_dic[obj_int]['upstream']) < n and line.count < 0:

                bed_genes_dic[obj_int]['upstream'].append(
                    re.search(r".*ID=([^;]+).*", line.fields[11]).group(1))

            # If the gene is downstream the object and we need more downstream
            # genes :
            elif len(bed_genes_dic[obj_int]['downstream']) < n and line.count > 0:

                bed_genes_dic[obj_int]['downstream'].append(
                    re.search(r".*ID=([^;]+).*", line.fields[11]).group(1))

            # If he gene was not selected because we have enough for this
            # object :
            elif len(bed_genes_dic[obj_int]['downstream']) >= n:

                read_on = False

    return (bed_genes_dic)


def closest_genes_caller(
        objs_list,
        objs_format,
        gff_path,
        n,
        ortho_list,
        out_path):
    """
    Takes a list of 'ID' values from a GFF, or a list of BED intervals, and for
    each object, gets the n closest genes upstream and downstream of it.

    Since this function must parse the entiere GFF, it also saves some time by
    building meanwhile a dictionnary with genes as key and their strand as
    value. It will be used by the mapping function to tell if two orthologous
    genes are on the same strand or not.
    """

    # This dictionnary will contain the genes' strand in prevision of the
    # orthologs mapping
    orthologs_strand = {}

    # Store problematic objects so they are displayed but with an error
    # message
    pb_objs_genes_dic = {}

    # This list will shrink to seepup the parsing
    ortho_remaining = list(ortho_list[:])
    
    # TODO
    # Try to uniq  objs_list
    objs_list = list(set(objs_list))

    if objs_format == "id":

        # This dictionnary will contain for each object, the list of matching
        # intervals.
        objs_intervals_dic = {obj: [] for obj in objs_list}

        # Parse the entiere GFF file
        with open(gff_path) as gff:

            for line in gff:

                # If line is not a comment
                if line[0] != "#":

                    line = line.strip().split('\t')

                    # Get the 'ID' value
                    results = re.search(r".*ID=([^;]+).*", line[8])
                    if results:

                        ID = results.group(1)

                        if ID in objs_list:

                            objs_intervals_dic[ID].append(line)

                        # In prevision of the orthologs mapping, get their
                        # strand
                        if line[2] == "gene" and ID in ortho_remaining:

                            orthologs_strand[ID] = line[6]
                            ortho_remaining.remove(ID)

        # Try to merge all the intervals
        for obj in objs_intervals_dic:

            if objs_intervals_dic[obj] == []:

                print(
                    "ERROR : there is no features in {} with '{}' as 'ID'.".format(
                        gff_path, obj))

                pb_objs_genes_dic[obj] = {'upstream': ["FEATURE"],
                                          'overlap': ["NOT"],
                                          'downstream': ["FOUND"]
                                          }

            else:

                merged_intervals = BedTool(
                    objs_intervals_dic[obj]).sort().merge(
                    d=1000000000)

                if len(merged_intervals) > 1:

                    print(
                        "ERROR : {} is found in multiple regions : {}".format(
                            obj, objs_intervals_dic[obj]))

                    pb_objs_genes_dic[obj] = {
                        'upstream': ["FEATURE"],
                        'overlap': ["FOUND"],
                        'downstream': ["IN MULTIPLE REGIONS"]}

                else:

                    objs_intervals_dic[obj] = merged_intervals[0]

        # Remove problematic objects
        objs_intervals_dic = {
            obj: objs_intervals_dic[obj] for obj in objs_intervals_dic if not isinstance(
                objs_intervals_dic[obj], list)}

        # Make a unique BED with all merged intervals for which we want to retrieve the
        # n closest genes upstream and downstream
        objs_intervals_bed = BedTool(list(objs_intervals_dic.values())).sort()

        # Retrieve these genes
        bed_genes_dic = closest_genes(objs_intervals_bed, gff_path, n)

        # Generate a wrap up dictionnary
        objs_genes_dic = {obj: bed_genes_dic['_'.join(
            objs_intervals_dic[obj])] for obj in objs_intervals_dic}

    elif objs_format == "bed":

        # This object is created if somme intervals need to be corrected
        new_objs_list = False

        try:

            # The BED verison of the intervals of interest
            ints_bed = BedTool(objs_list).sort()

        except BaseException:

            print("WARNING : some BED-style objects seems malformed. Maybe some starts are greater than stop.\nTrying to correct that...")

            new_objs_list = []
            for obj in objs_list:

                splitted = obj.split("\t")

                # If indeed, the start is greater that the stop, invert them
                if int(splitted[1]) > int(splitted[2]):

                    new_splitted = [splitted[0], splitted[2], splitted[1]]
                    new_objs_list.append("\t".join(new_splitted))

                else:

                    new_objs_list.append(obj)

            # The BED verison of the intervals of interest
            ints_bed = BedTool(new_objs_list).sort()

            # Store the new objects to a file so the user can visualize the
            # changes that were made.
            corrected_objs_path = out_path + "_corrected_inputs.tsv"

#            corrected_objs_lines = [ ['\t'.join([objs_list[i], '-->', new_objs_list[i]])] for i in range(len(objs_list))]
            corrected_objs_lines = [ [ objs_list[i], '-->', new_objs_list[i] ] for i in range(len(objs_list))]

            with open(corrected_objs_path, mode='a', newline='') as cor:

                fieldnames = [
                    "original seq",
                    "original start",
                    "original stop",
                    "new_seq",
                    "new_start",
                    "new_stop"]

                writer = csv.writer(cor, delimiter="\t", lineterminator="\n")

                writer.writerow(fieldnames)

                for line in corrected_objs_lines:

                    writer.writerow(line)

                # Add a trailling line
                cor.write("\n")

            print('Done. Modified objects can be visualized in "{}".'.format(
                corrected_objs_path))

        # Parse the entiere GFF file
        gff_seqs = []

        with open(gff_path) as gff:

            for line in gff:

                # If line is not a comment
                if line[0] != "#":

                    line = line.strip().split('\t')

                    if line[0] not in gff_seqs:

                        gff_seqs.append(line[0])

                    # Get the 'ID' value
                    results = re.search(r".*ID=([^;]+).*", line[8])
                    if results:

                        ID = results.group(1)

                        # In prevision of the orthologs mapping, get their
                        # strand
                        if line[2] == "gene" and ID in ortho_remaining:

                            orthologs_strand[ID] = line[6]
                            ortho_remaining.remove(ID)

        # Test if intervals' sequence is indeed in the GFF file.
        for interval in ints_bed:

            if interval.chrom not in gff_seqs:

                print(
                    "ERROR : {} is not found in : {}".format(
                        interval.chrom, gff_path))

                pb_objs_genes_dic['_'.join(interval)] = {'upstream': ["SEQUENCE"],
                                                         'overlap': ["NOT"],
                                                         'downstream': ["FOUND"]
                                                         }

        # Retrieve these genes
        objs_genes_dic = closest_genes(ints_bed, gff_path, n)


        # In case some objects were corrected (start greater than stop),
        # retrieve their original name
        if new_objs_list:
            
            for i in range(len(objs_list)):

                old_int = objs_list[i].replace("\t", "_")
                new_int = new_objs_list[i].replace("\t", "_")

                if old_int != new_int and new_int in objs_genes_dic :
                    print(objs_genes_dic.keys())
                    print("old : {}, new : {}".format(old_int, new_int))
                    uncorrected_int = objs_genes_dic[new_int]
                    objs_genes_dic[old_int] = uncorrected_int
                    del objs_genes_dic[new_int]


    # Merge the two dictionnaries
    objs_genes_dic = {**objs_genes_dic, **pb_objs_genes_dic}

    return ([objs_genes_dic, orthologs_strand])


def input_formatter(list_path, gffB_path):
    """
    Detects how to deal with the input_list file.
    """

    try :
        input_list = pd.read_csv(list_path, sep="\t", header=None,  )
    
    except pd.errors.EmptyDataError:
        print("WARNING : '{}' is empty !".format(list_path))
        sys.exit()
    
    nb_col = len(input_list.columns)

    if nb_col == 6:

        print("Detected format for {} : both A objects and B objects are in a BED-style format.".format(list_path))
        input_list.columns = [
            "A_seq",
            "A_start",
            "A_end",
            "B_seq",
            "B_start",
            "B_end"]
        input_list["A"] = input_list["A_seq"].astype(
            str) + "_" + input_list["A_start"].astype(str) + "_" + input_list["A_end"].astype(str)
        input_list["B"] = input_list["B_seq"].astype(
            str) + "_" + input_list["B_start"].astype(str) + "_" + input_list["B_end"].astype(str)

    elif nb_col == 2:

        print("Detected format for {} : both A objects and B objects are 'ID' values.".format(
            list_path))
        input_list.columns = ["A", "B"]

    elif nb_col == 4:

        # Columns types
        type2 = input_list[1].dtype.kind
        type3 = input_list[2].dtype.kind
        type4 = input_list[3].dtype.kind

        if type2 == 'i' and type3 == 'i' and type4 != 'i':

            print("Detected format for {} : A objects are in a BED-style format, and B objects are 'ID' values.".format(list_path))
            input_list.columns = ["A_seq", "A_start", "A_end", "B"]
            input_list["A"] = input_list["A_seq"].astype(
                str) + "_" + input_list["A_start"].astype(str) + "_" + input_list["A_end"].astype(str)

        elif type2 != 'i' and type3 == 'i' and type4 == 'i':

            print("Detected format for {} : A objects are 'ID' values, and B objects are in a BED-style format.".format(list_path))
            input_list.columns = ["A", "B_seq", "B_start", "B_end"]
            input_list["B"] = input_list["B_seq"].astype(
                str) + "_" + input_list["B_start"].astype(str) + "_" + input_list["B_end"].astype(str)

        # In this case, either sequence names of B  or 'ID' values of B are
        # intergers.
        elif type2 == 'i' and type3 == 'i' and type4 == 'i':

            print("WARNING : as cols 2,3, and 4 are all integers, {}'s format is not obvious.\nEither genome B's objects are in BED-style and sequence names are intergers, or genome B's objects are 'ID' values but also integers.\nPerforming an additional analysis to guess the right format...".format(list_path))

            # Try to guess :

            # Get the list of sequences and 'ID' values in GFFB :
            gffB_seqs = []
            gffB_IDs = []

            gffB = BedTool(gffB_path)
            for interval in gffB:

                gffB_seqs.append(interval.fields[0])

                results = re.search(r".*ID=([^;]+).*", interval.fields[8])
                if results:

                    ID = results.group(1)
                    gffB_IDs.append(ID)

            gffB_seqs = list(set(gffB_seqs))
            gffB_IDs = list(set(gffB_IDs))

            # Test column 2
            found_seqs = len(
                [value for value in input_list[1] if str(value) in gffB_seqs])
            seq_prct = int(100 * found_seqs / len(input_list[1]))
            # Test column 4
            found_IDs = int(
                100 * len([value for value in input_list[3] if str(value) in gffB_IDs]))
            ID_prct = found_IDs / len(input_list[3])

            print(
                "Done.\n{}% of values in the second column of {} was found among genome B's sequences.\n{}% of values in the fourth column were found among genome B's 'ID' values.\nAs a consequence : ".format(
                    seq_prct,
                    list_path,
                    ID_prct))

            if found_seqs > found_IDs:

                print(
                    "Detected format for {} : A objects are 'ID' values, and B objects are in a BED-style format.".format(list_path))
                input_list.columns = ["A", "B_seq", "B_start", "B_end"]
                input_list["B"] = input_list["B_seq"].astype(
                    str) + "_" + input_list["B_start"].astype(str) + "_" + input_list["B_end"].astype(str)

            else:

                print(
                    "Detected format for {} : A objects are in a BED-style format, and B objects are 'ID' values.".format(list_path))
                input_list.columns = ["A_seq", "A_start", "A_end", "B"]
                input_list["A"] = input_list["A_seq"].astype(
                    str) + "_" + input_list["A_start"].astype(str) + "_" + input_list["A_end"].astype(str)

        else:

            sys.exit(
                "Format problem for {} : with 4 columns, columns 2 and 3 or columns 3 and 4 should be intergers.".format(list_path))

    else:

        sys.exit(
            "The number of columns provided in {} is not correct. It should have 2, 4, or 6 columns.".format(list_path))

    return (input_list)


def confirm_choice(message="Do you want to continue?"):
    """
    A wrapup function to prompt a confirm_choice question.
    """

    while True:
        answer = input(message)
        if answer.lower() in ["y", "yes"]:
            return (True)
        elif answer.lower() in ["n", "no"]:
            return (False)
        else:
            print("Enter either yes/no")


def check_synteny():
    """
    This is the main function of this program.

    - It first gets the user arguments and load the required files,

    - Then it formats the objects from genome A and B with input_formatter()
    (they can be 'ID' value from the GFF3 or BED-style genomic regions),

    - It then gives the object list to closest_genes_caller() which return a
    dictionnary per genome with all the genes that are close to the input
    objects,

    - Next every pair of objects A-B is given to micro_synteny_mapper(), along
    with their closest genes. The function returns for each pair a dictionnary
    with the name of the two objects, their 'mapping', and a bunch of numbers
    about this mapping (stats).

    - A the end, the main function writes two output file :
        - "x.out" which contains the mappings for each pair of objects
        (this is convenient to visualize two genomic regions face to face, and
        be can used for a manul checking),
        - "x.tsv" with one row per pair, with its stats. This is were the user
        can quickly check which pairs of objects are indeed in synteny.
    """

    logo = r"""  ____  _                  _      ____                 _
 / ___|| |__    ___   ___ | | __ / ___|  _   _  _ __  | |_   ___  _ __   _   _
| |    | '_ \  / _ \ / __|| |/ / \___ \ | | | || '_ \ | __| / _ \| '_ \ | | | |
| |___ | | | ||  __/| (__ |   <   ___) || |_| || | | || |_ |  __/| | | || |_| |
 \____||_| |_| \___| \___||_|\_\ |____/  \__, ||_| |_| \__| \___||_| |_| \__, |
                                         |___/                           |___/
"""
    print(logo)

    # Argument parsing
    p = configargparse.ArgParser(default_config_files=['check_synteny.config', 'check_synteny.yaml'])
    p.add('--erase',                help="if an output file already exists, erase it", action='store_true')
    p.add('--ncpus',     type=int,  help="number of cpus to use", default=1)
    p.add('--flankA',    type=int,  help="number of flanking genes to consider for A (max 26)", default=2)
    p.add('--flankB',    type=int,  help="number of flanking genes to consider for B (max 26)", default=4)
    p.add('--up_min',    type=int,  help="minimum number of upstream orthologs in B", default=1)
    p.add('--ov_min',    type=int,  help="minimum number of upstream orthologs in B", default=0)
    p.add('--do_min',    type=int,  help="minimum number of upstream orthologs in B", default=1)
    p.add('--insertion', type=int,  help="insertion threshold", default=0)
    p.add('--ins_max',   type=int,  help="maximum number of insertion in B", default=-1)
    p.add('--del_max',   type=int,  help="maximum number of insertion in B", default=-1)
    p.add('--rea_max',   type=int,  help="maximum number of insertion in B", default=-1)

    # DEV MODE
    if "DEV_MODE" in locals() and DEV_MODE:
        print("\nWARNING : DEV_MODE = True ! \n")
        args = p.parse_args()
        args.gffA = HOME + "PIPELINE/SYNTENY/TEST_FILES/Hsap_reduced.gff3"
        args.gffB = HOME + "PIPELINE/SYNTENY/TEST_FILES/Ppan_reduced.gff"
        args.ortho = HOME + "PIPELINE/SYNTENY/TEST_FILES/ortho.tsv"
        args.list = HOME + "PIPELINE/SYNTENY/TEST_FILES/input_list.tsv"
        args.out = HOME + "PIPELINE/SYNTENY/output"

    else:
        # Argument parsing   
        p.add('--config', is_config_file=True, help="a config file with any parameter with .ini or .yaml style syntax (eg. key=value or key: value). Command line parameters OVERRIDE config file parameters.")
        p.add('--gffA',  required=True, help="GFF3 for genomeA")  # this option can be set in a config file because it starts with '--'
        p.add('--gffB',  required=True, help="GFF3 for genomeB")
        p.add('--list',  required=True, help="a .tsv file with the genomic objects to test (first column(s) : object from A, next column(s) : object from B. An object can be any feature with an 'ID' attribute (represented by a single column), or an genomic interval in BED-style format (3 columns): sequence, start(ZERO-BASED), end")
        p.add('--ortho', required=True, help="a .tsv file with the pairs of orthologous genes (first column : gene from A, second column : gene from B)")
        p.add('--out',   required=True, help="output name without extension (a 'X.out' and 'X_summary.tsv' will be produced.")

        args = p.parse_args()

    # Load the input files

    # GFF files
    gffA_path = args.gffA

    gffB_path = args.gffB

    # List of orthologous genes from A and B.
    orthologs_path = args.ortho
    try : 
        orthologs_df = pd.read_csv(
            orthologs_path,
            sep="\t",
            header=None,
            names=[
                "A",
                "B"],
            dtype=str)
    except :
        print("ERROR : could not parse {}.".format(orthologs_path))
        sys.exit()


    # List of objects from genome A and their supposed syntenic object from
    # genome B on the same row.
    list_path = args.list

    # Path without extension used to create the ouput files.
    out_path = args.out

    if not args.erase:

        for file in [out_path + '.out', out_path + '.tsv']:

            if os.path.exists(file):

                erase = confirm_choice(
                    "WARNING : '{}' already exists. Do you want to continue? ".format(file))

                if not erase:
                    sys.exit()
    else :
        erase = args.erase

    # Minus one stands as infinity
    if args.ins_max == -1:
        args.ins_max = math.inf
    if args.del_max == -1:
        args.del_max = math.inf
    if args.rea_max == -1:
        args.rea_max = math.inf

    # Detect the input_list format.
    # gffB is for dealing with exceptions.
    input_list = input_formatter(list_path, gffB_path)
        
    # Send A and B objects to closest_genes_caller, according to their format
    if "A_seq" in input_list.columns:  # .i.e A objects are BED intervals

        A_format = 'bed'
        input_list = input_list.astype(
            {"A_seq": str, "A_start": str, "A_end": str})
        A_objs = input_list.loc[:, [
            'A_seq', 'A_start', 'A_end']].values.tolist()
        A_objs = ['\t'.join(ele) for ele in A_objs]

    else:  # .i.e A objects are 'ID' values from gffA

        A_format = 'id'
        A_objs = list(input_list['A'])

    if "B_seq" in input_list.columns:

        B_format = 'bed'
        input_list = input_list.astype(
            {"B_seq": str, "B_start": str, "B_end": str})
        B_objs = input_list.loc[:, [
            'B_seq', 'B_start', 'B_end']].values.tolist()
        B_objs = ['\t'.join(ele) for ele in B_objs]

    else:

        B_format = 'id'
        B_objs = list(input_list['B'])

    # For each objects from genome A and genome B, get a dictionnary with the n
    # upstream and downstream closest genes, plus genes overlapped by the
    # object.
    # Also get a dictionnay with genes and their DNA strand for the mapping.
    A_genes, A_orthologs_strand = closest_genes_caller(
        A_objs, A_format, gffA_path, args.flankA, orthologs_df['A'].values, out_path)
    B_genes, B_orthologs_strand = closest_genes_caller(
        B_objs, B_format, gffB_path, args.flankB, orthologs_df['B'].values, out_path)

    # For every pair of objects from genome A and B, make a mapping dictionnary
    # .
    # TODO : more elegant way to get a master dictionnary at the end
    valid_pairs_df = input_list[input_list['A'].isin(
        A_genes) & input_list['B'].isin(B_genes)]
    master_dict = valid_pairs_df.apply(
        lambda x: micro_synteny_mapper(
            A_genes,
            x['A'],
            B_genes,
            x['B'],
            orthologs_df,
            A_orthologs_strand,
            B_orthologs_strand,
            args.ins_max,
            args.del_max,
            args.rea_max,
            args.up_min,
            args.ov_min,
            args.do_min),
        axis=1).to_dict()

    # Write the output files

    # .out file
    with open(out_path + '.out', mode='w', newline='') as out:

        for pair in master_dict:

            out.write(
                ">{} vs {}\n".format(
                    master_dict[pair]['A'],
                    master_dict[pair]['B']))

            mapping = master_dict[pair]['mapping']

            table = pt.PrettyTable()

            for key in mapping:

                table.add_column(key, mapping[key])

            print(table, file=out)

            # Add a line break after each mapping
            out.write("\n")

    # TSV file
    with open(out_path + ".tsv", mode='w', newline='') as tsv:

        fieldnames = ['object from genome A',
                      'object from genome B',
                      '#insertions',
                      '#deletions',
                      '#rearrangments',
                      '#ortho upstream',
                      '#ortho overlap',
                      '#ortho downstream',
                      'isSyntenic']

        writer = csv.writer(tsv, delimiter="\t", lineterminator="\n")

        writer.writerow(fieldnames)

        for pair in master_dict:

            writer.writerow([master_dict[pair]['A'],
                             master_dict[pair]['B']] + master_dict[pair]['stats'])

    print(
        "\nThe results can be found in '{}' (visual mappings) and '{}' (stats).".format(
            out_path +
            '.out',
            out_path +
            '.tsv'))

    return ([A_genes, B_genes, A_orthologs_strand,
            B_orthologs_strand, master_dict])


A_genes, B_genes, A_orthologs_strand, B_orthologs_strand, master_dict = check_synteny()
