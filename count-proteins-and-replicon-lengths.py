#!/usr/bin/env python

"""
count-proteins-and-replicon-lengths.py by Rohan Maddamsetti.

Usage: python count-proteins-and-replicon-lengths.py

This script tabulates the length and number of proteins for all chromosomes and plasmids in these complete genomes.

This is particularly important, to filter out genomes with erroneous plasmid and chromosome annotation in RefSeq!
I found that some 'plasmids' in RefSeq, like  NZ_CP104448.1, are actually mis-annotated chromosomes, and
the annotated chromosome is actually a plasmid!

These data are needed to filter out genomes where the chromosome is smaller than some plasmids-- this is most likely
caused by misannotation in the RefSeq database.

"""

import os
from Bio import SeqIO
from tqdm import tqdm


def tabulate_replicon_lengths_and_proteins(genomes_dir, replicon_length_csv_file):
    with open(replicon_length_csv_file, 'w') as outfh:
        header = "AnnotationAccession,SeqID,SeqType,replicon_length,protein_count\n"
        outfh.write(header)
        for gbk in tqdm(os.listdir(genomes_dir)):
            if not gbk.endswith(".gbff"): continue
            annotation_accession = gbk.split("_genomic.gbff")[0]
            infile = os.path.join(genomes_dir, gbk)
            with open(infile, "r") as genome_fh:
                for replicon in SeqIO.parse(genome_fh, "gb"):
                    SeqID = replicon.id
                    if "plasmid" in replicon.description:
                        SeqType = "plasmid"
                    elif "complete" in replicon.description or "chromosome" in replicon.description:
                        SeqType = "chromosome"
                    else:
                        SeqType = "unknown"
                    replicon_length = str(len(replicon))
                    ## count the number of proteins in the replicon.
                    protein_count = 0
                    for feature in replicon.features:
                        if feature.type == "CDS" and "translation" in feature.qualifiers:
                            protein_count += 1
                    ## now write out the data for the replicon.
                    row = ','.join([annotation_accession, SeqID, SeqType, str(replicon_length), str(protein_count)])
                    outfh.write(row + "\n")
    return


def main():
    genomes_dir = "../results/gbk-annotation/"
    replicon_length_csv_file = "../results/replicon-lengths-and-protein-counts.csv"
    tabulate_replicon_lengths_and_proteins(genomes_dir, replicon_length_csv_file)
    return


if __name__ == "__main__":
    main()
