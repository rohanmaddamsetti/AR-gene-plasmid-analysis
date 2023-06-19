#!/usr/bin/env python

'''
count-cds.py by Rohan Maddamsetti.

This script goes through the gbk annotation files,
makes a dictionary of NCBI_Nucleotide_Accessions to number of sequences
in the database with that accession.

Usage: python count-cds.py

'''

import os
from tqdm import tqdm
from Bio import SeqIO

## some genomes have artifactual duplicate CDS in their annotation
## (at least 10 Xanthomonas strains as of 12/24/2021).
## To handle these cases, make sure that the location of each cds per accession
## has not been seen before.

cds_counts = {}
outf = "../results/protein_db_CDS_counts.csv"

gbk_annotation_dir = "../results/gbk-annotation/"
gbk_files = [x for x in os.listdir(gbk_annotation_dir) if x.endswith("_genomic.gbff")]

for x in tqdm(gbk_files):
    gbk_path = os.path.join(gbk_annotation_dir, x)
    annotation_accession = x.split("_genomic.gbff")[0]
    genome_CDS_count = 0
    with open(gbk_path,'rt') as gbk_fh:
        for replicon in SeqIO.parse(gbk_fh, "gb"):
            ## check for artifactual duplicates
            ## (same location, found multiple times).
            location_set = set()
            for feat in replicon.features:
                if feat.type != "CDS": continue
                prot_location = str(feat.location)
                if prot_location in location_set:
                    continue
                else:
                    location_set.add(prot_location)
                    genome_CDS_count += 1
    cds_counts[annotation_accession] = genome_CDS_count


with open(outf, "w") as out_fh:
    header = "NCBI_Nucleotide_Accession,CDS_count"
    out_fh.write(header + "\n")
    for k,v in cds_counts.items():
        row = ','.join([k,str(v)])
        out_fh.write(row + "\n")


