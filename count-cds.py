#!/usr/bin/env python

'''
count-cds.py by Rohan Maddamsetti.

This script goes through the headers in ../results/protein_db.faa,
makes a dictionary of NCBI_Nucleotide_Accessions to number of sequences
in the database with that accession.

Usage: python count-cds.py > ../results/protein_db_CDS_counts.csv

'''

cds_counts = {}

with open("../results/protein_db.faa","r") as protein_db_fh:
    for line in protein_db_fh:
        ## skip unless it's a sequence header
        if not line.startswith('>'): continue
        fields = line.split()
        ncbi_nucleotide_accession = fields[0].split('|')[-1].split('_prot')[0]
        if ncbi_nucleotide_accession in cds_counts:
            cds_counts[ncbi_nucleotide_accession] += 1
        else:
            cds_counts[ncbi_nucleotide_accession] = 1

print("NCBI_Nucleotide_Accession,CDS_count")
for k,v in cds_counts.items():
    print(','.join([k,str(v)]))


