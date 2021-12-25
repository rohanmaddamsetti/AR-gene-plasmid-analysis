#!/usr/bin/env python

'''
count-cds.py by Rohan Maddamsetti.

This script goes through the headers in ../results/protein_db.faa,
makes a dictionary of NCBI_Nucleotide_Accessions to number of sequences
in the database with that accession.

Usage: python count-cds.py > ../results/protein_db_CDS_counts.csv

'''

cds_counts = {}

## some genomes have artifactual duplicate CDS in their annotation
## (at least 10 Xanthomonas strains as of 12/24/2021).
## To handle these cases, make sure that the location of each cds per accession
## has not been seen before.

accession_to_CDS_locs = {} ## This is a dict of dicts for O(1) lookups.

with open("../results/protein_db.faa","r") as protein_db_fh:
    for line in protein_db_fh:
        ## skip unless it's a sequence header
        if not line.startswith('>'): continue
        fields = line.split()
        ncbi_nucleotide_accession = fields[0].split('|')[-1].split('_prot')[0]
        prot_location = fields[-2].strip("[]").split('=')[-1]

        ## check for artifactual duplicates
        ## (same location, found multiple times).
        if ncbi_nucleotide_accession in accession_to_CDS_locs:
            episome_dict = accession_to_CDS_locs[ncbi_nucleotide_accession]
            if prot_location in episome_dict: ## then we have seen this protein
                continue ## in this episome already-- skip this artifact.
            else: ## mark that we are seeing this protein for the first time.
                episome_dict.update({prot_location:True})
        else:
            accession_to_CDS_locs[ncbi_nucleotide_accession] = {}

        ## after passing error checking, update the number of cds in this episome.
        if ncbi_nucleotide_accession in cds_counts:
            cds_counts[ncbi_nucleotide_accession] += 1
        else:
            cds_counts[ncbi_nucleotide_accession] = 1

print("NCBI_Nucleotide_Accession,CDS_count")
for k,v in cds_counts.items():
    print(','.join([k,str(v)]))


