#!/usr/bin/env python

"""
cross-check-Hawkey2022-accessions.py by Rohan Maddamsetti.

This script prints out any accessions in the Hawkey data that match to analyzed genomes.

"""

hawkey_ids = []
with open("../data/Hawkey2022-Hybrid-Assemblies-NCBI-BioProject-PRJNA646837/PRJNA646837_AssemblyDetails.txt", "r") as hawkey_accessions_fh:
    for line in hawkey_accessions_fh:
        if line.startswith("#"): continue
        my_id = line.split()[0].replace("GCA", "GCF")
        hawkey_ids.append(my_id)
        

genome_id_to_accession = dict()
genome_ids = []
with open("../results/gbk-annotation-of-analyzed-complete-genomes.csv", "r") as main_genomes_fh:
    for i,line in enumerate(main_genomes_fh):
        if i == 0: continue ## skip the header
        accession = line.split(",")[0]
        genome_id = accession.rpartition("_")[0]
        genome_ids.append(genome_id)
        genome_id_to_accession[genome_id] = accession


for x in hawkey_ids:
    if x in genome_ids:
        print(genome_id_to_accession[x])
