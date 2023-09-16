#!/usr/bin/env python

"""
filter-MOBTyper-results.py by Rohan Maddamsetti

This script filters the MOBTyper results that Hye-in Son generated based on the complete genomes
examined in this project.

Usage: python filter-MOBTyper-results.py > ../results/FileS5-MOBTyper-plasmid-annotations.tsv
"""

genome_accessions = []

with open("../results/chromosome-plasmid-table.csv") as genome_accession_fh:
    for i, line in enumerate(genome_accession_fh):
        line = line.strip()
        if i == 0: continue ## skip the header.
        annotation_accession = line.split(",")[-1]
        genome_accessions.append(annotation_accession)

## now remove duplicates and sort.
genome_accessions = list(set(genome_accessions))
genome_accessions.sort()

## now filter the MOBTyper results that Hye-in Son generated.
with open("../data/unique_mob_results_from_Hyein_Son.txt") as mobtyper_results_fh:
    for i, line in enumerate(mobtyper_results_fh):
        line = line.strip()
        if i == 0: ## print the header
            print(line)
        else:
            for annotation_accession in genome_accessions:
                if annotation_accession in line:
                    print(line)
