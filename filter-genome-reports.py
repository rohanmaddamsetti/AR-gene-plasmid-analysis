#!/usr/bin/env python

''' 
filter-genome-reports.py by Rohan Maddamsetti.

This script goes through the genome report in data/prokaryotes.txt, and filters
for lines with either complete genomes or both plasmids and chromosomes.

Usage: python filter-genome-reports.py > ../results/best-prokaryotes.txt

'''

with open("../data/GENOME_REPORTS/prokaryotes.txt","r") as g_report:
    for i, line in enumerate(g_report):
        line = line.strip()
        if i == 0: ## print the header.
            print(line)
        else:
            if ("Complete Genome" in line) or (('plasmid' in line) and ('chromosome' in line)):
                print(line)
