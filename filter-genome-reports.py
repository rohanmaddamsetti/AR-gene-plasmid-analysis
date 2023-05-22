#!/usr/bin/env python

''' 
filter-genome-reports.py by Rohan Maddamsetti.

This script goes through the genome report in data/prokaryotes.txt, and filters
for lines with complete genomes.

Usage: python filter-genome-reports.py > ../results/best-prokaryotes.txt
'''

with open("../data/GENOME_REPORTS/prokaryotes.txt","r") as g_report:
    for i, line in enumerate(g_report):
        line = line.strip()
        if i == 0: ## print the header.
            print(line)
        else:
            if ("Complete Genome" in line):
                print(line)
