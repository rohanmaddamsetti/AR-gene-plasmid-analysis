#!/usr/bin/env python

''' 
filter-genome-reports.py by Rohan Maddamsetti.

This script goes through the genome report in data/prokaryotes.txt, and filters
for lines with plasmids.

Usage: python filter-genome-reports.py > ../results/prokaryotes-with-plasmids.txt

'''

with open("../data/prokaryotes.txt","r") as g_report:
    for i, line in enumerate(g_report):
        line = line.strip()
        if i == 0: ## print the header.
            print(line)
        else:
            if ('plasmid' in line) and ('chromosome' in line):
                print(line)
