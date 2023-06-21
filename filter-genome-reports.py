#!/usr/bin/env python

''' 
filter-genome-reports.py by Rohan Maddamsetti.

This script goes through the genome report in data/prokaryotes.txt,
filters for lines with complete genomes, and IMPORTANTLY,
replaces "GCA" with "GCF" throughout, so that only data deposted in NCBI RefSeq
is analyzed.

See the documentation here for details:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/troubleshooting/faq/#what-is-the-difference-between-a-genbank-gca-and-refseq-gcf-genome-assembly

Usage: python filter-genome-reports.py > ../results/best-prokaryotes.txt
'''

with open("../data/GENOME_REPORTS/prokaryotes.txt","r") as g_report:
    for i, line in enumerate(g_report):
        line = line.strip()
        if i == 0: ## print the header.
            print(line)
        else:
            if ("Complete Genome" in line):
                updated_line = line.replace("GCA", "GCF")
                fields = line.split('\t')
                ftp_path = fields[20]
                if ftp_path == "-": continue ## skip rows that don't have an FTP URL.
                print(updated_line)
