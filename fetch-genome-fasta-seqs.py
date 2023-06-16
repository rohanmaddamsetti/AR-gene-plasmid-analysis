#!/usr/bin/env python

'''
fetch-genome-fasta-seqs.py by Rohan Maddamsetti.

This script reads in ../results/best-prokaryotes.txt,
and downloads the corresponding fasta seqs.

NOTE: for path names to be processed properly, this script must be run
from the src/ directory as:  "python fetch-genome-fasta-seqs.py".

NOTE: This script should not be run directly. Instead, run fetch-and-dereplicate-seqs.sh,
which runs this script.
'''

import urllib.request
from os.path import basename, exists
import gzip
import os
from tqdm import tqdm


def download_FASTA_genomes(genome_report_file):
    ## open the genome report file, and parse line by line.
    with open(genome_report_file, "r") as genome_report_fh:
        for i, line in enumerate(tqdm(genome_report_fh)):
            line = line.strip()
            if i == 0: ## get the names of the columns from the header.
                column_names_list = line.split('\t')
                continue ## don't process the header further.
            fields = line.split('\t')
            ftp_path = fields[20]

            annotation_accession = basename(ftp_path)
            ## Now download the fasta files  if it doesn't exist on disk.
            fasta_ftp_path = ftp_path + '/' + annotation_accession + "_genomic.fna.gz"
            fasta_fname = "../results/genome-fasta-files/" + annotation_accession + "_genomic.fna.gz"
            if exists(fasta_fname): continue ## no need to get it if we already have it.
            fetch_attempts = 5
            not_fetched = True
            while not_fetched and fetch_attempts:
                try:
                    urllib.request.urlretrieve(fasta_ftp_path, filename=fasta_fname)
                    not_fetched = False ## assume success if the previous line worked,
                    fetch_attempts = 0 ## and don't try again.
                except urllib.error.URLError: ## if some problem happens, try again.
                    fetch_attempts -= 1
    return None


genome_report_file = "../results/best-prokaryotes.txt"
download_FASTA_genomes(genome_report_file)
