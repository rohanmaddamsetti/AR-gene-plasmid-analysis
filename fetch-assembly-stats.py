#!/usr/bin/env python

'''
fetch-assembly-stats.py by Rohan Maddamsetti.

This script reads in ../results/best-prokaryotes.txt.

NOTE: for path names to be processed properly, this script must be run
from the src/ directory as:  "python fetch-assembly-stats.py".
'''

import urllib.request
from os.path import basename, exists
import gzip
import os
from tqdm import tqdm

## open the genome report file, and parse line by line.
with open("../results/best-prokaryotes.txt", "r") as genome_report_fh:
    for i, line in enumerate(tqdm(genome_report_fh)):
        line = line.strip()
        if i == 0: ## get the names of the columns from the header.
            column_names_list = line.split('\t')
            continue ## don't process the header further.
        fields = line.split('\t')
        ftp_path = fields[20]

        ## Now download the assembly stats file  if it doesn't exist on disk.
        QCstats_ftp_path = ftp_path + '/' + basename(ftp_path) + "_assembly_stats.txt"
        QCstats_fname = "../results/refseq-assembly-statistics/" + basename(ftp_path) + "_assembly_stats.txt"
        if exists(QCstats_fname): continue ## no need to get it if we already have it.
        fetch_attempts = 5
        not_fetched = True
        while not_fetched and fetch_attempts:
            try:
                urllib.request.urlretrieve(QCstats_ftp_path, filename=QCstats_fname)
                not_fetched = False ## assume success if the previous line worked,
                fetch_attempts = 0 ## and don't try again.
            except urllib.error.URLError: ## if some problem happens, try again.
                fetch_attempts -= 1


