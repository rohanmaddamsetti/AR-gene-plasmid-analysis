#!/usr/bin/env python

'''
fetch-gbk-annotation.py by Rohan Maddamsetti.

This script reads in ../results/best-prokaryotes.txt.

NOTE: for path names to be processed properly, this script must be run
from the src/ directory as python fetch-gbk-annotation.py.
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

        ## Now download the Genbank annotation if it doesn't exist on disk.
        gbk_ftp_path = ftp_path + '/' + basename(ftp_path) + "_genomic.gbff.gz"
        gbff_gz_fname = "../results/gbk-annotation/" + basename(ftp_path) + "_genomic.gbff.gz"
        if exists(gbff_gz_fname): continue ## no need to get it if we already have it.
        gbk_fetch_attempts = 5
        gbk_not_fetched = True
        while gbk_not_fetched and gbk_fetch_attempts:
            try:
                urllib.request.urlretrieve(gbk_ftp_path, filename=gbff_gz_fname)
                gbk_not_fetched = False ## assume success if the previous line worked,
                gbk_fetch_attempts = 0 ## and don't try again.
            except urllib.error.URLError: ## if some problem happens, try again.
                gbk_fetch_attempts -= 1


