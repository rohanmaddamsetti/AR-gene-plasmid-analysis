#!/usr/bin/env python

'''
fetch-gbk-annotation.py by Rohan Maddamsetti.

This script reads in ../results/AR-gene-duplication/prokaryotes-with-plasmids.txt.

NOTE: for path names to be processed properly, this script must be run
from the src/ directory as python fetch-gbk-annotation.py.
'''

import urllib.request
from os.path import basename, exists
import gzip
import os
from tqdm import tqdm

## open the genome report file, and parse line by line.
with open("../results/prokaryotes-with-plasmids.txt", "r") as genome_report_fh:
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
        if not exists(gbff_gz_fname):
            try:
                urllib.request.urlretrieve(gbk_ftp_path, filename=gbff_gz_fname)
            except urllib.error.URLError: ## some problem happens
                try: ## try one more time, in case there was some connection problem.
                    urllib.request.urlretrieve(gbk_ftp_path,filename=gbff_gz_fname)
                except: ## skip if still can't get it.
                    continue


