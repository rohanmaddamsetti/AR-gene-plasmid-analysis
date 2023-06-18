#!/usr/bin/env python

'''
fetch-genome-and-plasmid-cds.py by Rohan Maddamsetti.

This script reads in ../results/best-prokaryotes.txt.

NOTE: for path names to be processed properly, this script must be run
from the src/ directory as python fetch-genome-and-plasmid-cds.py.

'''

import urllib.request
from os.path import basename
import gzip
import os
from tqdm import tqdm

## create an empty protein database file.
with open("../results/protein_db.faa", "w") as protein_db_fh:
    ## open the genome report file, and parse line by line.
    with open("../results/best-prokaryotes.txt", "r") as genome_report_fh:
        for i, line in enumerate(tqdm(genome_report_fh)):
            line = line.strip()
            if i == 0: ## get the names of the columns from the header.
                column_names_list = line.split('\t')
                continue ## don't process the header further.
            fields = line.split('\t')
            ftp_path = fields[20]
            print("CURRENT FTP PATH: ", ftp_path)
            ## NOW GET PROTEINS.
            ## download the gzipped translated cds (amino acid sequences).
            tr_cds_ftp_path = ftp_path + '/' + basename(ftp_path) + "_translated_cds.faa.gz"
            protein_fetch_attempts = 5
            proteins_not_fetched = True
            while proteins_not_fetched and protein_fetch_attempts: 
                try:
                    urllib.request.urlretrieve(tr_cds_ftp_path,filename="../results/temp_tr_cds.faa.gz")
                    ## if the previous line completed successfully, then unzip,
                    ## concatenate to the protein database file, and delete the temporary file.
                    with gzip.open("../results/temp_tr_cds.faa.gz", "rt") as temp_tr_cds_infile:
                        protein_db_fh.write(temp_tr_cds_infile.read())
                    proteins_not_fetched = False ## assume success if the previous lines worked,
                    protein_fetch_attempts = 0 ## and don't try again.
                except: ## if some problem occurs during the download or decompression, try again.
                    protein_fetch_attempts -= 1
                ## if the temporary file exists, clean up before trying again or moving on.
                if os.path.exists("../results/temp_tr_cds.faa.gz"):
                    os.remove("../results/temp_tr_cds.faa.gz")
                
