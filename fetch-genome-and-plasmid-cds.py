#!/usr/bin/env python

'''
fetch-genome-and-plasmid-cds.py by Rohan Maddamsetti.

This script reads in ../results/prokaryotes-with-plasmids.txt.

NOTE: for path names to be processed properly, this script must be run
from the src/ directory as python fetch-genome-and-plasmid-cds.py.

'''

import urllib.request
from os.path import basename
import gzip
import os
from tqdm import tqdm

## create an empty gene database file.
with open("../results/gene_db.fna", "w") as gene_db_fh:
    ## create an empty protein database file.
    with open("../results/protein_db.faa", "w") as protein_db_fh:
        ## open the genome report file, and parse line by line.
        with open("../results/prokaryotes-with-plasmids.txt", "r") as genome_report_fh:
            for i, line in enumerate(tqdm(genome_report_fh)):
                line = line.strip()
                if i == 0: ## get the names of the columns from the header.
                    column_names_list = line.split('\t')
                    continue ## don't process the header further.
                fields = line.split('\t')
                ftp_path = fields[20]
                ## NOW GET PROTEINS.
                ## download the gzipped translated cds (amino acid sequences).
                tr_cds_ftp_path = ftp_path + '/' + basename(ftp_path) + "_translated_cds.faa.gz"
                try:
                    urllib.request.urlretrieve(tr_cds_ftp_path,filename="../results/temp_tr_cds.faa.gz")
                except urllib.error.URLError: ## some problem happens
                    try: ## try one more time, in case there was some connection problem.
                        urllib.request.urlretrieve(tr_cds_ftp_path,filename="../results/temp_tr_cds.faa.gz")
                    except: ## skip if the translated cds don't exist.
                        continue
                ## unzip, concatenate to the protein database file, and delete the temporary file.
                with gzip.open("../results/temp_tr_cds.faa.gz", "rt") as temp_tr_cds_infile:
                    protein_db_fh.write(temp_tr_cds_infile.read())
                    os.remove("../results/temp_tr_cds.faa.gz")
                ## NOW GET GENES.
                ## download the cds (nucleotide sequences).
                cds_ftp_path = ftp_path + '/' + basename(ftp_path) + "_cds_from_genomic.fna.gz"
                try:
                    urllib.request.urlretrieve(cds_ftp_path,filename="../results/temp_cds.fna.gz")
                except urllib.error.URLError: ## some problem happens
                    try: ## try one more time, in case there was some connection problem.
                        urllib.request.urlretrieve(cds_ftp_path,filename="../results/temp_cds.fna.gz")
                    except: ## skip if the cds don't exist.
                        continue
                ## unzip, concatenate to the gene database file, and delete the temporary file.
                with gzip.open("../results/temp_cds.fna.gz", "rt") as temp_cds_infile:
                    gene_db_fh.write(temp_cds_infile.read())
                    os.remove("../results/temp_cds.fna.gz")

