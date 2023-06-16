#!/usr/bin/env python

'''
fetch-genome-fasta-seqs.py by Rohan Maddamsetti.

This script first reads in ../results/FileS3-Complete-Genomes-in-with-Duplicated-ARG-annotation.csv
to get a list of the Annotation_Accessions that are actually included in the analysis (aka for filtering).

This script then reads in ../results/best-prokaryotes.txt, and downloads the fasta seqs for the
Annotation_Accessions found in FileS3.

NOTE: for path names to be processed properly, this script must be run
from the src/ directory as:  "python fetch-genome-fasta-seqs.py".
'''

import urllib.request
from os.path import basename, exists
import gzip
import os
from tqdm import tqdm


def get_acc_list(FileS3):
    acc_list = []
    with open(FileS3, "r") as FileS3_fh:
        for i, line in enumerate(FileS3_fh):
            if i == 0: continue ## skip the header
            annotation_accession = line.split(',')[0] ## get the first column
            acc_list.append(annotation_accession)
    return(acc_list)


def download_FASTA_genomes(genome_report_file, acc_list=None):
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
            ## If acc_list is not None, then skip genomes that are not in the acc_list.
            if ((acc_list is not None) and (annotation_accession not in acc_list)):
                continue
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


##FileS3="../results/FileS3-Complete-Genomes-in-with-Duplicated-ARG-annotation.csv"
##acc_list = get_acc_list(FileS3)
genome_report_file = "../results/best-prokaryotes.txt"
download_FASTA_genomes(genome_report_file, acc_list=None)
