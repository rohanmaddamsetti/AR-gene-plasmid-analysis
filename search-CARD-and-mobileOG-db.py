#!/usr/bin/env python

""" search-CARD-and-mobileOG-db.py by Rohan Maddamsetti.

This script uses DIAMOND (https://github.com/bbuchfink/diamond)
to find ARGs (by comparison to the CARD database)
and MGE-associated proteins (by comparison to the mobileOG-db database).

The search parameters are >80% identity over >85% of the target sequence, following
the protocol in the paper "Improved annotation of antibiotic resistance
determinants reveals microbial resistomes cluster by ecology" by Gibson, Forsberg, Dantas
in ISME (2015).
"""

import subprocess
import os


def makeDIAMONDsearchDB(fastapath, db_outfile):
    if not os.path.exists(db_outfile + ".dmnd"):
        makedb_args = ["diamond", "makedb", "--in", fastapath, "-d", db_outfile]
        subprocess.run(makedb_args)        
    return


def runDIAMONDblastp(referenceDB, queryfile, outputfile):
    if not os.path.exists(outputfile):
        DIAMOND_blastp_args = ["diamond", "blastp", "--fast", "--id", "80", "--subject-cover", "85", "--salltitles", "--max-target-seqs", "1", "--outfmt", "6", "qseqid", "pident", "evalue", "bitscore", "-d", referenceDB, "-q", queryfile, "-o", outputfile]
        subprocess.run(DIAMOND_blastp_args)
    return


def main():
    ## make a DIAMOND database for CARD.
    DIAMOND_CARD_db = "../results/CARD"
    makeDIAMONDsearchDB("../data/card-data/protein_fasta_protein_homolog_model.fasta", DIAMOND_CARD_db)
    ## make a DIAMOND database for mobileOG-db.
    DIAMONDmobileOGdb = "../results/mobileOG-db"
    makeDIAMONDsearchDB("../data/mobileOG-db_beatrix-1-6_v1_all/mobileOG-db_beatrix-1.6.All.faa", DIAMONDmobileOGdb)
    
    ## search duplicate-proteins.faa against CARD
    runDIAMONDblastp(DIAMOND_CARD_db, "../results/duplicate-proteins.faa", "../results/duplicate-proteins-in-CARD.tsv")
    ## search duplicate-proteins.faa against mobileOG-db
    runDIAMONDblastp(DIAMONDmobileOGdb, "../results/duplicate-proteins.faa", "../results/duplicate-proteins-in-mobileOG-db.tsv")
    
    ## search all-proteins.faa against CARD 
    runDIAMONDblastp(DIAMOND_CARD_db, "../results/all-proteins.faa", "../results/all-proteins-in-CARD.tsv")
    ## search duplicate-proteins.faa against mobileOG-db
    runDIAMONDblastp(DIAMONDmobileOGdb, "../results/all-proteins.faa", "../results/all-proteins-in-mobileOG-db.tsv")
    return


main()
