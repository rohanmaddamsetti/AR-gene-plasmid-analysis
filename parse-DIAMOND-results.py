#!/usr/bin/env python

"""
parse-DIAMOND-results.py by Rohan Maddamsetti

Take the outputs of search-CARD-and-mobileOG-db.py, and make csv
files for analysis in ARG-duplication-analysis.R.

This is a control analysis to show that the results don't specifically depend on genbank annotation.

"""

def ParseHeader(header):
    fields = header.split("|")
    seq_id = fields[0].split("=")[-1]
    annotation_accession = fields[1].split("=")[-1]
    count = fields[2].split("=")[-1]
    chromosome_count = fields[3].split("=")[-1]
    plasmid_count = fields[4].split("=")[-1]
    unassembled_count = fields[5].split("=")[-1]
    ## the join("|") is to take care of the odd case where there is a "|" in the product annotation field!
    product = "|".join(fields[6:]).split('=')[-1].replace("_", " ") 
    return (seq_id, annotation_accession, count, chromosome_count, plasmid_count, unassembled_count, product)


def PrintCSV(results_file, output_file):
    with open(output_file,"w") as output_fh:
        output_fh.write("SeqID,Annotation_Accession,count,chromosome_count,plasmid_count,unassembled_count,product,percent_identity,evalue,bitscore\n")
        with open(results_file, "r") as results_fh:
            for line in results_fh:
                header, percent_identity, evalue, bitscore = line.split("\t")
                seq_id, annotation_accession, count, chromosome_count, plasmid_count, unassembled_count, product = ParseHeader(header)
                row = ",".join([seq_id, annotation_accession, count, chromosome_count, plasmid_count, unassembled_count, product, percent_identity, evalue, bitscore]) + "\n"
                output_fh.write(row)
    return


def main():
    ARG_results_file = "../results/all-proteins-in-CARD.tsv"
    ARG_output_file = "../results/all-proteins-in-CARD.csv"
    MGE_results_file = "../results/all-proteins-in-mobileOG-db.tsv"
    MGE_output_file = "../results/all-proteins-in-mobileOG-db.csv"
    PrintCSV(ARG_results_file, ARG_output_file)
    PrintCSV(MGE_results_file, MGE_output_file)
    return


main()

