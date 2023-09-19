#!/usr/bin/env python

## protein_csv_to_fasta.py by Rohan Maddamsetti.

## Usage: python protein_csv_to_fasta.py
## Usage: python protein_csv_to_fasta.py --ignore-singletons

import argparse


parser = argparse.ArgumentParser(description="make FASTA query files of proteins for CARD and mobileOG-db searches.")
parser.add_argument('--ignore-singletons', dest='ignore_singletons', action='store_const',
                    const=True, default=False,
                    help="only work with multi-copy proteins (default: tabulate all proteins)")
args = parser.parse_args()
print("ignore singletons? =>", args.ignore_singletons)
if args.ignore_singletons:
    input_csv = "../results/duplicate-proteins.csv"
    outf = "../results/duplicate-proteins.faa"
else:
    input_csv = "../results/all-proteins.csv"
    outf = "../results/all-proteins.faa"
print("infile =>", input_csv)
print("outfile =>", outf)


with open(outf, "w") as outfh:
    with open(input_csv, "r") as csv_fh:
        for i, line in enumerate(csv_fh):
            if i == 0: continue ## skip the header
            line = line.strip() ## remove leading and lagging whitespace.
            fields = line.split(",")
            seq_id = fields[0]
            annotation_accession = fields[1]
            count = fields[2]
            chromosome_count = fields[3]
            plasmid_count = fields[4]
            unassembled_count = fields[5]
            product = fields[6]
            seq = fields[7]
            fasta_header = ">" + "|".join(["SeqID="+seq_id, "annotation_accession="+annotation_accession, "count="+count, "chromosome_count="+chromosome_count, "plasmid_count="+plasmid_count, "unassembled_count="+unassembled_count, "product="+product]).replace(" ","_")
            outfh.write(fasta_header + "\n")
            outfh.write(seq + "\n")
