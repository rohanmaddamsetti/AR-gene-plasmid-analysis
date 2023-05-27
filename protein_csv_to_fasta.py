#!/usr/bin/env python

## protein_csv_to_fasta.py by Rohan Maddamsetti.

## Usage: python protein_csv_to_fasta.py > ../results/duplicate-proteins.faa
## Usage: python protein_csv_to_fasta.py > ../results/Duplicated-ARGs.faa
## Usage: python protein_csv_to_fasta.py > ../results/all-proteins.faa

#input_csv = "../results/duplicate-proteins.csv"
#input_csv = "../results/FileS4-Duplicated-ARGs.csv"
input_csv = "../results/all-proteins.csv"
 
with open(input_csv, "r") as csv_fh:
    for i, line in enumerate(csv_fh):
        if i == 0: continue ## skip the header
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
        print(fasta_header)
        print(seq)
