#!/usr/bin/env python

'''
tabulate-proteins.py by Rohan Maddamsetti.

This script tabulates all proteins in a strain that have multiple
identical copies across all chromosomes and plasmids.

Usage: python tabulate-proteins.py ## for tabulating all proteins

Usage: python tabulate-proteins.py --ignore-singletons ## for tabulating multi-copy proteins only
'''

import os
import gzip
from Bio import SeqIO
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description="Tabulate protein copies across all chromosomes and plasmids of the bacterial strains in these data.")
parser.add_argument('--ignore-singletons', dest='ignore_singletons', action='store_const',
                    const=True, default=False,
                    help="only tabulate multi-copy proteins (default: tabulate all proteins)")
args = parser.parse_args()
print("ignore singletons? =>", args.ignore_singletons)
if args.ignore_singletons:
    outf = "../results/duplicate-proteins.csv"
else:
    outf = "../results/all-proteins.csv"
print("outfile =>", outf)
    
## use the data in chromosome-plasmid-table.csv to look up replicon type,
## based on Annotation_Accession and then NCBI Nucleotide ID.

replicon_type_lookup_table = {}
chromosome_plasmid_tbl = "../results/chromosome-plasmid-table.csv"
with open(chromosome_plasmid_tbl, 'r') as chromosome_plasmid_fh:
    for i, line in enumerate(chromosome_plasmid_fh):
        if i == 0: continue ## skip the header
        line = line.strip()
        fields = line.split(',')
        my_annot_accession = fields[-1]
        rep_type = fields[-2] 
        rep_id = fields[-3]
        if my_annot_accession in replicon_type_lookup_table:
            replicon_type_lookup_table[my_annot_accession][rep_id] = rep_type
        else:
            replicon_type_lookup_table[my_annot_accession] = {rep_id : rep_type}

with open(outf, 'w') as outfh:
    header = "Annotation_Accession,count,chromosome_count,plasmid_count,unassembled_count,product,sequence\n"
    outfh.write(header)
    for gbk_gz in tqdm(os.listdir("../results/gbk-annotation")):
        if not gbk_gz.endswith(".gbff.gz"): continue
        annotation_accession = gbk_gz.split("_genomic.gbff.gz")[0]
        ## IMPORTANT TODO: make sure chromosome-plasmid-table.csv
        ## and the data in ../results/gbk-annotation are consistent.
        ## The next line is a temporary consistency check.
        if annotation_accession not in replicon_type_lookup_table: continue
        infile = "../results/gbk-annotation/" + gbk_gz
        with gzip.open(infile, "rt") as genome_fh:
            protein_dict = {}
            for replicon in SeqIO.parse(genome_fh, "gb"):
                replicon_id = replicon.id
                replicon_type = "NA"
                if replicon_id in replicon_type_lookup_table[annotation_accession]:
                    replicon_type = replicon_type_lookup_table[annotation_accession][replicon_id]
                else: ## replicon is not annotated as a plasmid or chromosome
                    ## in the NCBI Genome report, prokaryotes.txt.
                    ## assume that this is an unassembled contig or scaffold.
                    replicon_type = "contig"
                for feat in replicon.features:
                    if feat.type != "CDS": continue
                    try:
                        prot_seq = feat.qualifiers['translation'][0]
                    except:
                        continue
                    prot_id = feat.qualifiers['protein_id'][0]
                    try: ## replace all commas with semicolons! otherwise csv format is messed up.
                        prot_product = feat.qualifiers['product'][0].replace(',',';')
                    except:
                        prot_product = "NA"
                    if prot_seq not in protein_dict:
                        protein_dict[prot_seq] = { "count":0,
                                                   "chromosome_count":0,
                                                   "plasmid_count":0,
                                                   "unassembled_count":0,
                                                   "product":prot_product,
                                                   "locations":[]}
                    ## check to see that the location has not been seen before.
                    prot_location = str(feat.location)
                    if prot_location in protein_dict[prot_seq]["locations"]:
                        continue
                    ## in all cases, add the location and update the counts.
                    protein_dict[prot_seq]["locations"].append(prot_location)
                    protein_dict[prot_seq]["count"] += 1
                    if replicon_type == "chromosome": ## keep track of copies on chromosomes and plasmids
                        protein_dict[prot_seq]['chromosome_count'] += 1
                    elif replicon_type == "plasmid":
                        protein_dict[prot_seq]['plasmid_count'] += 1
                    elif replicon_type == "contig":
                        protein_dict[prot_seq]['unassembled_count'] += 1
            if args.ignore_singletons:
                ## then throw away all single copy entries.
                filtered_prot_dict = {k:v for (k,v) in protein_dict.items() if v['count'] > 1}
            else:
                filtered_prot_dict = protein_dict
            for seq, v in filtered_prot_dict.items():
                row = ','.join([annotation_accession,str(v["count"]), str(v["chromosome_count"]), str(v["plasmid_count"]), str(v["unassembled_count"]), v["product"], seq])
                row = row + '\n'
                outfh.write(row)
