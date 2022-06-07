#!/usr/bin/env python

'''
tabulate-Duke-ESBL-proteins.py by Rohan Maddamsetti.

This script tabulates all proteins in a strain that have multiple
identical copies across all chromosomes and plasmids.

Usage: python tabulate-Duke-ESBL-proteins.py ## for tabulating all proteins

Usage: python tabulate-Duke-ESBL-proteins.py --ignore-singletons ## for tabulating multi-copy proteins only
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
    outf = "../results/Duke-ESBL-duplicate-proteins.csv"
else:
    outf = "../results/Duke-ESBL-all-proteins.csv"
print("outfile =>", outf)
    

inputdir = "../data/LongRead-Assemblies-NCBI-BioProject-PRJNA290784/"

with open(outf, 'w') as outfh:
    header = "Annotation_Accession,count,product,sequence\n"
    outfh.write(header)
    for gbk_gz in tqdm(os.listdir(inputdir)):
        if not gbk_gz.endswith(".gbff.gz"): continue
        annotation_accession = gbk_gz.split("_genomic.gbff.gz")[0]
        infile = inputdir + gbk_gz
        print(infile)
        with gzip.open(infile, "rt") as genome_fh:
            protein_dict = {}
            for replicon in SeqIO.parse(genome_fh, "gb"):
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
                                                   "product":prot_product,
                                                   "locations":[]}
                    ## check to see that the location has not been seen before.
                    prot_location = str(feat.location)
                    if prot_location in protein_dict[prot_seq]["locations"]:
                        continue
                    ## in all cases, add the location and update the counts.
                    protein_dict[prot_seq]["locations"].append(prot_location)
                    protein_dict[prot_seq]["count"] += 1
            if args.ignore_singletons:
                ## then throw away all single copy entries.
                filtered_prot_dict = {k:v for (k,v) in protein_dict.items() if v['count'] > 1}
            else:
                filtered_prot_dict = protein_dict
            for seq, v in filtered_prot_dict.items():
                row = ','.join([annotation_accession,str(v["count"]), v["product"], seq])
                row = row + '\n'
                outfh.write(row)
