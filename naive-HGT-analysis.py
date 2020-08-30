#!/usr/bin/env python

'''
naive-HGT-analysis.py by Rohan Maddamsetti.

This script counts and tabulates all proteins in a strain that have multiple
identical copies across all chromosomes and plasmids.

Usage: python naive-HGT-analysis.py
'''

import os
import gzip
from Bio import SeqIO

outf = "../results/naive-HGT.csv"
with open(outf, 'w') as outfh:
    header = "Annotation_Accession,count,chromosome_count,plasmid_count,product,sequence\n"
    outfh.write(header)
    print(header)
    for gbk_gz in os.listdir("../results/gbk-annotation"):
        if not gbk_gz.endswith(".gbff.gz"): continue
        annotation_accession = gbk_gz.split("_genomic.gbff.gz")[0]
        infile = "../results/gbk-annotation/" + gbk_gz
        with gzip.open(infile, "rt") as genome_fh:
            protein_dict = {}
            for replicon in SeqIO.parse(genome_fh, "gb"):
                replicon_id = replicon.id
                replicon_type = "NA"
                replicon_description = replicon.description
                if "chromosome" in replicon.description:
                    replicon_type = "chromosome"
                elif "plasmid" in replicon.description:
                    replicon_type = "plasmid"
                ''' note: I can get all the nice genbank annotations from the replicon.annotations object.
                    leave alone for now.'''
                ##print(replicon.annotations)
                for feat in replicon.features:
                    if feat.type != "CDS": continue
                    try:
                        prot_seq = feat.qualifiers['translation'][0]
                    except: continue
                    prot_id = feat.qualifiers['protein_id'][0]
                    try:
                        prot_product = feat.qualifiers['product'][0]
                    except:
                        prot_product = "NA"
                    if prot_seq not in protein_dict:
                        protein_dict[prot_seq] = { "count":0,
                                                   "chromosome_count":0,
                                                   "plasmid_count":0,
                                                   "product":prot_product}
                    ## in all cases, update the counts.
                    protein_dict[prot_seq]["count"] += 1
                    if replicon_type == "chromosome": ## keep track of copies on chromosomes and plasmids
                        protein_dict[prot_seq]['chromosome_count'] += 1
                    elif replicon_type == "plasmid":
                        protein_dict[prot_seq]['plasmid_count'] += 1
            ## now throw away all single copy entries.
            filtered_prot_dict = {k:v for (k,v) in protein_dict.items() if v['count'] > 1}
            for seq, v in filtered_prot_dict.items():
                row = ','.join([annotation_accession,str(v["count"]), str(v["chromosome_count"]), str(v["plasmid_count"]), v["product"], seq])
                row = row + '\n'
                outfh.write(row)
                print(row)
