#!/usr/bin/env python

'''
tabulate-proteins-in-clinical-genomes.py by Rohan Maddamsetti.

Tabulate the proteins in the independent datasets of clinical isolates.
Usage: python tabulate-proteins-in-clinical-genomes.py
'''

import os
import gzip
from Bio import SeqIO
from tqdm import tqdm
import argparse


def tabulate_proteins(inputdir, outf, ignore_singletons=False):
    
    with open(outf, 'w') as outfh:
        header = "Annotation_Accession,count,chromosome_count,plasmid_count,unassembled_count,product,sequence\n"
        outfh.write(header)
        for gbk_gz in tqdm(os.listdir(inputdir)):
            if not gbk_gz.endswith(".gbff.gz"): continue
            annotation_accession = gbk_gz.split("_genomic.gbff.gz")[0]
            infile = inputdir + gbk_gz
            print(infile)
            with gzip.open(infile, "rt") as genome_fh:
                protein_dict = {}
                for replicon in SeqIO.parse(genome_fh, "gb"):
                    replicon_id = replicon.id
                    replicon_type = "contig"
                    if "chromosome" in replicon.description:
                        replicon_type = "chromosome"
                    elif "plasmid" in replicon.description:
                        replicon_type = "plasmid"
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
                if ignore_singletons:
                    ## then throw away all single copy entries.
                    filtered_prot_dict = {k:v for (k,v) in protein_dict.items() if v['count'] > 1}
                else:
                    filtered_prot_dict = protein_dict
                for seq, v in filtered_prot_dict.items():
                    row = ','.join([annotation_accession, str(v["count"]), str(v["chromosome_count"]), str(v["plasmid_count"]), str(v["unassembled_count"]), v["product"], seq])
                    row = row + '\n'
                    outfh.write(row)


def main():

    tabulate_proteins("../data/Mahmud2022-Hybrid-Assemblies-NCBI-BioProject-PRJNA824420/",
                      "../results/Mahmud2022-duplicate-proteins.csv", ignore_singletons=True)
    tabulate_proteins("../data/Mahmud2022-Hybrid-Assemblies-NCBI-BioProject-PRJNA824420/",
                      "../results/Mahmud2022-all-proteins.csv", ignore_singletons=False)

    tabulate_proteins("../data/Hawkey2022-Hybrid-Assemblies-NCBI-BioProject-PRJNA646837/",
                      "../results/Hawkey2022-duplicate-proteins.csv", ignore_singletons=True)
    tabulate_proteins("../data/Hawkey2022-Hybrid-Assemblies-NCBI-BioProject-PRJNA646837/",
                      "../results/Hawkey2022-all-proteins.csv", ignore_singletons=False)

    tabulate_proteins("../data/BARNARDS-LongRead-Assemblies-NCBI-BioProject-PRJNA767644/",
                      "../results/BARNARDS-duplicate-proteins.csv", ignore_singletons=True)
    tabulate_proteins("../data/BARNARDS-LongRead-Assemblies-NCBI-BioProject-PRJNA767644/",
                      "../results/BARNARDS-all-proteins.csv", ignore_singletons=False)

    tabulate_proteins("../data/LongRead-Assemblies-NCBI-BioProject-PRJNA290784/",
                      "../results/Duke-ESBL-duplicate-proteins.csv", ignore_singletons=True)
    tabulate_proteins("../data/LongRead-Assemblies-NCBI-BioProject-PRJNA290784/",
                      "../results/Duke-ESBL-all-proteins.csv", ignore_singletons=False)


main()
