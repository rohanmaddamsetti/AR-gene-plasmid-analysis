#!/usr/bin/env python

'''
count-ARG-MGE-adjacencies.py by Rohan Maddamsetti

This script counts the number of duplicated ARGs adjacent to MGE-associated genes,
the number of duplicated ARGs not adjacent to MGE-associated genes,
the number of singleton ARGs adjacent to MGE-associated genes,
and the number of singleton ARGs not adjacent to MGE-associated genes.

Since "tabulate-proteins.py" makes a table of duplicated genes for each genome,
we can easily check whether a given protein is duplicated or not.

The basic idea is to iterate over each gene in each genome, and keep track of its
"left" and "right" neighbors. If the gene is an ARG, then check whether it is a duplicated ARG,
check its neighbors to see whether either one is an MGE-associated gene, and then
update the relevant count.

We need to have special cases to examine the left and right neighbors of the very first
and the very last gene in each replicon.

Pseudocode:

for each genome:
    for each replicon:
        left_neighbor = ""
        right_neighbor = ""

        ## note that the final value of "gene" is the last gene,
        ## and the final value of left_neighbor" is the left neighbor of the last gene.
        ## save info for the two special cases to treat last.
        first_gene = ""
        first_gene_right_neighbor = ""
        for each gene:
            if first_gene == "":
                 first_gene = gene ## save for last
                 left_neighbor = gene
                 continue ## skip, and treat last.
            else if first_gene_right_neighbor is blank
                 ## don't skip in this case, this is the first case for the loop.
                 first_gene_right_neighbor = gene

            if gene is an ARG: ## then examine and update counts.
                      is_duplicated = TRUE if gene is duplicated, else FALSE.
            end

            left_neighbor = gene

'''


import os
from os.path import basename
import re
import gzip
from Bio import SeqIO
from tqdm import tqdm

def make_replication_type_lookup_tbl(infile="../results/chromosome-plasmid-table.csv"):
    replicon_type_lookup_table = {}
    with open(infile, 'r') as chromosome_plasmid_fh:
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
    return replicon_type_lookup_table


def make_duplicated_proteins_lookup_tbl(infile="../results/duplicate-proteins.csv"):
    duplicated_proteins_lookup_table = {}
    with open(infile, 'r') as duplicated_proteins_fh:
        for i, line in enumerate(duplicated_proteins_fh):
            if i == 0: continue ## skip the header
            line = line.strip()
            fields = line.split(',')
            my_annot_accession = fields[0]
            my_product_annot = fields[-2]
            my_sequence = fields[-1]
            if my_annot_accession in duplicated_proteins_lookup_table:
                duplicated_proteins_lookup_table[my_annot_accession][my_sequence] = my_product_annot
            else:
                duplicated_proteins_lookup_table[my_annot_accession] = {my_sequence : my_product_annot}
        return duplicated_proteins_lookup_table


    
## use the data in chromosome-plasmid-table.csv to look up replicon type,
## based on Annotation_Accession and then NCBI Nucleotide ID.
replication_type_lookup_table = make_replicon_type_lookup_tbl()

## populate a dictionary of {Annotation_Accession : {dup_sequence : dup_annotation}},
## from reading in duplicate-proteins.csv.
duplicated_proteins_lookup_table = make_duplicated_proteins_lookup_tbl()


gbk_annotation_dir = "../results/gbk-annotation/"
outf = "../results/ARG-MGE-adjacency-counts.csv"

## the data we care about.
dARGs_next_to_MGEs = 0
dARGs_not_next_to_MGEs = 0
sARGs_next_to_MGEs = 0
sARGs_not_next_to_MGEs = 0


with open(outf, "w") as out_fh:
    header = "dARGs_next_to_MGEs,dARGs_not_next_to_MGEs,sARGs_next_to_MGEs,sARGs_not_next_to_MGEs\n"
    out_fh.write(header)

    gbk_gz_files = [x for x in os.listdir(gbk_annotation_dir) if x.endswith("_genomic.gbff.gz")]
    for gbk_gz in tqdm(gbk_gz_files):        
        infile = gbk_annotation_dir + gbk_gz
        annotation_accession = basename(infile).split("_genomic.gbff.gz")[0]
        ## IMPORTANT TODO: make sure chromosome-plasmid-table.csv
        ## and the data in ../results/gbk-annotation are consistent.
        ## The next line is a temporary consistency check.
        if annotation_accession not in replicon_type_lookup_table: continue
        
        with gzip.open(infile,'rt') as genome_fh:
            for replicon in SeqIO.parse(genome_fh, "gb"):
                replicon_id = replicon.id
                replicon_type = "NA"
                replicon_length = len(replicon.seq)
        
                if replicon_id in replicon_type_lookup_table[annotation_accession]:
                    replicon_type = replicon_type_lookup_table[annotation_accession][replicon_id]
                else: ## replicon is not annotated as a plasmid or chromosome
                    ## in the NCBI Genome report, prokaryotes.txt.
                    ## assume that this is an unassembled contig or scaffold.
                    replicon_type = "contig"

                observed_locations = set() ## check for artifactual duplication due to duplicated annotation.

                if annotation_accession not in duplicated_proteins_lookup_table:
                    ## then there are no duplications in this genome.
                    break
                
                my_dup_dict = duplicated_proteins_lookup_table[annotation_accession]
                hit_first_prot = False
                first_prot = ""
                last_prot = ""
                            
                for feat in replicon.features:
                    if feat.type != "CDS": continue
                    try:
                        prot_seq = feat.qualifiers['translation'][0]
                    except:
                        continue

                    try: ## replace all commas with semicolons for csv formatting.
                        prot_product = feat.qualifiers['product'][0].replace(',',';')
                    except:
                        prot_product = "NA"
                    
                    prot_location = str(feat.location)
                    if prot_location in observed_locations:
                        continue
                    else: ## if not seen before, then add to observed_locations.
                        observed_locations.add(prot_location)

                    prot_id = feat.qualifiers['protein_id'][0]
                    last_prot = prot_id ## this is the last prot we have seen.
                    if not hit_first_prot:
                        my_first_prot = prot_id
                        hit_first_prot = True
                                                
                    cur_prot = { "seq" : prot_seq,
                                 "id" : prot_id,
                                 "product" : prot_product,
                                 "location" : feat.location }

                    ## check if the gene is duplicated.
                    is_duplicated = True if prot_seq in my_dup_dict else False

                

                    ## write out the final counts to file.
                    final_data = [dARGs_next_to_MGEs, dARGs_not_next_to_MGEs,
                                  sARGs_next_to_MGEs, sARGs_not_next_to_MGEs]
                        row = ','.join([str(x) for x in final_data]) + '\n'
                        out_fh.write(row)
                        
