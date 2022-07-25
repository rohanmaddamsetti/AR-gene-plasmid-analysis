#!/usr/bin/env python

''' count-ARG-MGE-adjacencies.py by Rohan Maddamsetti

IMPORTANT NOTE: This script runs on ALL downloaded genomes, and not just those
which are annotated in a particular category. I should figure out whether this is what we
want, or whether I could restrict to the annotated genomes for consistency with the rest
of the data analysis.

TODO: adjust this behavior, so that only annotated genomes are counted.

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
'''

import os
from os.path import basename
import re
import gzip
from itertools import islice
from Bio import SeqIO
from tqdm import tqdm


## This is from the iterools documentation here:
## https://docs.python.org/release/2.3.5/lib/itertools-example.html
## the recipe still works in python3.
def window(seq, n=3):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result


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


def get_prot_data(feature):
    try:
        prot_seq = feature.qualifiers['translation'][0]
    except:
        prot_seq = "NA"
        
    try: ## replace all commas with semicolons for csv formatting.
        prot_product = feature.qualifiers['product'][0].replace(',',';')
    except:
        prot_product = "NA"
        
    prot_location = str(feature.location)
    prot_id = feature.qualifiers['protein_id'][0]
                        
    cur_prot = { "seq" : prot_seq,
                 "id" : prot_id,
                 "product" : prot_product,
                 "location" : prot_location }
    return cur_prot


def get_ARG_adjacency_type(left_prot, cur_prot, right_prot, dup_dict):
    '''
    return 0 if the protein is not an ARG.
    return 1 if the protein is a duplicated ARG next to an MGE.
    return 2 if the protein is a duplicated ARG not next to an MGE.
    return 3 if the protein is a singleton ARG next to an MGE.
    return 4 if the protein is a singleton ARG not next to an MGE.
    '''

    ARG_regex = "chloramphenicol|Chloramphenicol|tetracycline|Tetracycline|macrolide|lincosamide|streptogramin|multidrug|lactamase|glycopeptide resistance|VanZ|bacitracin|polymyxin B|trimethoprim-resistant|sulfonamide-resistant|quinolone|Quinolone|oxacin|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythomycin|antibiotic resistance"

    MGE_regex = "IS|transposon|Transposase|transposase|Transposable|transposable|virus|Phage|phage|integrase|Integrase|baseplate|tail|intron|Mobile|mobile|antitoxin|toxin|capsid|plasmid|Plasmid|conjug|Tra"
    
    ## check if the gene is an ARG.
    is_ARG = True if re.search(ARG_regex, cur_prot["product"]) else False
    if not is_ARG:
        return 0 ## skip if we're not looking at an ARG.
    
    ## check if any neighbors are MGE proteins.
    left_neighbor_is_MGE = True if re.search(MGE_regex, left_prot["product"]) else False
    right_neighbor_is_MGE = True if re.search(MGE_regex, right_prot["product"]) else False
    neighbor_is_MGE = True if left_neighbor_is_MGE or right_neighbor_is_MGE else False
    
    ## check if the ARG is duplicated.
    is_duplicated = True if cur_prot["seq"] in dup_dict else False
    
    ## now update the counters.
    if is_duplicated and neighbor_is_MGE:
        adj_type = 1
        print("D-ARG AND NEIGHBOR IS MGE")
    elif is_duplicated and not neighbor_is_MGE:
        adj_type = 2
        print("D-ARG AND NO MGE NEIGHBOR")
    elif not is_duplicated and neighbor_is_MGE:
        adj_type = 3
        print("S-ARG AND NEIGHBOR IS MGE")
    elif not is_duplicated and not neighbor_is_MGE:
        adj_type = 4
        print("S-ARG AND NO MGE NEIGHBOR")
    else:
        raise AssertionError("this line should never run.")

    ## This is for debugging-- make sure we are getting the right output.
    print(left_prot["product"])
    print(cur_prot["product"])
    print(right_prot["product"])

    return adj_type

    
def count_ARG_MGE_adjacencies(gbk_annotation_dir, duplicated_proteins_lookup_table):

    ## the data we care about.
    dARGs_next_to_MGEs = 0
    dARGs_not_next_to_MGEs = 0
    sARGs_next_to_MGEs = 0
    sARGs_not_next_to_MGEs = 0

    
    gbk_gz_files = [x for x in os.listdir(gbk_annotation_dir) if x.endswith("_genomic.gbff.gz")]
    for gbk_gz in tqdm(gbk_gz_files):        
        infile = gbk_annotation_dir + gbk_gz
        annotation_accession = basename(infile).split("_genomic.gbff.gz")[0]
        
        with gzip.open(infile,'rt') as genome_fh:
            for replicon in SeqIO.parse(genome_fh, "gb"):
                ## skip if there are no duplications in this genome.
                if annotation_accession not in duplicated_proteins_lookup_table:
                    break
                
                observed_locations = set() ## check for artifactual duplication due to duplicated annotation.
                dup_dict = duplicated_proteins_lookup_table[annotation_accession]
                
                first_prot = dict()
                first_prot_right_neighbor = dict()
                ## filter the features for proteins.
                replicon_proteins = (feature for feature in replicon.features if feature.type == "CDS")
                
                ## iterate over windows of 3 proteins.
                for (left_neighbor, feat, right_neighbor) in window(replicon_proteins, n=3):

                    left_prot = get_prot_data(left_neighbor)
                    cur_prot = get_prot_data(feat)
                    right_prot = get_prot_data(right_neighbor)

                    ## initialize the first_prot data. we will handle this after the loop.
                    if first_prot == "":
                        first_prot = left_prot
                        first_prot_right_neighbor = cur_prot
                        
                    if cur_prot["location"] in observed_locations:
                        continue
                    else: ## if not seen before, then add to observed_locations.
                        observed_locations.add(cur_prot["location"])

                    ## figure out what counter to update.
                    adj_type = get_ARG_adjacency_type(left_prot, cur_prot, right_prot, dup_dict)
                    if adj_type == 0:
                        continue
                    elif adj_type == 1:
                        dARGs_next_to_MGEs += 1
                    elif adj_type == 2:
                        dARGs_not_next_to_MGEs += 1
                    elif adj_type == 3:
                        sARGs_next_to_MGEs += 1
                    elif adj_type == 4:
                        sARGs_not_next_to_MGEs += 1
                    else:
                        raise AssertionError("this line shouldn't run too.")       
                ## now we have finished iterating through the replicon,
                ## handle the first and last protein cases.
                last_prot_left_neighbor = cur_prot
                last_prot = right_prot
                last_prot_right_neighbor = first_prot

                last_adj_type = get_ARG_adjacency_type(last_prot_left_neighbor, last_prot, last_prot_right_neighbor, dup_dict)
                if last_adj_type == 0:
                    break
                elif last_adj_type == 1:
                    dARGs_next_to_MGEs += 1
                elif last_adj_type == 2:
                    dARGs_not_next_to_MGEs += 1
                elif last_adj_type == 3:
                    sARGs_next_to_MGEs += 1
                elif last_adj_type == 4:
                    sARGs_not_next_to_MGEs += 1
                else:
                    raise AssertionError("this line shouldn't run too.")
    
                first_prot_left_neighbor = last_prot
                
                first_adj_type = get_ARG_adjacency_type(first_prot_left_neighbor, first_prot, first_prot_right_neighbor, dup_dict)
                if first_adj_type == 0:
                    break
                elif first_adj_type == 1:
                    dARGs_next_to_MGEs += 1
                elif first_adj_type == 2:
                    dARGs_not_next_to_MGEs += 1
                elif first_adj_type == 3:
                    sARGs_next_to_MGEs += 1
                elif first_adj_type == 4:
                    sARGs_not_next_to_MGEs += 1
                else:
                    raise AssertionError("this line shouldn't run three")
                
    return (dARGs_next_to_MGEs, dARGs_not_next_to_MGEs,
            sARGs_next_to_MGEs, sARGs_not_next_to_MGEs)


def main():
    ## populate a dictionary of {Annotation_Accession : {dup_sequence : dup_annotation}},
    ## from reading in duplicate-proteins.csv.
    duplicated_proteins_lookup_table = make_duplicated_proteins_lookup_tbl()

    gbk_annotation_dir = "../results/gbk-annotation/"
    outf = "../results/ARG-MGE-adjacency-counts.csv"

    with open(outf, "w") as out_fh:
        header = "dARGs_next_to_MGEs,dARGs_not_next_to_MGEs,sARGs_next_to_MGEs,sARGs_not_next_to_MGEs\n"
        out_fh.write(header)
        ## write out the final counts to file.
        final_data = count_ARG_MGE_adjacencies(gbk_annotation_dir, duplicated_proteins_lookup_table)
        row = ','.join([str(x) for x in final_data]) + '\n'
        out_fh.write(row)

## run the script.
main()
