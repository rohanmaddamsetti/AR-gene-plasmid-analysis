#!/usr/bin/env python

""" count-ARG-MGE-adjacencies.py by Rohan Maddamsetti

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

There are two special cases to examine the left and right neighbors of the very first
and the very last gene in each replicon.

Usage on DCC: sbatch -p scavenger -t 24:00:00 --mem=4G --wrap="python count-ARG-MGE-adjacencies.py"
"""

import os
from os.path import basename
import re
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

        
def make_replicon_type_lookup_tbl(infile="../results/chromosome-plasmid-table.csv"):
    ## use the data in chromosome-plasmid-table.csv to look up replicon type,
    ## based on Annotation_Accession and then NCBI Nucleotide ID.
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
            my_annot_accession = fields[1]
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
                        
    cur_prot = { "seq" : prot_seq,
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

    chloramphenicol_keywords = "chloramphenicol|Chloramphenicol"
    tetracycline_keywords = "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
    MLS_keywords = "macrolide|lincosamide|streptogramin"
    multidrug_keywords = "Multidrug resistance|multidrug resistance|antibiotic resistance"
    beta_lactam_keywords = "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\S*"
    glycopeptide_keywords = "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
    polypeptide_keywords = "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
    diaminopyrimidine_keywords = "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
    sulfonamide_keywords = "sulfonamide|Sul1|sul1|sulphonamide"
    quinolone_keywords = "quinolone|Quinolone|oxacin|qnr|Qnr"
    aminoglycoside_keywords = "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
    macrolide_keywords = "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
    antimicrobial_keywords = "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\S*"
    ARG_regex = "|".join([chloramphenicol_keywords, tetracycline_keywords,
                          MLS_keywords, multidrug_keywords, beta_lactam_keywords,
                          glycopeptide_keywords, polypeptide_keywords, diaminopyrimidine_keywords,
                          sulfonamide_keywords, quinolone_keywords, aminoglycoside_keywords,
                          macrolide_keywords, antimicrobial_keywords])

    
    transposon_keywords = "IS|transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|tra[A-Z]|conjugate transposon|Transpos\S*|Tn[0-9]|tranposase|Tnp|Ins|ins"
    plasmid_keywords = "relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation|Mob\S*|Plasmid|Rep|Conjug\S*"
    phage_keywords = "capsid|phage|Tail|tail|head|tape measure|antiterminatio|Phage|virus|Baseplate|baseplate|coat|entry exclusion"
    other_HGT_keywords = "Integrase|integrase|excision\S*|exonuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip|intron|antitoxin|toxin|Toxin|Reverse transcriptase|hok|Hok|competence|addiction"
    MGE_regex = "|".join([transposon_keywords, plasmid_keywords, phage_keywords, other_HGT_keywords])
    
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

    
def count_ARG_MGE_adjacencies(gbk_annotation_dir, duplicated_proteins_lookup_table, replicon_type_lookup_table):

    ## the data we care about.
    dARGs_next_to_MGEs = 0
    dARGs_not_next_to_MGEs = 0
    sARGs_next_to_MGEs = 0
    sARGs_not_next_to_MGEs = 0
    
    gbk_files = [x for x in os.listdir(gbk_annotation_dir) if x.endswith("_genomic.gbff")]
    for gbk in tqdm(gbk_files):        
        infile = os.path.join(gbk_annotation_dir, gbk)
        annotation_accession = basename(infile).split("_genomic.gbff")[0]
        ## skip if there are no duplications in this genome.
        if annotation_accession not in duplicated_proteins_lookup_table:
            continue
        ## only examine genomes listed in chromosome-plasmid-table.csv
        ## for consistency in the data analysis.
        if annotation_accession not in replicon_type_lookup_table:
            continue
        
        with open(infile,'r') as genome_fh:
            for replicon in SeqIO.parse(genome_fh, "gb"):
                observed_locations = set() ## check for artifactual duplication due to duplicated annotation.
                dup_dict = duplicated_proteins_lookup_table[annotation_accession]
                
                first_prot = dict()
                first_prot_right_neighbor = dict()
                ## filter the features for proteins.
                replicon_proteins = (feature for feature in replicon.features if feature.type == "CDS")

                loop_ran = False
                ## iterate over windows of 3 proteins.
                for (left_neighbor, feat, right_neighbor) in window(replicon_proteins, n=3):
                    loop_ran = True
                    
                    left_prot = get_prot_data(left_neighbor)
                    cur_prot = get_prot_data(feat)
                    right_prot = get_prot_data(right_neighbor)

                    ## initialize the first_prot data. we will handle this after the loop.
                    if not first_prot: ## the dictionary is empty.
                        first_prot = left_prot
                        first_prot_right_neighbor = cur_prot
                        
                    if cur_prot["location"] in observed_locations:
                        continue
                    else: ## if not seen before, then add to observed_locations.
                        observed_locations.add(cur_prot["location"])

                    ## figure out what counter to update.
                    adj_type = get_ARG_adjacency_type(left_prot, cur_prot, right_prot, dup_dict)
                    if adj_type == 1:
                        dARGs_next_to_MGEs += 1
                    elif adj_type == 2:
                        dARGs_not_next_to_MGEs += 1
                    elif adj_type == 3:
                        sARGs_next_to_MGEs += 1
                    elif adj_type == 4:
                        sARGs_not_next_to_MGEs += 1
                        
                if (loop_ran):
                ## now we have finished iterating through the replicon,
                    ## handle the first and last protein cases.
                    last_prot_left_neighbor = cur_prot
                    last_prot = right_prot
                    last_prot_right_neighbor = first_prot

                    last_adj_type = get_ARG_adjacency_type(last_prot_left_neighbor, last_prot, last_prot_right_neighbor, dup_dict)
                    if last_adj_type == 1:
                        dARGs_next_to_MGEs += 1
                    elif last_adj_type == 2:
                        dARGs_not_next_to_MGEs += 1
                    elif last_adj_type == 3:
                        sARGs_next_to_MGEs += 1
                    elif last_adj_type == 4:
                        sARGs_not_next_to_MGEs += 1

                    first_prot_left_neighbor = last_prot

                    first_adj_type = get_ARG_adjacency_type(first_prot_left_neighbor, first_prot, first_prot_right_neighbor, dup_dict)
                    if first_adj_type == 1:
                        dARGs_next_to_MGEs += 1
                    elif first_adj_type == 2:
                        dARGs_not_next_to_MGEs += 1
                    elif first_adj_type == 3:
                        sARGs_next_to_MGEs += 1
                    elif first_adj_type == 4:
                        sARGs_not_next_to_MGEs += 1

    return (dARGs_next_to_MGEs, dARGs_not_next_to_MGEs,
            sARGs_next_to_MGEs, sARGs_not_next_to_MGEs)


def main():
    ## populate a dictionary of {Annotation_Accession : {dup_sequence : dup_annotation}},
    ## from reading in duplicate-proteins.csv.
    duplicated_proteins_lookup_table = make_duplicated_proteins_lookup_tbl()
    replicon_type_lookup_table = make_replicon_type_lookup_tbl()
    
    gbk_annotation_dir = "../results/gbk-annotation/"
    outf = "../results/ARG-MGE-adjacency-counts.csv"

    with open(outf, "w") as out_fh:
        header = "dARGs_next_to_MGEs,dARGs_not_next_to_MGEs,sARGs_next_to_MGEs,sARGs_not_next_to_MGEs\n"
        out_fh.write(header)
        ## write out the final counts to file.
        final_data = count_ARG_MGE_adjacencies(gbk_annotation_dir, duplicated_proteins_lookup_table, replicon_type_lookup_table)
        row = ','.join([str(x) for x in final_data]) + '\n'
        out_fh.write(row)

## run the script.
main()
