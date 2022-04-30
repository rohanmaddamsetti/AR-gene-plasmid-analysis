#!/usr/bin/env python

'''
join-duplications.py by Rohan Maddamsetti

The goal of this script is to find larger regions of duplicated genes
within each genome in my dataset.

The basic idea is to score each gene as being in one of two states:
either within a duplicated region, or outside a duplicated region.

Since "tabulate-proteins.py" makes a table of duplicated genes for each genome,
we can use this information to find contiguous regions of duplicated genes.

Pseudocode:

for each genome:
    list_of_duplications = []
    for each amplicon:
        cur_duplication = []
        for each gene:
            if gene is duplicated:
                append to cur_duplication
            else:
               if cur_duplication is not empty
                   flush/print cur_duplication to list_of_duplications.
                   cur_duplication = []
        if cur_duplication is not empty:
            // this this duplication might wrap around the start of the chromosome.
            if cur_duplication is connected to the first duplication:
                merge cur_duplication to the first duplication.
            else:
                flush/print cur_duplication to list_of_duplications.
'''


import os
from os.path import basename
import re
import gzip
from Bio import SeqIO
from tqdm import tqdm


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

            
## populate a dictionary of {Annotation_Accession : {dup_sequence : dup_annotation}},
## from reading in duplicate-proteins.csv.
duplicated_proteins_lookup_table = {}
duplicated_proteins_tbl = "../results/duplicate-proteins.csv"
with open(duplicated_proteins_tbl, 'r') as duplicated_proteins_fh:
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

            
gbk_annotation_dir = "../results/gbk-annotation/"
outf = "../results/joined-duplicate-proteins.csv"

with open(outf, "w") as out_fh:
    header = "Annotation_Accession,Replicon_Accession,Replicon_type,region_index,region_length,num_proteins_in_region,region_start,region_end,protein_index,protein_id,protein_start,protein_end,protein_length,product,sequence\n"
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
                duplicated_regions = []
                cur_duplication = []
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
                        
                    ## check if the gene is duplicated.
                    if prot_seq in my_dup_dict: ## if so, add to the current dup.
                        cur_prot = { "seq" : prot_seq,
                                     "id" : prot_id,
                                     "product" : prot_product,
                                     "location" : feat.location }
                        cur_duplication.append(cur_prot)
                    else:
                        if len(cur_duplication):
                            duplicated_regions.append(cur_duplication)
                            cur_duplication = []
                if len(cur_duplication):
                    if len(duplicated_regions) == 0: ## then this is the only duplication.
                        duplicated_regions.append(cur_duplication)
                    else:
                        first_duplication = duplicated_regions[0]
                        first_dup_prot = first_duplication[0]["id"]
                        last_dup_prot = cur_duplication[-1]["id"]
                        if (first_dup_prot == first_prot) and (last_dup_prot == last_prot):
                            ## then merge these duplicated regions.
                            wrapped_duplication = cur_duplication + first_duplication
                            duplicated_regions[0] = wrapped_duplication
                        else:
                            duplicated_regions.append(cur_duplication)
                
                ## write out the duplicated regions in this replicon.
                ## use 1-based indexing for R.
                for my_dup_region_index, dup_region in enumerate(duplicated_regions,start=1):
                    ## IMPORTANT NOTE: unfortunately, CompoundLocations that span the
                    ## beginning/end of the replicon will have a start of 0 and end
                    ## that is the size of the replicon. The length should be correct though.
                    my_region_start = dup_region[0]["location"].start
                    my_region_end = dup_region[-1]["location"].end
                    ## for simplicity, define the length of the region as
                    ## the sum of the lengths of the duplicated proteins.
                    my_region_length = sum([len(prot["location"]) for prot in dup_region])
                    for my_dup_index, dup in enumerate(dup_region, start=1):
                        row_data = [annotation_accession,replicon_id,replicon_type,
                                    my_dup_region_index, my_region_length,
                                    len(dup_region), ## number of proteins in the region.
                                    my_region_start, my_region_end,
                                    my_dup_index,dup["id"],dup["location"].start,
                                    dup["location"].end, len(dup["location"]),
                                    dup["product"],dup["seq"]]
                        row = ','.join([str(x) for x in row_data]) + '\n'
                        out_fh.write(row)

