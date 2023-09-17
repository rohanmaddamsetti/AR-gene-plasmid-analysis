#!/usr/bin/env python

"""
find-bad-replicons.py by Rohan Maddamsetti.

Usage: python find-bad-replicons.py > ../results/bad-replicons.csv


This script prints out bad replicons (chromosomes and plasmids) where some plasmid is
larger than the chromosome, based on the data in
../results/replicon-lengths-and-protein-counts.csv 

"""


def get_replicons_from_genomes_with_chromosomes_smaller_than_a_plasmid(replicon_length_file):
    bad_replicon_tuple_list  = []
    cur_annotation_accession = None
    cur_chromosome_length = -1
    with open(replicon_length_file, "r") as replicon_length_fh:
        for i, line in enumerate(replicon_length_fh):
            if i == 0: continue ## skip the header
            line = line.strip() ## remove leading and lagging whitespace.
            AnnotationAccession, SeqID, SeqType, replicon_length, protein_count = line.split(",")
            replicon_length = int(replicon_length)
            protein_count = int(protein_count)
            ## truncate the NCBI 'NC_' or 'NZ_' prefixes for compatibility with other files.
            NCBI_Nucleotide_Accession = SeqID.replace("NC_", "").replace("NZ_", "")
            my_replicon_tuple = (AnnotationAccession, NCBI_Nucleotide_Accession)
            if AnnotationAccession != cur_annotation_accession and SeqType == "chromosome":
                cur_annotation_accession = AnnotationAccession
                cur_chromosome_length = replicon_length
            else: ## AnnotationAccession == cur_annotation_accession
                if (SeqType == "plasmid") and (replicon_length > cur_chromosome_length) and (my_replicon_tuple not in bad_replicon_tuple_list):
                    bad_replicon_tuple_list.append(my_replicon_tuple)
    return bad_replicon_tuple_list


def main():    
    replicon_length_file = "../results/replicon-lengths-and-protein-counts.csv"
    """ make a list of genomes in which the chromosome is smaller than some plasmid.
    These are often misannotated genomes in which the real chromosome is switched up with a plasmid."""
    bad_replicon_tuple_list = get_replicons_from_genomes_with_chromosomes_smaller_than_a_plasmid(replicon_length_file)

    print("Annotation_Accession,NCBI_Nucleotide_Accession")
    for my_tuple in bad_replicon_tuple_list:
        print(",".join(my_tuple))
    return


if __name__ == "__main__":
    main()
