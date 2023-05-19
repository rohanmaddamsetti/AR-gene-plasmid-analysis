#!/usr/bin/env python

"""
process-generality-expt-gdiffs.py by Rohan Maddamsetti.

I take the evolved genomes, and write out evolved-mutations.csv for
downstream analysis in R.

Usage: python process-one-day-expt-gdiffs.py

IMPORTANT NOTE: miniIS1 transpositions are not called properly,
## probably because the boundaries of the repeat don't match the realized increases in coverage.
## if I care about this, (aka goes into publication figures) then I need to re-annotate synthetic MGE boundaries
## and re-run breseq and downstream analysis.

"""

from os.path import join, dirname, realpath
from os import listdir, walk
import pandas as pd

## add genomediff parser to path.
import sys
sys.path.append("genomediff-python")
import genomediff


def write_evolved_pop_mutations(evol_pop_labels, mixed_pop_paths, outf, FALSE_POSITIVE_BLACKLIST):
    outfh = open(outf,"w")
    outfh.write("Sample,Transposon,Plasmid,Antibiotic,AntibioticConcentration,Population,Mutation,Mutation_Category,Gene,Position,Allele,Frequency\n")
    
    for input_dir in mixed_pop_paths:
        for path, subdirs, files in walk(input_dir):
            for f in [x for x in files if x.endswith('annotated.gd')]:
                full_f = join(path, f)
                infh = open(full_f, 'r', encoding='utf-8')
                sample = full_f.split("generality-expts-genome-analysis/")[1].split("/output/evidence/")[0]
                my_row = evol_pop_labels[evol_pop_labels.Sample == sample].iloc[0]
                transposon = my_row["Transposon"]
                plasmid = my_row["Plasmid"]
                antibiotic = my_row["Antibiotic"]
                antibiotic_concentration = str(int(my_row["AntibioticConcentration"]))
                population = str(int(my_row["Population"]))
                

                gd = genomediff.GenomeDiff.read(infh)
                muts = gd.mutations
                muttype = ""
                allele = ""
                for rec in muts:
                    pos = str(rec.attributes['position'])
                    mut_category = rec.attributes['mutation_category']
                    gene = rec.attributes['gene_name']
                    ## skip over false positive blacklisted genes.
                    if gene in FALSE_POSITIVE_BLACKLIST: continue
                    ## annotate transpositions of the mini-transposons.
                    if "repeat_name" in rec.attributes and rec.attributes["repeat_name"].startswith( "mini"):
                        transposon_name = rec.attributes["repeat_name"]
                        gene = "-".join([transposon_name, gene])
                    freq = str(rec.attributes['frequency'])
                    if 'new_seq' in rec.attributes:
                        allele = rec.attributes['ref_seq'] + "->" + rec.attributes['new_seq']
                    else:
                        allele = rec.type
                    if rec.type == 'SNP':
                        muttype = rec.attributes['snp_type']
                    else:
                        muttype = rec.type
                    mutation_row_data = ','.join([sample, transposon, plasmid, antibiotic, antibiotic_concentration, population, muttype, mut_category, gene, pos, allele, freq])
                    outfh.write(mutation_row_data + "\n")


def main():
    srcdir = dirname(realpath(__file__))
    assert srcdir.endswith("src")
    projdir = dirname(srcdir)
    assert projdir.endswith("ARG-duplications")

    genome_results_dir = join(projdir,"results","generality-expts-genome-analysis")

    evolved_pops_label_f = join(projdir,"data","generality-expts-sample-metadata.csv")
    evolved_pop_labels = pd.read_csv(evolved_pops_label_f)

    evolved_pops = [x for x in listdir(genome_results_dir) if x in evolved_pop_labels.Sample.values]
    evolved_pop_paths = [join(genome_results_dir, x, "output", "evidence") for x in evolved_pops]
    
    ''' tabulate the mutations in the evolved populations.'''
    mut_table_outf = join(genome_results_dir,"evolved_mutations.csv")
    FALSE_POSITIVE_BLACKLIST = []
    write_evolved_pop_mutations(evolved_pop_labels, evolved_pop_paths, mut_table_outf, FALSE_POSITIVE_BLACKLIST)

    
main()
    
