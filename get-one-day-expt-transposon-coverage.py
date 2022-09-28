#!/usr/bin/env python

""" get-one-day-expt-transposon-coverage.py by Rohan Maddamsetti.

This script calculates the coverage for the actual B30 and B59 transposons
in the ancestral and evolved samples from my first one-day evolution experiment genome data.

Usage: python get-one-day-expt-transposon-coverage.py > ../results/one-day-expt-genome-analysis/transposon-coverage.csv

"""

import os

## Global constants, from the reference transposon files.
B30_TRANSPOSON_COORDS = (3032, 4409)
B59_TRANSPOSON_COORDS = (1502, 2879)


def get_B59_or_B30_transposon_coverage(transposon_coverage_f):
    """ This function can only handle B30 and B59 coverage files. """
    if "B30" in transposon_coverage_f:
        my_transposon = "B30"
        Tn_start, Tn_end = B30_TRANSPOSON_COORDS
    elif "B59" in transposon_coverage_f:
        my_transposon = "B59"
        Tn_start, Tn_end = B59_TRANSPOSON_COORDS
    else:
        raise AssertionError("Unknown transposon!")

    total_coverage = 0.0
    positions_examined = 0
    with open(transposon_coverage_f, "r") as transposon_coverage_fh:
        ## the header is line number 0! so line 1 is the first line with data.
        for i, line in enumerate(transposon_coverage_fh):
            if (i < Tn_start): continue
            if (i >= Tn_end): break ## don't count the last position of the transposon. some coverage oddity?
            positions_examined += 1
            fields = line.split("\t") ## tab-delimited data.
            top_coverage = float(fields[0])
            bottom_coverage = float(fields[1])
            total_position_coverage = top_coverage + bottom_coverage
            total_coverage += total_position_coverage

        ## round to 2 decimal places.
        my_transposon_coverage = str(round(float(total_coverage)/float(positions_examined),2))
    return (my_transposon, my_transposon_coverage)


def get_transposon_coverage_for_sample(breseq_outpath):
    ## awkward hack to get rid of the "remapped-" prefix for the ancestral clones.
    my_sample = "RM7" + os.path.basename(breseq_outpath).split("RM7")[1]
    coverage_dir = os.path.join(breseq_outpath, "08_mutation_identification")
    transposon_coverage_f = [os.path.join(coverage_dir, x) for x in os.listdir(coverage_dir) if (x.startswith("B59") or x.startswith("B30"))][0]
    my_transposon, my_transposon_coverage = get_B59_or_B30_transposon_coverage(transposon_coverage_f)
    return (my_sample, my_transposon, my_transposon_coverage)


def main():

    breseq_dir = "../results/one-day-expt-genome-analysis"

    ## use the data remapped to the clone GFF3 references.
    ancestral_clone_paths = [os.path.join(breseq_dir, x) for x in os.listdir(breseq_dir) if x.startswith("remapped")]

    evolved_pop_dir = os.path.join(breseq_dir, "one-day-expt-evolved-pops")
    evolved_pop_paths = [os.path.join(evolved_pop_dir, x) for x in os.listdir(evolved_pop_dir) if x.startswith("RM")]


    sample_vec = []
    transposon_vec = []
    transposon_coverage_vec = []
    
    ## walk through the file structure for the ancestral samples.
    for p in ancestral_clone_paths:
        my_sample, my_transposon, my_transposon_coverage = get_transposon_coverage_for_sample(p)
        sample_vec.append(my_sample)
        transposon_vec.append(my_transposon)
        transposon_coverage_vec.append(my_transposon_coverage)
        
    ## walk through the file structure for the evolved samples.
    for p in evolved_pop_paths:
        my_sample, my_transposon, my_transposon_coverage = get_transposon_coverage_for_sample(p)
        sample_vec.append(my_sample)
        transposon_vec.append(my_transposon)
        transposon_coverage_vec.append(my_transposon_coverage)

    ## now print the data to file.
    rownum = len(sample_vec)
    assert len(transposon_vec) == rownum
    assert len(transposon_coverage_vec) == rownum

    print("Sample,Transposon,TransposonCoverage")
    for i in range(rownum):
        myrow = ",".join([sample_vec[i], transposon_vec[i], transposon_coverage_vec[i]])
        print(myrow)

main()
