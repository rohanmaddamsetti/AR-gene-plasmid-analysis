#!/usr/bin/env python

""" get-generality-expts-transposon-coverage.py by Rohan Maddamsetti.

This script calculates the coverage for the actual B30 and B59 transposons
in the ancestral and evolved samples from my first one-day evolution experiment genome data.

Usage: python get-generality-expts-transposon-coverage.py > ../results/generality-expts-genome-analysis/generality-expts-transposon-coverage.csv

"""

import os
import pandas as pd

## Global constants, from the reference transposon files.
B107_TRANSPOSON_COORDS = (2198, 3710)
B111_TET_COORDS = (2607, 3797)
B123_TET_COORDS = (2896, 4086)
B109_TET_COORDS = (3017, 4207)
B110_TRANSPOSON_COORDS = (1440, 4485)
B90_TRANSPOSON_COORDS = (2990, 3969)
B91_TRANSPOSON_COORDS = (2990, 3993)
B92_TRANSPOSON_COORDS = (2990, 4038)
B95_TRANSPOSON_COORDS = (3140, 3987)


def get_generality_transposon_coverage(transposon_coverage_f):
    if "B107" in transposon_coverage_f:
        my_transposon = "B107"
        Tn_start, Tn_end = B107_TRANSPOSON_COORDS
    elif "B111" in transposon_coverage_f:
        my_transposon = "B111"
        Tn_start, Tn_end = B111_TET_COORDS
    elif "B123" in transposon_coverage_f:
        my_transposon = "B123"
        Tn_start, Tn_end = B123_TET_COORDS
    elif "B109" in transposon_coverage_f:
        my_transposon = "B109"
        Tn_start, Tn_end = B109_TET_COORDS
    elif "B110" in transposon_coverage_f:
        my_transposon = "B110"
        Tn_start, Tn_end = B110_TRANSPOSON_COORDS
    elif "B90" in transposon_coverage_f:
        my_transposon = "B90"
        Tn_start, Tn_end = B90_TRANSPOSON_COORDS
    elif "B91" in transposon_coverage_f:
        my_transposon = "B91"
        Tn_start, Tn_end = B91_TRANSPOSON_COORDS
    elif "B92" in transposon_coverage_f:
        my_transposon = "B92"
        Tn_start, Tn_end = B92_TRANSPOSON_COORDS
    elif "B95" in transposon_coverage_f:
        my_transposon = "B95"
        Tn_start, Tn_end = B95_TRANSPOSON_COORDS
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
    transposon_coverage_f = [os.path.join(coverage_dir, x) for x in os.listdir(coverage_dir) if (x.startswith("B107") or x.startswith("B111") or x.startswith("B123") or x.startswith("B109") or x.startswith("B110") or x.startswith("B90") or x.startswith("B91") or x.startswith("B92") or x.startswith("B95"))][0]
    my_transposon, my_transposon_coverage = get_generality_transposon_coverage(transposon_coverage_f)
    return (my_sample, my_transposon, my_transposon_coverage)


def main():

    breseq_dir = "../results/generality-expts-genome-analysis"

    ## get the evolved pops metadata.
    ev_pops = pd.read_csv("../data/generality-expts-sample-metadata.csv")
    
    ## use the data remapped to the clone GFF3 references.
    evolved_pop_paths = [os.path.join(breseq_dir, x) for x in os.listdir(breseq_dir) if x in ev_pops.Sample.values]

    sample_vec = []
    transposon_vec = []
    transposon_coverage_vec = []
            
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
