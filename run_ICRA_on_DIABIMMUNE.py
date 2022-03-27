#===============================================================================
# run_ICRA.py by Rohan Maddamsetti.
# I modified linear_example.py by the SGVFinder authors, for my purposes.
# Run ICRA to take a folder full of fastq files and generate input for SGVFinder.
# NOTE: Memory requirements are at least 20GB per file, and even then - it takes time, depending on the
# sample. So it is highly not recommended to use this file as is, but rather to edit it to work with your
# HPC environment. 
# NOTE2: See the README for requirements etc.
#===============================================================================
from glob import glob
from os.path import join, splitext, basename
from ICRA import single_file
import os

INPUT_FOLDER = "../data/DIABIMMUNE-filtered-reads/E010481"
OUTPUT_FOLDER = "../results/DIABIMMUNE/E010481"
fastqs_1 = glob(join(INPUT_FOLDER, '*.1.fastq'))
for f in fastqs_1: #Parallelize on your cluster
    print f
    ##single_file(f, f.replace('.1.fastq', '.2.fastq'), OUTPUT_FOLDER, 8, True, 1e-6, 100, 10, 100, 100, 60, 1e5, 2e7, 'genomes', False)

