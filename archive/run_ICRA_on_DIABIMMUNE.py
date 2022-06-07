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

INPUT_FOLDER = "/work/rm431/SGVFinder/data/DIABIMMUNE-filtered-reads"
OUTPUT_FOLDER = "/work/rm431/SGVFinder/results/DIABIMMUNE"
## make the folder if it doesn't exist.
## see this stackoverflow answer-- this code avoids a possible race condition.
## https://stackoverflow.com/questions/273192/how-can-i-safely-create-a-nested-directory
try:
    os.makedirs(OUTPUT_FOLDER)
except OSError:
    if not os.path.isdir(OUTPUT_FOLDER):
        raise

fastqs_1 = glob(join(INPUT_FOLDER, '*.1.fastq'))
idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
## This trick depends on running this script using a slurm jobarray.
f = fastqs_1[idx]
single_file(f, f.replace('.1.fastq', '.2.fastq'), OUTPUT_FOLDER, 8, True, 1e-6, 100, 10, 100, 100, 60, 1e5, 2e7, 'genomes', False)

