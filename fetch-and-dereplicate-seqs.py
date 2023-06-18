#!/usr/bin/env python

## fetch-and-dereplicate-seqs.py by Rohan Maddamsetti

import os
import subprocess

subprocess.call(["mkdir", "../results/genome-fasta-files"])
subprocess.call(["mkdir", "../results/dereplicated-genome-fasta-files"])
## fetch FASTA sequences for the genomes
subprocess.call(["python", "fetch-genome-fasta-seqs.py"])

## shell command used to run Assembly-deplicator.
subprocess.call(["./external/Assembly-dereplicator/dereplicator.py", "../results/genome-fasta-files", "../results/dereplicated-genome-fasta-files","--distance", "0.005"])

## now write the dereplicated strains to file.
subprocess.call("ls ../results/dereplicated-genome-fasta-files | grep \"_genomic.fna.gz\" > ../results/dereplicated-genomes.txt", shell=True)

