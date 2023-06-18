#!/bin/bash

## fetch-and-dereplicate-seqs.sh by Rohan Maddamsetti

mkdir ../results/genome-fasta-files  
mkdir ../results/dereplicated-genome-fasta-files  
## fetch FASTA sequences for the genomes
python fetch-genome-fasta-seqs.py  

## shell command used to run Assembly-deplicator.
./external/Assembly-dereplicator/dereplicator.py ../results/genome-fasta-files ../results/dereplicated-genome-fasta-files --distance 0.005  

## now write the dereplicated strains to file.
ls ../results/dereplicated-genome-fasta-files | grep "_genomic.fna.gz" > ../results/dereplicated-genomes.txt  
