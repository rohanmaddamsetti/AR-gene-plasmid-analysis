#!/bin/bash
#SBATCH --job-name=test-metaspades
#SBATCH --output=../results/meta-spades-E010481.out
#SBATCH --mem=350G #350 GB RAM

spades.py --meta -1 ../data/DIABIMMUNE-filtered-reads/E010481/E010481_33.0.late.trimmo.60.um.1.fastq.gz -2 ../data/DIABIMMUNE-filtered-reads/E010481/E010481_33.0.late.trimmo.60.um.2.fastq.gz -o ../results/meta-spades/E010481
