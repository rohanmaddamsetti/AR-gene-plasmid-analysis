#!/bin/bash
#SBATCH --job-name=test-metaplasmidspades
#SBATCH --output=../results/meta-plasmid-spades-E010481.out
#SBATCH --mem=300G #300 GB RAM

source ~/.bash_profile

spades.py --meta --plasmid -1 ../data/DIABIMMUNE-filtered-reads/E010481/E010481_33.0.late.trimmo.60.um.1.fastq.gz -2 ../data/DIABIMMUNE-filtered-reads/E010481/E010481_33.0.late.trimmo.60.um.2.fastq.gz -o ../results/meta-plasmid-spades/E010481
