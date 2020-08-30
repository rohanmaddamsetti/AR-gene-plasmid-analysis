#!/bin/bash
#SBATCH --job-name=test-hmmsearch
#SBATCH --output=../results/Lactamase_B_hits.out
#SBATCH -e ../results/hmmsearch-Lactamase_B.log
#SBATCH --mem=100G #100 GB RAM

source ~/.bash_profile

hmmsearch ../data/Lactamase_B.hmm ../results/protein_db.fa

