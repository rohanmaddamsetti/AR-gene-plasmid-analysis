#!/bin/bash
#SBATCH --job-name=test-run-phmmer
#SBATCH --output=../results/test-run-phmmer-hits.out
#SBATCH -e ../results/test-run-phmmer-hits.log
#SBATCH --mem=100G #100 GB RAM

source ~/.bash_profile

phmmer  --noali -E 1e-20 --tblout ../results/test-run-phmmer-hits.out ../data/protein-queries-from-Yi.faa ../results/protein_db.faa

