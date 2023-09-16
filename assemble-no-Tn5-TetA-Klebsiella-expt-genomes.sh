#!/usr/bin/env bash

## assemble-no-Tn5-TetA-Klebsiella-expt-genomes.sh by Rohan Maddamsetti.
## The strains in this evolution experiment did not have the B30 miniTn5-TetA transposon in them.
## If I run breseq with the B30-miniTn5-TetA.gb transposon reference sequence, no reads are mapped to it,
## and breseq calls a deletion over the whole sequence.

## This shell script does a bunch of genome assembly tasks using breseq and gdtools.
## COMMENT OUT all tasks that don't need to be done.

################################################################################
## first, assemble the ancestral genomes using the K. oxytoca M5a1 reference genome
## sequenced by Sternberg and Wang labs for the INTEGRATE paper.

#breseq -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1 -r ../data/Klebsiella-expt-reference-genome/Klebsiella_grimontii_M5a1.gb ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_1/*.fastq.gz

#breseq -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2 -r ../data/Klebsiella-expt-reference-genome/Klebsiella_grimontii_M5a1.gb ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_2/*.fastq.gz

################################################################################
## Apply mutations in the ancestral genomes to the K. oxytoca M5a1 reference genome.

#gdtools APPLY -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1.gff3 -f GFF3 -r ../data/Klebsiella-expt-reference-genome/Klebsiella_grimontii_M5a1.gb ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1/output/output.gd

#gdtools APPLY -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2.gff3 -f GFF3 -r ../data/Klebsiella-expt-reference-genome/Klebsiella_grimontii_M5a1.gb ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2/output/output.gd

################################################################################
## now, test the new references, by re-mapping reads.

#breseq -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/remapped-RM7-154-1 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_1/*.fastq.gz

#breseq -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/remapped-RM7-154-2 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_2/*.fastq.gz

################################################################################
## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 30, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.
## Only allow 5 max mismatches between the read and the reference.

######
## RM7.154.3-5 are RM7.154.1 Tet0 evolved pops 1-3.
breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-3 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_3/*.fastq.gz

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-4 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_4/*.fastq.gz 

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-5 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_5/*.fastq.gz

## RM7.154.6-8 are RM7.154.2 Tet0 evolved pops 1-3.

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-6 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_6/*.fastq.gz

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-7 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_7/*.fastq.gz

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-8 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_8/*.fastq.gz

## RM7.154.9-11 are RM7.154.1 Tet5 evolved pops 1-3.

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-9 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_9/*.fastq.gz

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-10 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_10/*.fastq.gz

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-11 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-1.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_11/*.fastq.gz

## RM7.154.12-14 are RM7.154.2 Tet5 evolved pops 1-3.

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-12 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_12/*.fastq.gz

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-13 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_13/*.fastq.gz

breseq -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-14 -r ../results/no-Tn5-TetA-Klebsiella-expt-genome-analysis/RM7-154-2.gff3 ../data/no-Tn5-TetA-Klebsiella-expt-genome-data/SeqCenter-2023-07-29/RM7_154_14/*.fastq.gz

