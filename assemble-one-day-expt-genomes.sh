#!/usr/bin/env bash

## assemble-one-day-expt-genomes.sh by Rohan Maddamsetti.
## This shell script does a bunch of genome assembly tasks using breseq and gdtools.
## COMMENT OUT all tasks that don't need to be done.

## IMPORTANT TODO: there may be a bug in which the synthetic transposon is annotated as " Synthetic Tn5" instead of "Synthetic Tn5"
## in the GFF3 files made by gdtools. This causes breseq 0.37 to fail. I fixed this my deleting the space by hand in my GFF files, but
## this bug may reoccur when automatically generating the GFF3 files.

################################################################################
## first, assemble the ancestral genomes using the K12-MG1655-NC_000913.gb reference genome.

## assemble the B59 ancestral strain using the K12-MG1655-NC_000913.gb reference genome.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/RM7-87-1 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B59-TetA.gb ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_1/*.fastq.gz"

## assemble the B59 + p15A ancestral strain using the K12-MG1655-NC_000913.gb reference genome.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/RM7-87-2 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B59-TetA.gb -r ../data/one-day-expt-reference-genome/A31-p15A.gb ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_2/*.fastq.gz"

## assemble the B59+ pUC ancestral strain using the K12-MG1655-NC_000913.gb reference genome.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/RM7-87-3 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B59-TetA.gb -r ../data/one-day-expt-reference-genome/A18-pUC.gb ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_3/*.fastq.gz"

## assemble the B30 ancestral strain RM6-200-6 using the K12-MG1655-NC_000913.gb reference genome.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/RM7-87-4 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B30-miniTn5-TetA.gb ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_4/*.fastq.gz"

## assemble the B30 + p15A ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/RM7-87-5 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B30-miniTn5-TetA.gb -r ../data/one-day-expt-reference-genome/A31-p15A.gb ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_5/*.fastq.gz"

## assemble the B30 + pUC ancestral strain RM6-176-18 using the K12-MG1655-NC_000913.gb reference genome.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/RM7-87-6 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B30-miniTn5-TetA.gb -r ../data/one-day-expt-reference-genome/A18-pUC.gb ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_6/*.fastq.gz"

################################################################################
## Apply mutations in the ancestral genomes to the K-12 MG1655 reference genome.
## NOTE: the output GFF3 files are combined references for the chromosome, the transposon, and the plasmid.

## apply mutations in the ancestral B59 genome, RM7-87-1, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/one-day-expt-genome-analysis/RM7-87-1.gff3 -f GFF3 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B59-TetA.gb ../results/one-day-expt-genome-analysis/RM7-87-1/output/output.gd

## apply mutations in the ancestral B59+A31 genome, RM7-87-2, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/one-day-expt-genome-analysis/RM7-87-2.gff3 -f GFF3 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B59-TetA.gb -r ../data/one-day-expt-reference-genome/A31-p15A.gb ../results/one-day-expt-genome-analysis/RM7-87-2/output/output.gd

## apply mutations in the ancestral B59+A18 genome, RM7-87-3, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/one-day-expt-genome-analysis/RM7-87-3.gff3 -f GFF3 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B59-TetA.gb -r ../data/one-day-expt-reference-genome/A18-pUC.gb ../results/one-day-expt-genome-analysis/RM7-87-3/output/output.gd

## apply mutations in the ancestral B30 genome, RM7-87-1, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/one-day-expt-genome-analysis/RM7-87-4.gff3 -f GFF3 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B30-miniTn5-TetA.gb ../results/one-day-expt-genome-analysis/RM7-87-4/output/output.gd

## apply mutations in the ancestral B30+A31 genome, RM7-87-5, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/one-day-expt-genome-analysis/RM7-87-5.gff3 -f GFF3 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B30-miniTn5-TetA.gb -r ../data/one-day-expt-reference-genome/A31-p15A.gb ../results/one-day-expt-genome-analysis/RM7-87-5/output/output.gd

## apply mutations in the ancestral B30+A18 genome, RM7-87-6, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/one-day-expt-genome-analysis/RM7-87-6.gff3 -f GFF3 -r ../data/one-day-expt-reference-genome/K12-MG1655-NC_000913.gb -r ../data/one-day-expt-reference-genome/B30-miniTn5-TetA.gb -r ../data/one-day-expt-reference-genome/A18-pUC.gb ../results/one-day-expt-genome-analysis/RM7-87-6/output/output.gd

################################################################################
## now, test the new references, by re-mapping reads.

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/remapped-RM7-87-1 -r ../results/one-day-expt-genome-analysis/RM7-87-1.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_1/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/remapped-RM7-87-2 -r ../results/one-day-expt-genome-analysis/RM7-87-2.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_2/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/remapped-RM7-87-3 -r ../results/one-day-expt-genome-analysis/RM7-87-3.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_3/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/remapped-RM7-87-4 -r ../results/one-day-expt-genome-analysis/RM7-87-4.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_4/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/remapped-RM7-87-5 -r ../results/one-day-expt-genome-analysis/RM7-87-5.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_5/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/one-day-expt-genome-analysis/remapped-RM7-87-6 -r ../results/one-day-expt-genome-analysis/RM7-87-6.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_6/*.fastq.gz"

################################################################################
## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 30, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.
## Only allow 5 max mismatches between the read and the reference.

############
## RM7.87.7 is K12 + B59 Tet0 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-7 -r ../results/one-day-expt-genome-analysis/RM7-87-1.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_7/*.fastq.gz"

## RM7.87.8 is K12 + B59 + A31 Tet0 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-8 -r ../results/one-day-expt-genome-analysis/RM7-87-2.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_8/*.fastq.gz"

## RM7.87.9 is K12 + B59 + A18 Tet0 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-9 -r ../results/one-day-expt-genome-analysis/RM7-87-3.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_9/*.fastq.gz"


######
## RM7.87.10 is K12 + B30 Tet0 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-10 -r ../results/one-day-expt-genome-analysis/RM7-87-4.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_10/*.fastq.gz"

## RM7.87.11 is K12 + B30 + A31 Tet0 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-11 -r ../results/one-day-expt-genome-analysis/RM7-87-5.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_11/*.fastq.gz"

## RM7.87.12 is K12 + B30 + A18 Tet0 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-12 -r ../results/one-day-expt-genome-analysis/RM7-87-6.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_12/*.fastq.gz"


############
## RM7.87.13 is K12 + B59 Tet0 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-13 -r ../results/one-day-expt-genome-analysis/RM7-87-1.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_13/*.fastq.gz"

## RM7.87.14 is K12 + B59 + A31 Tet0 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-14 -r ../results/one-day-expt-genome-analysis/RM7-87-2.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_14/*.fastq.gz"

## RM7.87.15 is K12 + B59 + A18 Tet0 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-15 -r ../results/one-day-expt-genome-analysis/RM7-87-3.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_15/*.fastq.gz"


######
## RM7.87.16 is K12 + B30 Tet0 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-16 -r ../results/one-day-expt-genome-analysis/RM7-87-4.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_16/*.fastq.gz"

## RM7.87.17 is K12 + B30 + A31 Tet0 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-17 -r ../results/one-day-expt-genome-analysis/RM7-87-5.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_17/*.fastq.gz"

## RM7.87.18 is K12 + B30 + A18 Tet0 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-18 -r ../results/one-day-expt-genome-analysis/RM7-87-6.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_18/*.fastq.gz"


############
## RM7.87.19 is K12 + B59 Tet0 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-19 -r ../results/one-day-expt-genome-analysis/RM7-87-1.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_19/*.fastq.gz"

## RM7.87.20 is K12 + B59 + A31 Tet0 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-20 -r ../results/one-day-expt-genome-analysis/RM7-87-2.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_20/*.fastq.gz"

## RM7.87.21 is K12 + B59 + A18 Tet0 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-21 -r ../results/one-day-expt-genome-analysis/RM7-87-3.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_21/*.fastq.gz"


######
## RM7.87.22 is K12 + B30 Tet0 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-22 -r ../results/one-day-expt-genome-analysis/RM7-87-4.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_22/*.fastq.gz"

## RM7.87.23 is K12 + B30 + A31 Tet0 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-23 -r ../results/one-day-expt-genome-analysis/RM7-87-5.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_23/*.fastq.gz"

## RM7.87.24 is K12 + B30 + A18 Tet0 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-24 -r ../results/one-day-expt-genome-analysis/RM7-87-6.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_24/*.fastq.gz"


############
## RM7.87.25 is K12 + B59 Tet5 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-25 -r ../results/one-day-expt-genome-analysis/RM7-87-1.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_25/*.fastq.gz"

## RM7.87.26 is K12 + B59 + A31 Tet5 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-26 -r ../results/one-day-expt-genome-analysis/RM7-87-2.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_26/*.fastq.gz"

## RM7.87.27 is K12 + B59 + A18 Tet5 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-27 -r ../results/one-day-expt-genome-analysis/RM7-87-3.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_27/*.fastq.gz"


######
## RM7.87.28 is K12 + B30 Tet5 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-28 -r ../results/one-day-expt-genome-analysis/RM7-87-4.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_28/*.fastq.gz"

## RM7.87.29 is K12 + B30 + A31 Tet5 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-29 -r ../results/one-day-expt-genome-analysis/RM7-87-5.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_29/*.fastq.gz"

## RM7.87.30 is K12 + B30 + A18 Tet5 evolved pop 1.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-30 -r ../results/one-day-expt-genome-analysis/RM7-87-6.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_30/*.fastq.gz"


############
## RM7.87.31 is K12 + B59 Tet5 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-31 -r ../results/one-day-expt-genome-analysis/RM7-87-1.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_31/*.fastq.gz"

## RM7.87.32 is K12 + B59 + A31 Tet5 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-32 -r ../results/one-day-expt-genome-analysis/RM7-87-2.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_32/*.fastq.gz"

## RM7.87.33 is K12 + B59 + A18 Tet5 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-33 -r ../results/one-day-expt-genome-analysis/RM7-87-3.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_33/*.fastq.gz"


######
## RM7.87.34 is K12 + B30 Tet5 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-34 -r ../results/one-day-expt-genome-analysis/RM7-87-4.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_34/*.fastq.gz"

## RM7.87.35 is K12 + B30 + A31 Tet5 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-35 -r ../results/one-day-expt-genome-analysis/RM7-87-5.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_35/*.fastq.gz"

## RM7.87.36 is K12 + B30 + A18 Tet5 evolved pop 2.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-36 -r ../results/one-day-expt-genome-analysis/RM7-87-6.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_36/*.fastq.gz"


############
## RM7.87.37 is K12 + B59 Tet5 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-37 -r ../results/one-day-expt-genome-analysis/RM7-87-1.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_37/*.fastq.gz"

## RM7.87.38 is K12 + B59 + A31 Tet5 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-38 -r ../results/one-day-expt-genome-analysis/RM7-87-2.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_38/*.fastq.gz"

## RM7.87.39 is K12 + B59 + A18 Tet5 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-39 -r ../results/one-day-expt-genome-analysis/RM7-87-3.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_39/*.fastq.gz"


######
## RM7.87.40 is K12 + B30 Tet5 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-40 -r ../results/one-day-expt-genome-analysis/RM7-87-4.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_40/*.fastq.gz"

## RM7.87.41 is K12 + B30 + A31 Tet5 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-41 -r ../results/one-day-expt-genome-analysis/RM7-87-5.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_41/*.fastq.gz"

## RM7.87.42 is K12 + B30 + A18 Tet5 evolved pop 3.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/one-day-expt-genome-analysis/RM7-87-42 -r ../results/one-day-expt-genome-analysis/RM7-87-6.gff3 ../data/one-day-expt-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_42/*.fastq.gz"
