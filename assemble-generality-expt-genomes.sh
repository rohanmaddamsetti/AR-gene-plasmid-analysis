#!/usr/bin/env bash

## assemble-generality-expt-genomes.sh by Rohan Maddamsetti.
## This shell script does a bunch of genome assembly tasks using breseq and gdtools.
## COMMENT OUT all tasks that don't need to be done.

## IMPORTANT TODO: there may be a bug in which the synthetic transposon is annotated as " Synthetic Tn5" instead of "Synthetic Tn5"
## in the GFF3 files made by gdtools. This causes breseq 0.37 to fail. I fixed this my deleting the space by hand in my GFF files, but
## this bug may reoccur when automatically generating the GFF3 files.

################################################################################
## first, assemble the ancestral genomes using the K12-MG1655-NC_000913.gb reference genome.

## here, assemble ancestors for the generality experiment varying the transposase.

## assemble the K-12 + B107 ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-1 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B107-IS85-TetA.gb ../data/generality-expts-genome-data/8336-S1_S1_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S1_S1_L001_R2_001.fastq.gz"

## assemble the K-12 + B111 ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-2 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B111-IS88-TetA.gb ../data/generality-expts-genome-data/8336-S2_S2_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S2_S2_L001_R2_001.fastq.gz"

## assemble the K-12 + B123 ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B123-IS101-TetA.gb ../data/generality-expts-genome-data/8336-S3_S3_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S3_S3_L001_R2_001.fastq.gz"

## assemble the K-12 + B134 ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-4 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B134-IS112-TetA.gb ../data/generality-expts-genome-data/8336-S4_S4_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S4_S4_L001_R2_001.fastq.gz"

## assemble the K-12 + B142 ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-5 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B142-IS120-TetA.gb ../data/generality-expts-genome-data/8336-S5_S5_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S5_S5_L001_R2_001.fastq.gz"

## Not a typo-- we omitted RM7-95-6 from this sequencing run.
## assemble the K-12 + B107 + p15A plasmid ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-7 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B107-IS85-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/8336-S6_S6_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S6_S6_L001_R2_001.fastq.gz"

## assemble the K-12 + B111 + p15A plasmid ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-8 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B111-IS88-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/8336-S7_S7_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S7_S7_L001_R2_001.fastq.gz"

## assemble the K-12 + B123 + p15A plasmid ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-9 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B123-IS101-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/8336-S8_S8_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S8_S8_L001_R2_001.fastq.gz"

## now, assemble ancestors for the generality experiment varying the antibiotic selection.
## assemble the K-12 + B109 ancestor.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-1 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B109-IS87-TetA.gb ../data/generality-expts-genome-data/8336-S9_S9_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S9_S9_L001_R2_001.fastq.gz"

## assemble the K-12 + B110 ancestor.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-2 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B110-IS88-TetA.gb ../data/generality-expts-genome-data/8336-S10_S10_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S10_S10_L001_R2_001.fastq.gz"

## assemble the K-12 + B109 + p15A plasmid ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B109-IS87-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/8336-S11_S11_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S11_S11_L001_R2_001.fastq.gz"

## assemble the K-12 + B110 + p15A plasmid ancestral strain.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-4 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B110-IS88-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/8336-S12_S12_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S12_S12_L001_R2_001.fastq.gz"

## assemble the K-12 + B90 ancestor.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-5 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B90-miniTn5-SmR.gb ../data/generality-expts-genome-data/8336-S13_S13_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S13_S13_L001_R2_001.fastq.gz"

## assemble the K-12 + B91 ancestor.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-6 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B91-miniTn5-KanR.gb ../data/generality-expts-genome-data/8336-S14_S14_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S14_S14_L001_R2_001.fastq.gz"

## assemble the K-12 + B92 ancestor.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-7 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B92-miniTn5-AmpR.gb ../data/generality-expts-genome-data/8336-S15_S15_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S15_S15_L001_R2_001.fastq.gz"

## assemble the K-12 + B95 ancestor.
sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-8 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B95-miniTn5-CmR.gb ../data/generality-expts-genome-data/8336-S16_S16_L001_R1_001.fastq.gz ../data/generality-expts-genome-data/8336-S16_S16_L001_R2_001.fastq.gz"


################################################################################
## Apply mutations in the ancestral genomes to the K-12 MG1655 reference genome.
## NOTE: the output GFF3 files are combined references for the chromosome, the transposon, and the plasmid.

## apply mutations in the ancestral B59 genome, RM7-87-1, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-87-1.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B59-TetA.gb ../results/generality-expts-genome-analysis/RM7-87-1/output/output.gd

## apply mutations in the ancestral B59+A31 genome, RM7-87-2, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-87-2.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B59-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../results/generality-expts-genome-analysis/RM7-87-2/output/output.gd

## apply mutations in the ancestral B59+A18 genome, RM7-87-3, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-87-3.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B59-TetA.gb -r ../data/generality-expts-reference-genome/A18-pUC.gb ../results/generality-expts-genome-analysis/RM7-87-3/output/output.gd

## apply mutations in the ancestral B30 genome, RM7-87-1, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-87-4.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B30-miniTn5-TetA.gb ../results/generality-expts-genome-analysis/RM7-87-4/output/output.gd

## apply mutations in the ancestral B30+A31 genome, RM7-87-5, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-87-5.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B30-miniTn5-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../results/generality-expts-genome-analysis/RM7-87-5/output/output.gd

## apply mutations in the ancestral B30+A18 genome, RM7-87-6, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-87-6.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B30-miniTn5-TetA.gb -r ../data/generality-expts-reference-genome/A18-pUC.gb ../results/generality-expts-genome-analysis/RM7-87-6/output/output.gd

################################################################################
## now, test the new references, by re-mapping reads.

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-87-1 -r ../results/generality-expts-genome-analysis/RM7-87-1.gff3 ../data/generality-expts-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_1/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-87-2 -r ../results/generality-expts-genome-analysis/RM7-87-2.gff3 ../data/generality-expts-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_2/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-87-3 -r ../results/generality-expts-genome-analysis/RM7-87-3.gff3 ../data/generality-expts-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_3/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-87-4 -r ../results/generality-expts-genome-analysis/RM7-87-4.gff3 ../data/generality-expts-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_4/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-87-5 -r ../results/generality-expts-genome-analysis/RM7-87-5.gff3 ../data/generality-expts-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_5/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-87-6 -r ../results/generality-expts-genome-analysis/RM7-87-6.gff3 ../data/generality-expts-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_6/*.fastq.gz"

################################################################################
## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 30, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.
## Only allow 5 max mismatches between the read and the reference.

############
## EXAMPLE TO USE FOR NEW DATA
## RM7.87.7 is K12 + B59 Tet0 evolved pop 1.
##sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-87-7 -r ../results/generality-expts-genome-analysis/RM7-87-1.gff3 ../data/generality-expts-genome-data/SeqCenter_RohanMaddamsetti220907/RM7_87_7/*.fastq.gz"
