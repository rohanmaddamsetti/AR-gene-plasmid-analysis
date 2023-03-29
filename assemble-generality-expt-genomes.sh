#!/usr/bin/env bash

## assemble-generality-expt-genomes.sh by Rohan Maddamsetti.
## This shell script does a bunch of genome assembly tasks using breseq and gdtools.
## COMMENT OUT all tasks that don't need to be done.

################################################################################
## first, assemble the ancestral genomes using the K12-MG1655-NC_000913.gb reference genome.
## assemble ancestors for the generality experiment varying the transposase.

## assemble the K-12 + B107 ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-1 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B107-IS85-TetA.gb ../data/generality-expts-genome-data/RM7-95-1/*.fastq.gz"

## assemble the K-12 + B111 ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-2 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B111-IS89-TetA.gb ../data/generality-expts-genome-data/RM7-95-2/*.fastq.gz"

## assemble the K-12 + B123 ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B123-IS101-TetA.gb ../data/generality-expts-genome-data/RM7-95-3/*.fastq.gz"

## assemble the K-12 + B134 ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-4 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B134-IS112-TetA.gb ../data/generality-expts-genome-data/RM7-95-4/*.fastq.gz"

## assemble the K-12 + B142 ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-5 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B142-IS120-TetA.gb ../data/generality-expts-genome-data/RM7-95-5/*.fastq.gz"

## Not a typo-- we omitted RM7-95-6 from this sequencing run.
## assemble the K-12 + B107 + p15A plasmid ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-7 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B107-IS85-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/RM7-95-7/*.fastq.gz"

## assemble the K-12 + B111 + p15A plasmid ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-8 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B111-IS89-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/RM7-95-8/*.fastq.gz"

## assemble the K-12 + B123 + p15A plasmid ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-95-9 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B123-IS101-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/RM7-95-9/*.fastq.gz"

## now, assemble ancestors for the generality experiment varying the antibiotic selection.
## assemble the K-12 + B109 ancestor.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-1 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B109-IS87-TetA.gb ../data/generality-expts-genome-data/RM7-114-1/*.fastq.gz"

## assemble the K-12 + B110 ancestor.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-2 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B110-IS88-TetA.gb ../data/generality-expts-genome-data/RM7-114-2/*.fastq.gz"

## assemble the K-12 + B109 + p15A plasmid ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B109-IS87-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/RM7-114-3/*.fastq.gz"

## assemble the K-12 + B110 + p15A plasmid ancestral strain.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-4 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B110-IS88-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../data/generality-expts-genome-data/RM7-114-4/*.fastq.gz"

## assemble the K-12 + B90 ancestor.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-5 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B90-miniTn5-SmR.gb ../data/generality-expts-genome-data/RM7-114-5/*.fastq.gz"

## assemble the K-12 + B91 ancestor.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-6 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B91-miniTn5-KanR.gb ../data/generality-expts-genome-data/RM7-114-6/*.fastq.gz"

## assemble the K-12 + B92 ancestor.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-7 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B92-miniTn5-AmpR.gb ../data/generality-expts-genome-data/RM7-114-7/*.fastq.gz"

## assemble the K-12 + B95 ancestor.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/RM7-114-8 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B95-miniTn5-CmR.gb ../data/generality-expts-genome-data/RM7-114-8/*.fastq.gz"


################################################################################
## Apply mutations in the ancestral genomes to the K-12 MG1655 reference genome.
## NOTE: the output GFF3 files are combined references for the chromosome, the transposon, and the plasmid.
## NOTE 2: RM7-95-2, RM7-95-4, and RM7-114-4 show the SAME intergenic (‑512/‑73) flhD ← / → insH1 mutation after remapping reads.
## to deal with this, I apply and remap mutations twice for these samples. The first pass gd files for these samples are called "prep-RM7-95-2", etc.

## apply mutations in the ancestral B107 genome, RM7-95-1, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-95-1.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B107-IS85-TetA.gb ../results/generality-expts-genome-analysis/RM7-95-1/output/output.gd

## apply mutations in the ancestral B111 genome, RM7-95-2, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/prep-RM7-95-2.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B111-IS89-TetA.gb ../results/generality-expts-genome-analysis/RM7-95-2/output/output.gd

## apply mutations in the ancestral B123 genome, RM7-95-3, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-95-3.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B123-IS101-TetA.gb ../results/generality-expts-genome-analysis/RM7-95-3/output/output.gd

## apply mutations in the ancestral B134 genome, RM7-95-4, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/prep-RM7-95-4.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B134-IS112-TetA.gb ../results/generality-expts-genome-analysis/RM7-95-4/output/output.gd

## apply mutations in the ancestral B142 genome, RM7-95-5, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-95-5.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B142-IS120-TetA.gb ../results/generality-expts-genome-analysis/RM7-95-5/output/output.gd

## apply mutations in the ancestral B107 + p15A genome, RM7-95-7, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-95-7.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B107-IS85-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../results/generality-expts-genome-analysis/RM7-95-7/output/output.gd

## apply mutations in the ancestral B111 + p15A genome, RM7-95-8, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-95-8.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B111-IS89-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../results/generality-expts-genome-analysis/RM7-95-8/output/output.gd

## apply mutations in the ancestral B123 + p15A genome, RM7-95-9, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-95-9.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B123-IS101-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../results/generality-expts-genome-analysis/RM7-95-9/output/output.gd

## now, apply mutations for the ancestors in the generality experiment varying the antibiotic selection.

## apply mutations in the ancestral B109 genome, RM7-114-1, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-114-1.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B109-IS87-TetA.gb ../results/generality-expts-genome-analysis/RM7-114-1/output/output.gd

## apply mutations in the ancestral B110 genome, RM7-114-2, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-114-2.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B110-IS88-TetA.gb ../results/generality-expts-genome-analysis/RM7-114-2/output/output.gd

## apply mutations in the ancestral B109 + p15A genome, RM7-114-3, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-114-3.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B109-IS87-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../results/generality-expts-genome-analysis/RM7-114-3/output/output.gd

## apply mutations in the ancestral B110 + p15A genome, RM7-114-4, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/prep-RM7-114-4.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B110-IS88-TetA.gb -r ../data/generality-expts-reference-genome/A31-p15A.gb ../results/generality-expts-genome-analysis/RM7-114-4/output/output.gd

## apply mutations in the ancestral B90 genome, RM7-114-5, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-114-5.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B90-miniTn5-SmR.gb ../results/generality-expts-genome-analysis/RM7-114-5/output/output.gd

## apply mutations in the ancestral B91 genome, RM7-114-6, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-114-6.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B91-miniTn5-KanR.gb ../results/generality-expts-genome-analysis/RM7-114-6/output/output.gd

## apply mutations in the ancestral B92 genome, RM7-114-7, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-114-7.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B92-miniTn5-AmpR.gb ../results/generality-expts-genome-analysis/RM7-114-7/output/output.gd

## apply mutations in the ancestral B95 genome, RM7-114-8, to the K12-MG1655-NC_000913.gb reference genome.
#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-114-8.gff3 -f GFF3 -r ../data/generality-expts-reference-genome/K12-MG1655-NC_000913.gb -r ../data/generality-expts-reference-genome/B95-miniTn5-CmR.gb ../results/generality-expts-genome-analysis/RM7-114-8/output/output.gd

################################################################################
## now, test the new references, by re-mapping reads.
## NOTE: prefixes for RM7-95-2, RM7-95-4, RM7-114-4 are "first-remapping", since we will APPLY and remap mutations a second time for these samples.
## The GD files for these samples are "prep-RM7-95-2" etc.

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-95-1 -r ../results/generality-expts-genome-analysis/RM7-95-1.gff3 ../data/generality-expts-genome-data/RM7-95-1/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/first-remapping-RM7-95-2 -r ../results/generality-expts-genome-analysis/prep-RM7-95-2.gff3 ../data/generality-expts-genome-data/RM7-95-2/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-95-3 -r ../results/generality-expts-genome-analysis/RM7-95-3.gff3 ../data/generality-expts-genome-data/RM7-95-3/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/first-remapping-RM7-95-4 -r ../results/generality-expts-genome-analysis/prep-RM7-95-4.gff3 ../data/generality-expts-genome-data/RM7-95-4/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-95-5 -r ../results/generality-expts-genome-analysis/RM7-95-5.gff3 ../data/generality-expts-genome-data/RM7-95-5/*.fastq.gz"

## again, not a typo. We omitted RM7-95-6 since not needed.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-95-7 -r ../results/generality-expts-genome-analysis/RM7-95-7.gff3 ../data/generality-expts-genome-data/RM7-95-7/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-95-8 -r ../results/generality-expts-genome-analysis/RM7-95-8.gff3 ../data/generality-expts-genome-data/RM7-95-8/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-95-9 -r ../results/generality-expts-genome-analysis/RM7-95-9.gff3 ../data/generality-expts-genome-data/RM7-95-9/*.fastq.gz"

## now remap reads to the generality over antibiotics genomes.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-114-1 -r ../results/generality-expts-genome-analysis/RM7-114-1.gff3 ../data/generality-expts-genome-data/RM7-114-1/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-114-2 -r ../results/generality-expts-genome-analysis/RM7-114-2.gff3 ../data/generality-expts-genome-data/RM7-114-2/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-114-3 -r ../results/generality-expts-genome-analysis/RM7-114-3.gff3 ../data/generality-expts-genome-data/RM7-114-3/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/first-remapping-RM7-114-4 -r ../results/generality-expts-genome-analysis/prep-RM7-114-4.gff3 ../data/generality-expts-genome-data/RM7-114-4/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-114-5 -r ../results/generality-expts-genome-analysis/RM7-114-5.gff3 ../data/generality-expts-genome-data/RM7-114-5/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-114-6 -r ../results/generality-expts-genome-analysis/RM7-114-6.gff3 ../data/generality-expts-genome-data/RM7-114-6/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-114-7 -r ../results/generality-expts-genome-analysis/RM7-114-7.gff3 ../data/generality-expts-genome-data/RM7-114-7/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-114-8 -r ../results/generality-expts-genome-analysis/RM7-114-8.gff3 ../data/generality-expts-genome-data/RM7-114-8/*.fastq.gz"

################################################################################
## reapply mutations to RM7-95-2, RM7-95-4, RM7-114-4.

#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-95-2.gff3 -f GFF3 -r ../results/generality-expts-genome-analysis/prep-RM7-95-2.gff3 ../results/generality-expts-genome-analysis/first-remapping-RM7-95-2/output/output.gd

#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-95-4.gff3 -f GFF3 -r ../results/generality-expts-genome-analysis/prep-RM7-95-4.gff3 ../results/generality-expts-genome-analysis/first-remapping-RM7-95-4/output/output.gd

#gdtools APPLY -o ../results/generality-expts-genome-analysis/RM7-114-4.gff3 -f GFF3 -r ../results/generality-expts-genome-analysis/prep-RM7-114-4.gff3 ../results/generality-expts-genome-analysis/first-remapping-RM7-114-4/output/output.gd

################################################################################
## remap reads one more time for RM7-95-2, RM7-95-4, RM7-114-4.

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-95-2 -r ../results/generality-expts-genome-analysis/RM7-95-2.gff3 ../data/generality-expts-genome-data/RM7-95-2/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-95-4 -r ../results/generality-expts-genome-analysis/RM7-95-4.gff3 ../data/generality-expts-genome-data/RM7-95-4/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/generality-expts-genome-analysis/remapped-RM7-114-4 -r ../results/generality-expts-genome-analysis/RM7-114-4.gff3 ../data/generality-expts-genome-data/RM7-114-4/*.fastq.gz"

################################################################################
## Transposon generality experiment.
## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 30, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.
## Only allow 5 max mismatches between the read and the reference.

## RM7-114-9-13:  Tet0 controls, no plasmid.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-9 -r ../results/generality-expts-genome-analysis/RM7-95-1.gff3 ../data/generality-expts-genome-data/RM7-114-9/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-10 -r ../results/generality-expts-genome-analysis/RM7-95-2.gff3 ../data/generality-expts-genome-data/RM7-114-10/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-11 -r ../results/generality-expts-genome-analysis/RM7-95-3.gff3 ../data/generality-expts-genome-data/RM7-114-11/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-12 -r ../results/generality-expts-genome-analysis/RM7-114-1.gff3 ../data/generality-expts-genome-data/RM7-114-12/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-13 -r ../results/generality-expts-genome-analysis/RM7-114-2.gff3 ../data/generality-expts-genome-data/RM7-114-13/*.fastq.gz"

## RM7-114-14-18:  Tet0 controls, p15A plasmid.
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-14 -r ../results/generality-expts-genome-analysis/RM7-95-7.gff3 ../data/generality-expts-genome-data/RM7-114-14/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-15 -r ../results/generality-expts-genome-analysis/RM7-95-8.gff3 ../data/generality-expts-genome-data/RM7-114-15/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-16 -r ../results/generality-expts-genome-analysis/RM7-95-9.gff3 ../data/generality-expts-genome-data/RM7-114-16/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-17 -r ../results/generality-expts-genome-analysis/RM7-114-3.gff3 ../data/generality-expts-genome-data/RM7-114-17/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-114-18 -r ../results/generality-expts-genome-analysis/RM7-114-4.gff3 ../data/generality-expts-genome-data/RM7-114-18/*.fastq.gz"

## RM7-115-1-10:  Tet5 replicates 1

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-1 -r ../results/generality-expts-genome-analysis/RM7-95-1.gff3 ../data/generality-expts-genome-data/RM7-115-1/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-2 -r ../results/generality-expts-genome-analysis/RM7-95-2.gff3 ../data/generality-expts-genome-data/RM7-115-2/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-3 -r ../results/generality-expts-genome-analysis/RM7-95-3.gff3 ../data/generality-expts-genome-data/RM7-115-3/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-4 -r ../results/generality-expts-genome-analysis/RM7-114-1.gff3 ../data/generality-expts-genome-data/RM7-115-4/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-5 -r ../results/generality-expts-genome-analysis/RM7-114-2.gff3 ../data/generality-expts-genome-data/RM7-115-5/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-6 -r ../results/generality-expts-genome-analysis/RM7-95-7.gff3 ../data/generality-expts-genome-data/RM7-115-6/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-7 -r ../results/generality-expts-genome-analysis/RM7-95-8.gff3 ../data/generality-expts-genome-data/RM7-115-7/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-8 -r ../results/generality-expts-genome-analysis/RM7-95-9.gff3 ../data/generality-expts-genome-data/RM7-115-8/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-9 -r ../results/generality-expts-genome-analysis/RM7-114-3.gff3 ../data/generality-expts-genome-data/RM7-115-9/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-10 -r ../results/generality-expts-genome-analysis/RM7-114-4.gff3 ../data/generality-expts-genome-data/RM7-115-10/*.fastq.gz"

## RM7-115-11-20: Tet5 replicates 2

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-11 -r ../results/generality-expts-genome-analysis/RM7-95-1.gff3 ../data/generality-expts-genome-data/RM7-115-11/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-12 -r ../results/generality-expts-genome-analysis/RM7-95-2.gff3 ../data/generality-expts-genome-data/RM7-115-12/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-13 -r ../results/generality-expts-genome-analysis/RM7-95-3.gff3 ../data/generality-expts-genome-data/RM7-115-13/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-14 -r ../results/generality-expts-genome-analysis/RM7-114-1.gff3 ../data/generality-expts-genome-data/RM7-115-14/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-15 -r ../results/generality-expts-genome-analysis/RM7-114-2.gff3 ../data/generality-expts-genome-data/RM7-115-15/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-16 -r ../results/generality-expts-genome-analysis/RM7-95-7.gff3 ../data/generality-expts-genome-data/RM7-115-16/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-17 -r ../results/generality-expts-genome-analysis/RM7-95-8.gff3 ../data/generality-expts-genome-data/RM7-115-17/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-18 -r ../results/generality-expts-genome-analysis/RM7-95-9.gff3 ../data/generality-expts-genome-data/RM7-115-18/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-19 -r ../results/generality-expts-genome-analysis/RM7-114-3.gff3 ../data/generality-expts-genome-data/RM7-115-19/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-20 -r ../results/generality-expts-genome-analysis/RM7-114-4.gff3 ../data/generality-expts-genome-data/RM7-115-20/*.fastq.gz"

################################################################################
## Antibiotic generality experiment.
## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 30, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.
## Only allow 5 max mismatches between the read and the reference.

## RM7-115-21-24:  No antibiotic controls (B90, B91, B92, B95)
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-21 -r ../results/generality-expts-genome-analysis/RM7-114-5.gff3 ../data/generality-expts-genome-data/RM7-115-21/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-22 -r ../results/generality-expts-genome-analysis/RM7-114-6.gff3 ../data/generality-expts-genome-data/RM7-115-22/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-23 -r ../results/generality-expts-genome-analysis/RM7-114-7.gff3 ../data/generality-expts-genome-data/RM7-115-23/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-24 -r ../results/generality-expts-genome-analysis/RM7-114-8.gff3 ../data/generality-expts-genome-data/RM7-115-24/*.fastq.gz"

## RM7-115-25-28: Spec 250, Kan 250, Carb 2000, Cm70 replicates 1 (B90, B91, B92, B95)
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-25 -r ../results/generality-expts-genome-analysis/RM7-114-5.gff3 ../data/generality-expts-genome-data/RM7-115-25/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-26 -r ../results/generality-expts-genome-analysis/RM7-114-6.gff3 ../data/generality-expts-genome-data/RM7-115-26/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-27 -r ../results/generality-expts-genome-analysis/RM7-114-7.gff3 ../data/generality-expts-genome-data/RM7-115-27/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-28 -r ../results/generality-expts-genome-analysis/RM7-114-8.gff3 ../data/generality-expts-genome-data/RM7-115-28/*.fastq.gz"

## RM7-115-29-32: Spec 250, Kan 250, Carb 2000, Cm70 replicates 2 (B90, B91, B92, B95)
sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-29 -r ../results/generality-expts-genome-analysis/RM7-114-5.gff3 ../data/generality-expts-genome-data/RM7-115-29/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-30 -r ../results/generality-expts-genome-analysis/RM7-114-6.gff3 ../data/generality-expts-genome-data/RM7-115-30/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-31 -r ../results/generality-expts-genome-analysis/RM7-114-7.gff3 ../data/generality-expts-genome-data/RM7-115-31/*.fastq.gz"

sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/generality-expts-genome-analysis/RM7-115-32 -r ../results/generality-expts-genome-analysis/RM7-114-8.gff3 ../data/generality-expts-genome-data/RM7-115-32/*.fastq.gz"

