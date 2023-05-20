#!/usr/bin/env bash

## assemble-DH5a-genomes.sh by Rohan Maddamsetti.
## This shell script does a bunch of genome assembly tasks using breseq and gdtools.
## COMMENT OUT all tasks that don't need to be done.
## Potential TODO: rewrite this task and re-script these tasks using SnakeMake.

## IMPORTANT TODO: there may be a bug in which the synthetic transposon is annotated as " Synthetic Tn5" instead of "Synthetic Tn5"
## in the GFF3 files made by gdtools. This causes breseq 0.37 to fail. I fixed this my deleting the space by hand in my GFF files, but
## this bug may reoccur when automatically generating the GFF3 files.

## IMPORTANT TODO: run the ancestral clone genomes in polymorphism mode, to see if mutations in "NEB5A_RS23175" show up.
## It is unclear whether the "evolved mutations" I seen in this pseudogene is in the ancestor, or if these are artifacts that show up
## in polymorphism mode. I should take a closer look here to clean up the downstream breseq runs. For now I just filter these mutations
## in process-gdiffs.py.

## IMPORTANT TODO: remake the reference GFF3 files, since I updated the Tn5 CDS annotation in my reference Genbank files.

################################################################################
## first, assemble the ancestral genomes using the NEB5-alpha-NZ_CP017100.gb reference genome.

## assemble the B30 + pUC ancestral strain RM6-176-18 using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM6-176-18 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B30-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A18-pUC.gb ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_18/*.fastq.gz"

## assemble the B30 ancestral strain RM6-200-6 using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM6-200-6 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B30-miniTn5-TetA.gb ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_6/*.fastq.gz"

## assemble the B30 + p15A ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM7-72-1 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B30-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A31-p15A.gb ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_1/*.fastq.gz"

## assemble the B20 + pUC ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM7-72-2 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B20-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A18-pUC.gb ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_2/*.fastq.gz"

## assemble the B20 ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM7-72-3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B20-miniTn5-TetA.gb ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_3/*.fastq.gz"

## assemble the B20 + p15A ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM7-72-4 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B20-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A31-p15A.gb ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_4/*.fastq.gz"

################################################################################
## apply mutations in the ancestral B30+A18 genome, RM6-176-18, to the NEB5-alpha-NZ_CP017100.gb reference genome.

## NOTE: the output GFF3 files are combined references for the chromosome, the transposon, and the plasmid.

## apply mutations in the ancestral B30+A18 (pUC) genome, RM6-176-18, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B30-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A18-pUC.gb ../results/DH5a-expt-genome-analysis/RM6-176-18/output/output.gd

## apply mutations in the ancestral B30+noplasmid genome, RM6-200-6, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B30-miniTn5-TetA.gb ../results/DH5a-expt-genome-analysis/RM6-200-6/output/output.gd

## apply mutations in the ancestral B30+A31 (p15A) genome, RM7-72-1, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B30-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A31-p15A.gb ../results/DH5a-expt-genome-analysis/RM7-72-1/output/output.gd

## apply mutations in the ancestral B20+pUC genome, RM7-72-2, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B20-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A18-pUC.gb ../results/DH5a-expt-genome-analysis/RM7-72-2/output/output.gd

## apply mutations in the ancestral B20+noplasmid genome, RM7-72-3, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B20-miniTn5-TetA.gb ../results/DH5a-expt-genome-analysis/RM7-72-3/output/output.gd

## apply mutations in the ancestral B20+p15A genome, RM7-72-4, to the NEB5-alpha-NZ_CP017100.gb reference genome.
##gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B20-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A31-p15A.gb ../results/DH5a-expt-genome-analysis/RM7-72-4/output/output.gd

################################################################################
## now, test the new references, by re-mapping reads.

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM6-176-18 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_18/*.fastq.gz"

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM6-200-6 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_6/*.fastq.gz"

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM7-72-1 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_1/*.fastq.gz"

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM7-72-2 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_2/*.fastq.gz"

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM7-72-3 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_3/*.fastq.gz"

##sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM7-72-4 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_4/*.fastq.gz"


################################################################################
## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 30, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.
## Only allow 5 max mismatches between the read and the reference.

## RM6.176.13-17 are DH5a + B30 + pUC Tet50 evolved pops 1-5.
#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-176-13 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_13/*.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-176-14 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_14/*.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-176-15 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_15/*.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-176-16 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_16/*.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-176-17 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_17/*.fastq.gz"


## RM6.200.1-5 are DH5a + B30 + Tet50 evolved pops 1-5.
#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-200-1 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_1/*.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-200-2 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_2/*.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-200-3 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_3/*.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-200-4 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_4/*.fastq.gz"

#sbatch -p scavenger --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM6-200-5 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_5/*.fastq.gz"

################################################################################
## now assemble the two clonal genomes from the RM6-200-3 mixed pop. sample.

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/RM6-147-1 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210630/063021_150/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/RM6-147-2 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210630/063021_151/*.fastq.gz"

################################################################################
## Let's assemble the DH5a + B30 + A18 pUC ancestor to see whether it has any secondary mutations (say the one mutation in all the evolved strains).
## D'oh! I already did this control! Let's compare this sample to that one, as an extra check.

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/RM7-60-5 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_60_5/*.fastq.gz"
################################################################################
## Now let's assemble the DH5a + B20 evolved populations.

## RM7.4.31-35 are DH5a + B20 evolved pops 1-5.
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-31 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_31/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-32 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_32/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-33 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_33/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-34 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_34/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-35 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_35/*.fastq.gz"

## RM7.4.36-40 are DH5a + B20 + A18 evolved pops 1-5.
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-36 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_36/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-37 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_37/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-38 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_38/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-39 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_39/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-40 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_40/*.fastq.gz"


## RM7.4.41-45 are DH5a + B30 + A31 evolved pops 1-5.
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-41 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_41/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-42 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_42/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-43 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_43/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-44 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_44/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-4-45 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_4_45/*.fastq.gz"


################################################################################
## RM7.72.5-9 is B20 + A31 (Tet 0 control populations).
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-5 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_5/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-6 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_6/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-7 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_7/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-8 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_8/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-9 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_9/*.fastq.gz"

## RM7.72.10-14 is B20 + A31 (Tet 50 populations).
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-10 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_10/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-11 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_11/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-12 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_12/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-13 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_13/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-14 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_14/*.fastq.gz"

## RM7.72.15-19 is B30+noplasmid (Tet 0 control populations).
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-15 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_15/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-16 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_16/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-17 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_17/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-18 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_18/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-19 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_19/*.fastq.gz"

## RM7.72.20-24 is B30+A18 (Tet 0 control populations).
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-20 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_20/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-21 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_21/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-22 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_22/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-23 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_23/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-24 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_24/*.fastq.gz"

## RM7.72.25-29 is B20 + A18 (Tet 0 control populations).
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-25 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_25/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-26 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_26/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-27 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_27/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-28 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_28/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-29 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_29/*.fastq.gz"

## RM7.72.30-34 is B20+noplasmid (Tet 0 control populations).
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-30 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_30/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-31 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_31/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-32 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_32/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-33 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_33/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-34 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_34/*.fastq.gz"

## RM7.72.35-39 is B30 + A31 (Tet 0 control populations).
#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-35 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_35/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-36 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_36/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-37 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_37/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-38 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_38/*.fastq.gz"

#sbatch -p scavenger -t 16:00:00 --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-72-39 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_39/*.fastq.gz"


################################################################################
## key control: do the ancestral clones with pUC show more false positives when assembled in polymorphism mode?
## idea instead of spike-in: assemble ancestral clones using polymorphism mode. if false positives are higher with pUC,
## then pUC ancestor clones should have more "polymorphism" than the other ancestors, despite being clonal.

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM7-60-5 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti220509/RM7_60_5/*.fastq.gz"

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM6-147-1 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210630/063021_150/*.fastq.gz"

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM6-147-2 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210630/063021_151/*.fastq.gz"

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM6-176-18 -r ../results/DH5a-expt-genome-analysis/RM6-176-18.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210805/RM6_176_18/*.fastq.gz"

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM6-200-6 -r ../results/DH5a-expt-genome-analysis/RM6-200-6.gff3 ../data/DH5a-genome-sequencing/MiGS_RohanMaddamsetti210917/RM6_200_6/*.fastq.gz"

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM7-72-1 -r ../results/DH5a-expt-genome-analysis/RM7-72-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_1/*.fastq.gz"

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM7-72-2 -r ../results/DH5a-expt-genome-analysis/RM7-72-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_2/*.fastq.gz"

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM7-72-3 -r ../results/DH5a-expt-genome-analysis/RM7-72-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_3/*.fastq.gz"

#sbatch --mem=8G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/clone-polymorphism-test/RM7-72-4 -r ../results/DH5a-expt-genome-analysis/RM7-72-4.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti220702/RM7_72_4/*.fastq.gz"

################################################################################
## new assemblies with B59. CRITICAL TODO: double-check the genome IDs with my notebook, before running.

## assemble the B59 ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM7-97-1 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B59-miniTn5-TetA.gb ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_1/*.fastq.gz"

## assemble the B59 + p15A ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM7-97-2 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B59-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A31-p15A.gb ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_2/*.fastq.gz"

## assemble the B20 + pUC ancestral strain using the NEB5-alpha-NZ_CP017100.gb reference genome.
#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/RM7-97-3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B59-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A18-pUC.gb ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_3/*.fastq.gz"

################################################################################

## apply mutations in the ancestral B59 genome to the NEB5-alpha-NZ_CP017100.gb reference genome.
#gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B59-miniTn5-TetA.gb ../results/DH5a-expt-genome-analysis/RM7-97-1/output/output.gd

## apply mutations in the ancestral B59+p15A genome to the NEB5-alpha-NZ_CP017100.gb reference genome.
#gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B59-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A31-p15A.gb ../results/DH5a-expt-genome-analysis/RM7-97-2/output/output.gd

## apply mutations in the ancestral B59+pUC genome to the NEB5-alpha-NZ_CP017100.gb reference genome.
#gdtools APPLY -o ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 -f GFF3 -r ../data/DH5a-genome-sequencing/NEB5-alpha-NZ_CP017100.gb -r ../data/DH5a-genome-sequencing/B59-miniTn5-TetA.gb -r ../data/DH5a-genome-sequencing/A18-pUC.gb ../results/DH5a-expt-genome-analysis/RM7-97-3/output/output.gd

################################################################################
## now, test the new references, by re-mapping reads.

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM7-97-1 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_1/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM7-97-2 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_2/*.fastq.gz"

#sbatch --mem=2G -c 1 --wrap="breseq -o ../results/DH5a-expt-genome-analysis/clones/remapped-RM7-97-3 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_3/*.fastq.gz"

################################################################################
## now assemble evolved mixed population samples with respect to their ancestor.
## set the minimum Phred score threshold for base quality to be 30, and polymorphic sites must have at least 4 reads on each strand supporting the mutation.
## Only allow 5 max mismatches between the read and the reference.

## RM7.97.4-8 are B59 + no plasmid Tet 0 evolved pops 1-5.
sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-4 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_4/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-5 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_5/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-6 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_6/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-7 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_7/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-8 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_8/*.fastq.gz"


## RM7.97.9-13 are B59 + no plasmid Tet 50 evolved pops 1-5.
sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-9 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_9/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-10 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_10/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-11 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_11/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-12 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_12/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-13 -r ../results/DH5a-expt-genome-analysis/RM7-97-1.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_13/*.fastq.gz"


## RM7.97.14-18 are B59 + p15A Tet 0 evolved pops 1-5.
sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-14 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_14/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-15 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_15/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-16 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_16/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-17 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_17/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-18 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_18/*.fastq.gz"


## RM7.97.19-23 are B59 + p15A Tet 50 evolved pops 1-5.
sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-19 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_19/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-20 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_20/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-21 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_21/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-22 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_22/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-23 -r ../results/DH5a-expt-genome-analysis/RM7-97-2.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_23/*.fastq.gz"


## RM7.97.24-28 are B59 + pUC Tet 0 evolved pops 1-5.
sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-24 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_24/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-25 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_25/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-26 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_26/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-27 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_27/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-28 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_28/*.fastq.gz"


## RM7.97.29-33 are B59 + pUC Tet 50 evolved pops 1-5.
sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-29 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_29/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-30 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_30/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-31 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_31/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-32 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_32/*.fastq.gz"

sbatch -p scavenger --mem=16G -c 8 --wrap="breseq -j 8 -p --polymorphism-minimum-variant-coverage-each-strand 4 -b 30 --maximum-read-mismatches 5 -o ../results/DH5a-expt-genome-analysis/mixed-pops/RM7-97-33 -r ../results/DH5a-expt-genome-analysis/RM7-97-3.gff3 ../data/DH5a-genome-sequencing/SeqCenter_RohanMaddamsetti221029/RM7_97_33/*.fastq.gz"
