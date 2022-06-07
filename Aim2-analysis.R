## Aim2-analysis.R by Rohan Maddamsetti.

## FOLLOW UP PAPER ANALYSIS TODO:
## look for evidence of recent diversification.
## In particular look more deeply at the result that AAA+ ATPases,
## and ATPases in general seem to be enriched in gene duplications.

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)
library(tidytext) ## for text mining with R.
library(forcats)


fancy_scientific <- function(x) {
    ## function for plotting better y-axis labels.
    ## see solution here for nice scientific notation on axes.
    ## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

## annotate source sequence as plasmid or chromosome.
genome.database <- read.csv("../results/chromosome-plasmid-table.csv")

## CRITICAL TODO: edit annotate-ecological-category.py to take care of the "blank" entries.
gbk.annotation <- as_tibble(read.csv("../results/computationally-annotated-gbk-annotation-table.csv")) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## get species name annotation from genome.database.
    left_join(genome.database) %>%
    ## CRITICAL STEP: remove the NCBI_Nucleotide_Accession and SequenceType columns.
    ## This is absolutely critical, otherwise each row is duplicated for every
    ## chromosome and plasmid, breaking the invariant that each row refers to one sequence,
    ## when we add this annotation to duplicate.proteins and singleton.proteins.
    select(-NCBI_Nucleotide_Accession, -SequenceType) %>%
    ## and we have to explicitly remove redundant rows now.
    distinct()

########################
## IMPORTANT NOTE: CHECK FOR CONSISTENCY IN protein_db_CDS_counts.csv!!!
## most importantly, examine replicons which are neither annotated as
## chromosomes or plasmids.
cds.counts <- read.csv("../results/protein_db_CDS_counts.csv")

protein.db.metadata <- genome.database %>%
    left_join(gbk.annotation) %>%
    left_join(cds.counts)    

## check out the different host and isolation source annotations.
chromosome.annotation <- protein.db.metadata %>%
    filter(SequenceType == "chromosome") %>%
    group_by(host, Annotation) %>%
    summarize(number = n()) %>%
    arrange(desc(number))

plasmid.annotation <- protein.db.metadata %>%
    filter(SequenceType == "plasmid") %>%
    group_by(host, Annotation) %>%
    summarize(number = n()) %>%
    arrange(desc(number))

## Check to see that the genomes and plasmids in the paper
## "Single-molecule sequencing to track plasmid diversity of
## hospital-associated carbapenemase-producing Enterobacteriaceae"
## by Conlan et al. (2014) in Science Advances are in AR.results.
## Strains with plasmids are in Table 2 of this paper.
Conlan.strains <- c("KPNIH1","KPNIH10","ECNIH3","ECNIH5","KPNIH27",
                    "CFNIH1","ECNIH2","KPNIH24","ECR091","KONIH1", 
                    "KPR0928","ECNIH4","PSNIH1","KPNIH32","PSNIH2",
                    "KPNIH33","ECONIH1","KPNIH30","KPNIH29","KPNIH31")

Conlan.strains %in% protein.db.metadata$Strain
################################################################################

## Simple analysis of recent protein duplications and HGT between
## chromosome and plasmids.

## remove all genes with the following keywords in the "product" annotation
IS.keywords <- "IS|transposon|Transposase|transposase|Transposable|transposable|hypothetical protein|Phage|phage|integrase|Integrase|tail|intron|Mobile|mobile|antitoxin|toxin|capsid|plasmid|Plasmid|conjug"

## remove all genes that match Elongation Factor Tu (2 copies in most bacteria).
EFTu.keywords <- "Tu | Tu|-Tu"

## now look at a few antibiotic-specific annotations.
antibiotic.keywords <- "lactamase|chloramphenicol|quinolone|antibiotic resistance|tetracycline|VanZ"

unknown.protein.keywords <- "unknown|Unknown|hypothetical|Hypothetical|Uncharacterized|Uncharacterised|uncharacterized|uncharacterised|DUF|unknow|putative protein in bacteria|Unassigned|unassigned"

## The regular expressions used by Zeevi et al. (2019).
## These are not used in this analysis, but nice to have on hand.
## Transposon: ‘transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon’
## plasmid: ‘relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation’
## phage: ‘capsid|phage|tail|head|tape measure|antiterminatio’
## other HGT mechanisms: ‘integrase|excision\S*|exo- nuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip’
## antibiotic resistance: ‘azole resistance|antibiotic resistance|TetR|tetracycline resistance|VanZ|betalactam\S*|beta-lactam|antimicrob\S*|lantibio\S*’.

## I want the sequence column for the duplicate genes,
## but not for the singletons, to save memory.

## import the 15GB file containing all proteins, including singletons.
## I can save a ton of memory if I don't import the sequence column,
## and by using the data.table package for import.
all.proteins <- data.table::fread("../results/all-proteins.csv",
                                  drop="sequence") %>%
    left_join(gbk.annotation)
## I am doing a left_join here, because I eventually want to
## predict where the unannotated strains come from.

all.singleton.proteins <- all.proteins %>% filter(count == 1)

## free up memory by deallocating all.proteins,
rm(all.proteins)
## and running garbage collection.
gc()

## read in duplicate proteins with sequences, using a separate file.
all.duplicate.proteins <- read.csv("../results/duplicate-proteins.csv") %>% left_join(gbk.annotation)
    ## now merge with gbk annotation.
    ## I am doing a left_join here, because I want the NA Manual_Accessions
    ## in order to predict where these unannotated strains come from.

## Some strains in chromosome-and-plasmid-table.csv and
## gbk-annotation-table.csv are missing from
## all-proteins.csv and duplicate-proteins.csv.
## These should be the genomes that do not have
## CDS annotated in their GFF annotation.
## list the 926 strains missing from the singletons data.
missing.ones <- gbk.annotation %>%
    filter(!(Annotation_Accession %in% all.singleton.proteins$Annotation_Accession))

#################
## For the first part of the data analysis, ignore genes from Unannotated isolates.
duplicate.proteins <- all.duplicate.proteins %>%
    filter(Annotation != "Unannotated") %>%
    ## TODO: in the future, there should not be any "blank" genomes.
    filter(Annotation != "blank")

singleton.proteins <- all.singleton.proteins %>%
    filter(Annotation != "Unannotated") %>%
    ## TODO: in the future, there should not be any "blank" genomes
    filter(Annotation != "blank")

## for now, remove these data structures from memory, since they are not used
## in any analyses yet.
rm(all.singleton.proteins)

## call garbage collector to free up memory.
gc()
