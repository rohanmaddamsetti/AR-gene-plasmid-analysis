## Aim1-analysis.R by Rohan Maddamsetti.

## analyse the distribution of AR genes on chromosomes versus plasmids in
## fully-sequenced genomes and plasmids in the NCBI Nucleotide database.

## CRITICAL TODO: re-annotate manually-curated-gbk-annotation-table.csv
## using the updated gbk-annotation table. This will add ~300 extra genomes
## to the analysis.

## CRITICAL ANALYSIS TODO: look for evidence of recent diversification.
## In particular look more deeply at the result that AAA+ ATPases,
## and ATPases in general seem to be enriched in gene duplications.

## CRITICAL ANALYSIS TODO: Make sure numbers in genome.database,
## gbk.annotation, and all.proteins, duplicate.proteins, and
## singleton.proteins, in terms of number of isolates in each
## category, are COMPLETELY consistent with each other.

## Potential TODO if reviewers or colleagues think necessary:
## annotate resistance genes using the CARD RGI tool
## https://github.com/arpcard/rgi
## and look at singleton and duplicate RGs identified by this workflow.

library(tidyverse)
library(cowplot)
library(data.table)

## annotate source sequence as plasmid or chromosome.
genome.database <- read.csv("../results/AR-gene-duplication/chromosome-plasmid-table.csv")

raw.gbk.annotation <- as_tibble(read.csv("../results/AR-gene-duplication/gbk-annotation-table.csv"))

length(unique(genome.database$Annotation_Accession))

## This uses the manually-annotated data that I have so far.
## IMPORTANT TODO: update manual annotation version,
## with a script that generates this file.
gbk.annotation <- as_tibble(read.csv("../data/manually-curated-gbk-annotation-table.csv")) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Manual_Annotation = replace_na(Manual_Annotation,"Unannotated"))

########################
## TEMPORARY HACK FOR SELF-CONSISTENCY:
gbk.annotation <- gbk.annotation %>% filter(Annotation_Accession %in% genome.database$Annotation_Accession)
########################
## IMPORTANT NOTE: CHECK FOR CONSISTENCY IN protein_db_CDS_counts.csv!!!
## most importantly, examine replicons which are neither annotated as
## chromosomes or plasmids.
cds.counts <- read.csv("../results/AR-gene-duplication/protein_db_CDS_counts.csv")

protein.db.metadata <- genome.database %>%
    left_join(gbk.annotation) %>%
    left_join(cds.counts)    

## check out the different host and isolation source annotations.
chromosome.annotation <- protein.db.metadata %>%
    filter(SequenceType == "chromosome") %>%
    group_by(host, Manual_Annotation) %>%
    summarize(number = n()) %>%
    arrange(desc(number))

plasmid.annotation <- protein.db.metadata %>%
    filter(SequenceType == "plasmid") %>%
    group_by(host, Manual_Annotation) %>%
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
## PSNIH2 is missing, but the others are in genome.database
## and AR.results. PSNIH2 is not in the original
## prokaryotes.txt genome reports file.
################################################################################

## Naive analysis of recent protein duplications and HGT between
## chromosome and plasmids.

## remove all genes with the following keywords in the "product" annotation
## TODO: plot distribution of JUST genes with these annotations.
IS.keywords <- "IS|transposon|Transposase|transposase|hypothetical protein|Phage|phage|integrase|Integrase|tail|intron|Mobile|mobile|antitoxin|toxin"

## now look at a few antibiotic-specific annotations.
antibiotic.keywords <- "lactamase|chloramphenicol|quinolone|antibiotic resistance|tetracycline|VanZ"

## IMPORTANT TODO: Use the same regular expressions used by Zeevi et al. (2019).
## Transposon: ‘transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon’
## plasmid: ‘relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation’
## phage: ‘capsid|phage|tail|head|tape measure|antiterminatio’
## other HGT mechanisms: ‘integrase|excision\S*|exo- nuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip’
## antibiotic resistance: ‘azole resistance|antibiotic resistance|TetR|tetracycline resistance|VanZ|betalactam\S*|beta-lactam|antimicrob\S*|lantibio\S*’.


## I want the sequence column for the duplicate genes,
## but not for the singletons, to save memory.

## import the 12GB file containing all proteins, including singletons.
## I can save a ton of memory if I don't import the sequence column,
## and by using the data.table package for import.
all.proteins <- data.table::fread("../results/AR-gene-duplication/all-proteins.csv",
                                  drop="sequence") %>%
    ## now merge with gbk annotation.
    ## I am doing a left_join here, because I want the NA Manual_Accessions
    ## in order to predict where these unannotated strains come from.
    left_join(gbk.annotation) ##%>%
    ## refer to NA annotations as "Unannotated".
    ##mutate(Manual_Annotation = replace_na(Manual_Annotation,"Unannotated"))

singleton.proteins <- all.proteins %>% filter(count == 1)

## free up memory by deallocating all.proteins,
rm(all.proteins)
## and running garbage collection.
gc()

## read in duplicate proteins with sequences, using a separate file.
duplicate.proteins <- read.csv("../results/AR-gene-duplication/duplicate-proteins.csv") %>%
    ## now merge with gbk annotation.
    ## I am doing a left_join here, because I want the NA Manual_Accessions
    ## in order to predict where these unannotated strains come from.
    left_join(gbk.annotation) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Manual_Annotation = replace_na(Manual_Annotation,"Unannotated"))

## Some strains in chromosome-and-plasmid-table.csv and
## gbk-annotation-table.csv are missing from
## all-proteins.csv and duplicate-proteins.csv.
## These should be the genomes that do not have
## CDS annotated in their GFF annotation.
## list the 634 strains missing from the singletons data.
missing.ones <- raw.gbk.annotation %>%
    filter(!(Annotation_Accession %in% singleton.proteins$Annotation_Accession))

#################
## TEMPORARY HACK FOR SELF-CONSISTENCY:
## I shouldn't have to do this, after updating the manual annotation.
duplicate.proteins <- duplicate.proteins %>%
    filter(Annotation_Accession %in% gbk.annotation$Annotation_Accession)

singleton.proteins <- singleton.proteins %>%
    filter(Annotation_Accession %in% gbk.annotation$Annotation_Accession)
    
##########################################
## CRITICAL BUG TO FIX:
## there are 7,910 isolates that have annotated proteins in their Genbank
## annotation. 634 are missing protein (CDS) annotation in their genome.
## This difference is acceptable, but I need to fix the discrepancy in counts
## caused by the manual annotation having been
## constructed on an older version of the gbk_annotation.csv, which was
## missing some strains.
## for now, I have restricted these data by filtering on Annotation_Accession--
## see the temporary fixes in the code above.
###########################################################################
## Table 1. Isolates with antibiotic resistance genes.

## First column: count the number of isolates in each category.
isolate.totals <- gbk.annotation %>%
    group_by(Manual_Annotation) %>%
    summarize(total_isolates = n()) %>%
    arrange(desc(total_isolates))

## Second column: count the number of isolates with duplications in each category.
isolates.with.duplicate.genes <- duplicate.proteins %>%
    ## next two lines is to count isolates rather than genes
    select(Annotation_Accession, Manual_Annotation) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(isolates_with_duplicate_genes = n()) %>%
    arrange(desc(isolates_with_duplicate_genes))

## Third col: count the number of isolates with duplicated AR genes in each category.
AR.category.counts <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(Annotation_Accession, Manual_Annotation) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(isolates_with_duplicated_AR_genes = n()) %>%
    arrange(desc(isolates_with_duplicated_AR_genes))

## join columns to make Table 1 with raw data.
raw.Table1 <- isolate.totals %>% left_join(isolates.with.duplicate.genes) %>%
    left_join(AR.category.counts) %>%
    mutate(isolates_with_duplicated_AR_genes = replace_na(isolates_with_duplicated_AR_genes,0)) %>%
    arrange(desc(isolates_with_duplicated_AR_genes))

calc.expected.isolates.with.AR.genes <- function(raw.Table1) {
    total.isolates.with.duplicated.genes <- sum(raw.Table1$isolates_with_duplicate_genes)
    total.isolates.with.duplicated.AR.genes <- sum(raw.Table1$isolates_with_duplicated_AR_genes)
    Table <- raw.Table1 %>%
        mutate(expected_isolates_with_duplicated_AR_genes = total.isolates.with.duplicated.AR.genes * isolates_with_duplicate_genes/total.isolates.with.duplicated.genes)
    return(Table)
}

calc.isolate.AR.gene.enrichment.pvals <- function(raw.Table1) {
    
    total.isolates.with.duplicated.genes <- sum(raw.Table1$isolates_with_duplicate_genes)
    total.isolates.with.duplicated.AR.genes <- sum(raw.Table1$isolates_with_duplicated_AR_genes)

    Table <- raw.Table1 %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = isolates_with_duplicated_AR_genes,
                   n = total.isolates.with.duplicated.AR.genes,
                   p = isolates_with_duplicate_genes/total.isolates.with.duplicated.genes
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        ## use Benjamini-Hochberg p-value correction.
        select(-binom.test.pval) ## drop original p-value after the correction.
    
    return(Table)
}

## Add a fourth column: expected number of isolates with duplicated AR genes,
## based on the percentage of isolates with duplicated genes.
Table1 <- calc.expected.isolates.with.AR.genes(raw.Table1) %>%
    ## Add a fifth column: p-values for deviation from
    ## expected number of duplicated AR genes, using binomial test,
    ## correcting for multiple tests.
    calc.isolate.AR.gene.enrichment.pvals()

## write Table 1 to file.
write.csv(x=Table1,file="../results/AR-gene-duplication/Table1.csv")
#############
## Figure 1: visual representation of the result in Table 1.

## remove Unannotated data.
Fig1.data <- Table1 %>%
    filter(Manual_Annotation != "Unannotated") %>%
    mutate(Manual_Annotation = factor(
               Manual_Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                          "Anthropogenic-environment", "Food", "Freshwater",
                          "Agriculture", "Sediment", "Soil", "Plant-host",
                          "Marine","Terrestrial", "Fungal-host"))))

Fig1A <- ggplot(Fig1.data, aes(x=Manual_Annotation,y=isolates_with_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,1800) +
    ggtitle("Isolates with duplicated genes")

Fig1B <- ggplot(Fig1.data, aes(x=Manual_Annotation,y=expected_isolates_with_duplicated_AR_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,1800) +
    ggtitle("Expected distribution of isolates with duplicated ARGs")

Fig1C <- ggplot(Fig1.data, aes(x=Manual_Annotation,y=isolates_with_duplicated_AR_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,1800) +
    ggtitle("Observed distribution of isolates with duplicated ARGs")

Fig1 <- plot_grid(Fig1A,Fig1B,Fig1C,labels=c('A','B','C'),ncol=1)
ggsave("../results/AR-gene-duplication/Fig1.pdf",Fig1)
###########################################################################
## Positive control 1: Make a version of Table 1, examining the distribution
## of AR genes that have NOT duplicated.

## Second column:
## count the number of isolates with singleton AR genes in each category.

AR.singleton.category.counts <- singleton.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(Annotation_Accession, Manual_Annotation) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(isolates_with_singleton_AR_genes = n()) %>%
    arrange(desc(isolates_with_singleton_AR_genes))

calc.expected.isolates.with.singleton.AR.genes <- function(raw.ControlTable1) {
    summed.isolates <- sum(raw.ControlTable1$total_isolates)
    total.isolates.with.singleton.AR.genes <- sum(raw.ControlTable1$isolates_with_singleton_AR_genes)
    ControlTable <- raw.ControlTable1 %>%
        mutate(expected_isolates_with_singleton_AR_genes = total.isolates.with.singleton.AR.genes * total_isolates/summed.isolates)
    return(ControlTable)
}

calc.isolate.singleton.AR.gene.enrichment.pvals <- function(raw.ControlTable1) {
    
    summed.isolates <- sum(raw.ControlTable1$total_isolates)
    total.isolates.with.singleton.AR.genes <- sum(raw.ControlTable1$isolates_with_singleton_AR_genes)

    ControlTable <- raw.ControlTable1 %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = isolates_with_singleton_AR_genes,
                   n = total.isolates.with.singleton.AR.genes,
                   p = total_isolates/summed.isolates
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        ## use Benjamini-Hochberg p-value correction.
        select(-binom.test.pval) ## drop original p-value after the correction.
    
    return(ControlTable)
}

## join columns to make Control Table 1 with raw data.
raw.ControlTable1 <- isolate.totals %>%
    left_join(AR.singleton.category.counts) %>%
    mutate(isolates_with_singleton_AR_genes = replace_na(isolates_with_singleton_AR_genes,0)) %>%
    arrange(desc(isolates_with_singleton_AR_genes))


## Nice result! No categories are enriched with singleton AR genes,
## as most isolates have a gene that matches an antibiotic keyword.
## Animal-host isolates are depleted (perhaps due to aphid bacteria isolates?)
ControlTable1 <- calc.expected.isolates.with.singleton.AR.genes(raw.ControlTable1) %>%
    calc.isolate.singleton.AR.gene.enrichment.pvals()

## write ControlTable 1 to file.
write.csv(x=ControlTable1,file="../results/AR-gene-duplication/ControlTable1.csv")

######################################################################################
## Figure 2: visual representation of the result in Control Table 1.

## remove Unannotated data.
Fig2.data <- ControlTable1 %>%
    filter(Manual_Annotation != "Unannotated") %>%
    mutate(Manual_Annotation = factor(
               Manual_Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                          "Anthropogenic-environment", "Food", "Freshwater",
                          "Agriculture", "Sediment", "Soil", "Plant-host",
                          "Marine","Terrestrial", "Fungal-host"))))

Fig2A <- ggplot(Fig2.data, aes(x=Manual_Annotation,y=total_isolates)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,1900) +
    ggtitle("Total isolates")

Fig2B <- ggplot(Fig2.data, aes(x=Manual_Annotation,y=expected_isolates_with_singleton_AR_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,1900) +
    ggtitle("Expected distribution of isolates with singleton ARGs")

Fig2C <- ggplot(Fig2.data, aes(x=Manual_Annotation,y=isolates_with_singleton_AR_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,1900) +
    ggtitle("Observed distribution of isolates with singleton ARGs")

Fig2 <- plot_grid(Fig2A,Fig2B,Fig2C,labels=c('A','B','C'),ncol=1)
ggsave("../results/AR-gene-duplication/Fig2.pdf",Fig2)
####################################################################
## Table 2. Enrichment/deletion analysis of AR genes using duplicated genes,
## rather than number of isolates as in Table 1.

## First column: the number of duplicated genes in each category.
duplicate.genes.count <- duplicate.proteins %>%
    group_by(Manual_Annotation) %>%
    summarize(duplicate_genes = sum(count)) %>%
    arrange(desc(duplicate_genes))

## Second column: the number of duplicated MGE genes.
duplicate.MGE.genes.count <- duplicate.proteins %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(MGE_duplicates = sum(count)) %>%
    arrange(desc(MGE_duplicates))

## Third column: the number of duplicated AR genes.
duplicate.AR.genes.count <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(AR_duplicates = sum(count)) %>%
    arrange(desc(AR_duplicates))

## Fourth column: the expected number of duplicated AR genes in each category.
calc.expected.AR.duplicates <- function(raw.Table2) {
    total.duplicated.genes <- sum(raw.Table2$duplicate_genes)
    total.AR.duplicates <- sum(raw.Table2$AR_duplicates)
    Table <- raw.Table2 %>%
        mutate(expected_AR_duplicates = total.AR.duplicates * duplicate_genes/total.duplicated.genes)
    return(Table)
}

## Fifth column: p-value for enrichment/depletion of duplicated AR genes in each category.
calc.AR.duplicate.enrichment.pvals <- function(raw.Table2) {

    total.duplicated.genes <- sum(raw.Table2$duplicate_genes)
    total.AR.duplicates <- sum(raw.Table2$AR_duplicates)

    Table <- raw.Table2 %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = AR_duplicates,
                   n = total.AR.duplicates,
                   p = duplicate_genes/total.duplicated.genes
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        ## use Benjamini-Hochberg p-value correction.
        select(-binom.test.pval) ## drop original p-value after the correction.
    
    return(Table)
}

Table2 <- duplicate.genes.count %>% ## first column
    left_join(duplicate.MGE.genes.count) %>% ## second column
    left_join(duplicate.AR.genes.count) %>% ## third column
    mutate(AR_duplicates = replace_na(AR_duplicates, 0)) %>%
    arrange(desc(AR_duplicates)) %>%
    calc.expected.AR.duplicates() %>% ## fourth column
    calc.AR.duplicate.enrichment.pvals() ## fifth column

## write Table 2 to file.
write.csv(x=Table2,file="../results/AR-gene-duplication/Table2.csv")

############################################################
## Figures 3 and 4: visualization of Table 2 results.

## remove Unannotated data.
Fig3and4.data <- Table2 %>%
    filter(Manual_Annotation != "Unannotated") %>%
    mutate(Manual_Annotation = factor(
               Manual_Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                          "Anthropogenic-environment", "Food", "Freshwater",
                          "Agriculture", "Sediment", "Soil", "Plant-host",
                          "Marine","Terrestrial", "Fungal-host"))))

Fig3A <- ggplot(Fig3and4.data, aes(x=Manual_Annotation,y=duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,120000) +
    ggtitle("Total duplicate genes")

Fig3B <- ggplot(Fig3and4.data, aes(x=Manual_Annotation,y=MGE_duplicates)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,120000) +
    ggtitle("Total MGE duplicate genes")

Fig4A <- ggplot(Fig3and4.data, aes(x=Manual_Annotation,y=expected_AR_duplicates)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,700) +
    ggtitle("Expected distribution of duplicate ARGs")

Fig4B <- ggplot(Fig3and4.data, aes(x=Manual_Annotation,y=AR_duplicates)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,700) +
    ggtitle("Observed distribution of duplicate ARGs")

Fig3 <- plot_grid(Fig3A,Fig3B,labels=c('A','B'),ncol=1)
ggsave("../results/AR-gene-duplication/Fig3.pdf",Fig3)
Fig4 <- plot_grid(Fig4A,Fig4B,labels=c('A','B'),ncol=1)
ggsave("../results/AR-gene-duplication/Fig4.pdf",Fig4)

############

## Table 2B: the number of types of duplicate genes in each category,
## and the average num.

## Columns 1 and 2:
duplicated.gene.type.count <- duplicate.proteins %>%
    group_by(Manual_Annotation) %>%
    ## each row corresponds to a type of duplicated gene.
    summarize(duplicate_gene_types = n(),
              mean.duplicate.num = sum(count)/n()) %>%
    arrange(desc(duplicate_gene_types))


## Cols 3 & 4: the number of duplicated MGE gene types
duplicated.MGE.type.count <- duplicate.proteins %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Manual_Annotation) %>%
    ## each row corresponds to a type of duplicated MGE.
    summarize(duplicate_MGE_types = n(),
              mean.MGE.duplicate.num = sum(count)/n()) %>%
    arrange(desc(duplicate_MGE_types))

## Cols 5 & 6: the number of duplicated AR genes.
duplicated.ARG.type.count <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    ## each row corresponds to a type of duplicated ARG.
    summarize(duplicate_ARG_types = n(),
              mean.ARG.duplicate.num = sum(count)/n()) %>%
    arrange(desc(duplicate_ARG_types))


Table2B <- duplicated.gene.type.count %>% 
    left_join(duplicated.MGE.type.count) %>%
    left_join(duplicated.ARG.type.count) %>%
    mutate(duplicate_ARG_types = replace_na(duplicate_ARG_types, 0)) %>%
    mutate(mean.ARG.duplicate.num = replace_na(mean.ARG.duplicate.num, 0)) %>%
    arrange(desc(duplicate_ARG_types))

## write Table 2B to file.
write.csv(x=Table2B, file="../results/AR-gene-duplication/Table2B.csv")

################################################################################

## Positive control 2: Make a version of Table 2, examining the distribution
## of AR genes that have NOT duplicated.

## First column: the number of singleton genes in each category.
singleton.genes.count <- singleton.proteins %>%
    group_by(Manual_Annotation) %>%
    summarize(singleton_genes = sum(count)) %>%
    arrange(desc(singleton_genes))

## Second column: the number of singleton MGE genes.
singleton.MGE.genes.count <- singleton.proteins %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(MGE_singletons = sum(count)) %>%
    arrange(desc(MGE_singletons))

## Third column: the number of singleton AR genes.
singleton.AR.genes.count <- singleton.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(AR_singletons = sum(count)) %>%
    arrange(desc(AR_singletons))

## Fourth column: the expected number of singleton AR genes in each category.
calc.expected.AR.singletons <- function(raw.ControlTable2) {
    total.singleton.genes <- sum(raw.ControlTable2$singleton_genes)
    total.AR.singletons <- sum(raw.ControlTable2$AR_singletons)
    ControlTable <- raw.ControlTable2 %>%
        mutate(expected_AR_singletons = total.AR.singletons * singleton_genes/total.singleton.genes)
    return(ControlTable)
}

## Fifth column: p-value for enrichment/depletion of singleton AR genes
## in each category.
calc.AR.singleton.enrichment.pvals <- function(raw.ControlTable2) {

    total.singleton.genes <- sum(raw.ControlTable2$singleton_genes)
    total.AR.singletons <- sum(raw.ControlTable2$AR_singletons)

    ControlTable <- raw.ControlTable2 %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = AR_singletons,
                   n = total.AR.singletons,
                   p = singleton_genes/total.singleton.genes
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        ## use Benjamini-Hochberg p-value correction.
        select(-binom.test.pval) ## drop original p-value after the correction.
    
    return(ControlTable)
}

## IMPORTANT: Unannotated strains are ridiculously enriched in singleton AR genes.
## the null distribution changes significantly depending on whether
## Unannotated strains are included in this table, or not.

## when including Unannotated strains, human isolates are highly depleted in
## singleton AR genes; when removing them, human isolates are significantly
## enriched in AR genes.

## In both cases, however, soil is significantly enriched in AR singletons.
## This shows that soil has more singleton AR genes than human isolates, no matter
## how this analysis is done.

## I will present the table that exclude the Unannotated strains, and I may
## put the table with the Unannotated strains included in the supplement.

raw.ControlTable2 <- singleton.genes.count %>% ## first column
    left_join(singleton.MGE.genes.count) %>% ## second column
    left_join(singleton.AR.genes.count) %>% ## third column
    mutate(AR_singletons = replace_na(AR_singletons, 0)) %>%
    arrange(desc(AR_singletons))

## CRITICAL STEP: keep Unannotated strains.
ControlTable2A <- raw.ControlTable2 %>%
    calc.expected.AR.singletons() %>% ## fourth column
    calc.AR.singleton.enrichment.pvals() ## fifth column

## CRITICAL STEP: remove Unannotated strains.
ControlTable2B <- raw.ControlTable2 %>%
    filter(Manual_Annotation != "Unannotated") %>%
    calc.expected.AR.singletons() %>% ## fourth column
    calc.AR.singleton.enrichment.pvals() ## fifth column

write.csv(x=ControlTable2A, file="../results/AR-gene-duplication/ControlTable2A.csv")
write.csv(x=ControlTable2B, file="../results/AR-gene-duplication/ControlTable2B.csv")

################################################################################
## Figures 5 and 6: visualization of control table 2B.

Fig5and6.data <- ControlTable2B %>%
    mutate(Manual_Annotation = factor(
               Manual_Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                          "Anthropogenic-environment", "Food", "Freshwater",
                          "Agriculture", "Sediment", "Soil", "Plant-host",
                          "Marine","Terrestrial", "Fungal-host"))))

Fig5A <- ggplot(Fig5and6.data, aes(x=Manual_Annotation,y=singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,7500000) +
    ggtitle("Total singleton genes")

Fig5B <- ggplot(Fig5and6.data, aes(x=Manual_Annotation,y=MGE_singletons)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,7500000) +
    ggtitle("Total MGE singleton genes")

Fig6A <- ggplot(Fig5and6.data, aes(x=Manual_Annotation,y=expected_AR_singletons)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  +
    ggtitle("Expected distribution of singleton ARGs")

Fig6B <- ggplot(Fig5and6.data, aes(x=Manual_Annotation,y=AR_singletons)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  +
    ggtitle("Observed distribution of singleton ARGs")

Fig5 <- plot_grid(Fig5A,Fig5B,labels=c('A','B'),ncol=1)
ggsave("../results/AR-gene-duplication/Fig5.pdf",Fig5)
Fig6 <- plot_grid(Fig6A,Fig6B,labels=c('A','B'),ncol=1)
ggsave("../results/AR-gene-duplication/Fig6.pdf",Fig6)
################################################################################

## Table 3. Show number of duplicated genes on chromosomes, and number of
## duplicate genes on plasmids, for each category, for duplicated genes
## and duplicated AR genes. This table of raw data goes into the text. Then,
## sum over all categories for a 2x2 contingency table and report the result of a
## Fisher's exact test for asssociation between duplicated AR genes and plasmids.

## Column 1
duplicate.chromosome.genes.count <- duplicate.proteins %>%
    group_by(Manual_Annotation) %>%
    summarize(chromosomal_duplicate_genes = sum(chromosome_count))

## Column 2
duplicate.plasmid.genes.count <- duplicate.proteins %>%
    group_by(Manual_Annotation) %>%
    summarize(plasmid_duplicate_genes = sum(plasmid_count))

## Column 3
duplicate.AR.chromosome.genes.count <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(chromosomal_AR_duplicate_genes = sum(chromosome_count))

## Column 4
duplicate.AR.plasmid.genes.count <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(plasmid_AR_duplicate_genes = sum(plasmid_count))

Table3 <- duplicate.chromosome.genes.count %>%
    left_join(duplicate.plasmid.genes.count) %>%
    left_join(duplicate.AR.chromosome.genes.count) %>%
    mutate(chromosomal_AR_duplicate_genes=replace_na(chromosomal_AR_duplicate_genes, 0)) %>%
    left_join(duplicate.AR.plasmid.genes.count) %>%
    mutate(plasmid_AR_duplicate_genes = replace_na(plasmid_AR_duplicate_genes, 0)) %>%
    arrange(desc(plasmid_AR_duplicate_genes))

## write Table 3 to file.
write.csv(x=Table3,file="../results/AR-gene-duplication/Table3.csv")


## get values for Fisher's exact test.
total.chr.AR.duplicates <- sum(Table3$chromosomal_AR_duplicate_genes)
total.plasmid.AR.duplicates <- sum(Table3$plasmid_AR_duplicate_genes)

total.chr.duplicates <- sum(Table3$chromosomal_duplicate_genes)
total.plasmid.duplicates <- sum(Table3$plasmid_duplicate_genes)

total.nonAR.chr.duplicates <- total.chr.duplicates - total.chr.AR.duplicates
total.nonAR.plasmid.duplicates <- total.plasmid.duplicates - total.plasmid.AR.duplicates

Table4.contingency.table <- matrix(c(total.chr.AR.duplicates,
                                     total.plasmid.AR.duplicates,
                                     total.nonAR.chr.duplicates,
                                     total.nonAR.plasmid.duplicates),nrow=2)
## label the rows and columns of the contingency table.
rownames(Table4.contingency.table) <- c("chromosome","plasmid")
colnames(Table4.contingency.table) <- c("AR duplicate genes","non-AR duplicate genes")

## write Table 4 contingency table to file.
write.csv(x=Table4.contingency.table,file="../results/AR-gene-duplication/Table4.csv")

## p < 1e-197
fisher.test(Table4.contingency.table)
fisher.test(Table4.contingency.table)$p.value
################################################################################
## Figure 7. Visualization of the data in Table 3.

Fig7.data <- Table3 %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    mutate(Manual_Annotation = factor(
               Manual_Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                          "Anthropogenic-environment", "Food", "Freshwater",
                          "Agriculture", "Sediment", "Soil", "Plant-host",
                          "Marine","Terrestrial", "Fungal-host"))))

Fig7A <- ggplot(Fig7.data, aes(x=Manual_Annotation,y=chromosomal_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,90000) +
    ggtitle("Chromosomal duplicate genes")

Fig7B <- ggplot(Fig7.data, aes(x=Manual_Annotation,y=plasmid_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,90000) +
    ggtitle("Plasmid duplicate genes")

Fig7C <- ggplot(Fig7.data, aes(x=Manual_Annotation,y=chromosomal_AR_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,500) +
    ggtitle("Chromosomal duplicate ARGs")

Fig7D <- ggplot(Fig7.data, aes(x=Manual_Annotation,y=plasmid_AR_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,500) +
    ggtitle("Plasmid duplicate ARGs")

Fig7 <- plot_grid(Fig7A, Fig7B, Fig7C, Fig7D,labels=c('A','B','C','D'),ncol=2)
ggsave("../results/AR-gene-duplication/Fig7.pdf",Fig7)
################################################################################
## Positive control 3: look at distribution of singleton AR genes on
## chromosomes and plasmids, to compare with Tables 3 and 4.

## Column 1
singleton.chromosome.genes.count <- singleton.proteins %>%
    group_by(Manual_Annotation) %>%
    summarize(chromosomal_singleton_genes = sum(chromosome_count))

## Column 2
singleton.plasmid.genes.count <- singleton.proteins %>%
    group_by(Manual_Annotation) %>%
    summarize(plasmid_singleton_genes = sum(plasmid_count))

## Column 3
singleton.AR.chromosome.genes.count <- singleton.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(chromosomal_AR_singleton_genes = sum(chromosome_count))

## Column 4
singleton.AR.plasmid.genes.count <- singleton.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(plasmid_AR_singleton_genes = sum(plasmid_count))

ControlTable3 <- singleton.chromosome.genes.count %>%
    left_join(singleton.plasmid.genes.count) %>%
    left_join(singleton.AR.chromosome.genes.count) %>%
    mutate(chromosomal_AR_singleton_genes=replace_na(chromosomal_AR_singleton_genes, 0)) %>%
    left_join(singleton.AR.plasmid.genes.count) %>%
    mutate(plasmid_AR_singleton_genes = replace_na(plasmid_AR_singleton_genes, 0)) %>%
    arrange(desc(plasmid_AR_singleton_genes))

## get values for Fisher's exact test.
total.chr.AR.singletons <- sum(ControlTable3$chromosomal_AR_singleton_genes)
total.plasmid.AR.singletons <- sum(ControlTable3$plasmid_AR_singleton_genes)

total.chr.singletons <- sum(ControlTable3$chromosomal_singleton_genes)
total.plasmid.singletons <- sum(ControlTable3$plasmid_singleton_genes)

total.nonAR.chr.singletons <- total.chr.singletons - total.chr.AR.singletons
total.nonAR.plasmid.singletons <- total.plasmid.singletons - total.plasmid.AR.singletons

ControlTable4.contingency.table <- matrix(c(total.chr.AR.singletons,
                                     total.plasmid.AR.singletons,
                                     total.nonAR.chr.singletons,
                                     total.nonAR.plasmid.singletons),nrow=2)
## label the rows and columns of the contingency table.
rownames(ControlTable4.contingency.table) <- c("chromosome","plasmid")
colnames(ControlTable4.contingency.table) <- c("AR singleton genes","non-AR singleton genes")

## This positive control shows singleton AR genes are highly enriched on plasmids,
## based on a comparison with the distribution of singleton genes overall.
## Therefore AR genes are generally associated with plasmids, regardless of
## status of being a duplication or not.

## This does NOT invalidate the main result of this analysis, that duplicate AR
## genes are more enriched on plasmids in comparison to the distribution of
## duplicate genes overall.

fisher.test(ControlTable4.contingency.table)
fisher.test(ControlTable4.contingency.table)$p.value

################################################################################
## Figure 8: visualization of result in ControlTable3

Fig8.data <- ControlTable3 %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    mutate(Manual_Annotation = factor(
               Manual_Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                          "Anthropogenic-environment", "Food", "Freshwater",
                          "Agriculture", "Sediment", "Soil", "Plant-host",
                          "Marine","Terrestrial", "Fungal-host"))))

Fig8A <- ggplot(Fig8.data, aes(x=Manual_Annotation,y=chromosomal_singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,7000000) +
    ggtitle("Chromosomal singleton genes")

Fig8B <- ggplot(Fig8.data, aes(x=Manual_Annotation,y=plasmid_singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,7000000) +
    ggtitle("Plasmid singleton genes")

Fig8C <- ggplot(Fig8.data, aes(x=Manual_Annotation,y=chromosomal_AR_singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,9000) +
    ggtitle("Chromosomal singleton ARGs")

Fig8D <- ggplot(Fig8.data, aes(x=Manual_Annotation,y=plasmid_AR_singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,9000) +
    ggtitle("Plasmid singleton ARGs")

Fig8 <- plot_grid(Fig8A, Fig8B, Fig8C, Fig8D,labels=c('A','B','C','D'),ncol=2)
ggsave("../results/AR-gene-duplication/Fig8.pdf",Fig8)

#######################################################################
## Make a figure for grant with Yi.
## Combine and reshape the data in Table 3 and Control Table 3.
grant.df <- Table3 %>%
    full_join(ControlTable3) %>%
    mutate(chromosomal_genes=chromosomal_duplicate_genes+chromosomal_singleton_genes) %>%
    select(-chromosomal_duplicate_genes,-chromosomal_singleton_genes) %>%
    mutate(plasmid_genes=plasmid_duplicate_genes+plasmid_singleton_genes) %>%
    select(-plasmid_duplicate_genes,-plasmid_singleton_genes) %>%
    mutate(chromosomal_AR_genes=chromosomal_AR_duplicate_genes+chromosomal_AR_singleton_genes) %>%
    select(-chromosomal_AR_duplicate_genes,-chromosomal_AR_singleton_genes) %>%
    mutate(plasmid_AR_genes=plasmid_AR_duplicate_genes+plasmid_AR_singleton_genes) %>%
    select(-plasmid_AR_duplicate_genes,-plasmid_AR_singleton_genes) %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    mutate(Manual_Annotation = factor(
               Manual_Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                              "Anthropogenic-environment", "Food", "Freshwater",
                              "Agriculture", "Sediment", "Soil", "Plant-host",
                              "Marine","Terrestrial", "Fungal-host"))))

grantFigA <- ggplot(grant.df, aes(x=Manual_Annotation,y=chromosomal_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") +
    ggtitle("Chromosomal genes")

grantFigB <- ggplot(grant.df, aes(x=Manual_Annotation,y=plasmid_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,7000000) +
    ggtitle("Plasmid genes")

grantFigC <- ggplot(grant.df, aes(x=Manual_Annotation,y=chromosomal_AR_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,9000) +
    ggtitle("Chromosomal AR genes")

grantFigD <- ggplot(grant.df, aes(x=Manual_Annotation,y=plasmid_AR_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,9000) +
    ggtitle("Plasmid AR genes")

grantFig <- plot_grid(grantFigA, grantFigB, grantFigC, grantFigD,labels=c('A','B','C','D'),ncol=2)
ggsave("../results/AR-gene-duplication/grantFig.pdf",grantFig)


## get values for Fisher's exact test.
total.chr.ARGs <- sum(grant.df$chromosomal_AR_genes)
total.plasmid.ARGs <- sum(grant.df$plasmid_AR_genes)

total.chr.genes <- sum(grant.df$chromosomal_genes)
total.plasmid.genes <- sum(grant.df$plasmid_genes)

total.nonAR.chr.genes <- total.chr.genes - total.chr.ARGs
total.nonAR.plasmid.genes <- total.plasmid.genes - total.plasmid.ARGs

grant.contingency.table <- matrix(c(total.chr.ARGs,
                                     total.plasmid.ARGs,
                                     total.nonAR.chr.genes,
                                     total.nonAR.plasmid.genes),nrow=2)
## label the rows and columns of the contingency table.
rownames(grant.contingency.table) <- c("chromosome","plasmid")
colnames(grant.contingency.table) <- c("AR genes","non-AR genes")

## This positive control shows singleton AR genes are highly enriched on plasmids,
## based on a comparison with the distribution of singleton genes overall.
## Therefore AR genes are generally associated with plasmids, regardless of
## status of being a duplication or not.

## This does NOT invalidate the main result of this analysis, that duplicate AR
## genes are more enriched on plasmids in comparison to the distribution of
## duplicate genes overall.

fisher.test(grant.contingency.table)
fisher.test(grant.contingency.table)$p.value


################################################################################
## Analysis of duplicate pairs found just on chromosome, just on plasmid, or
## on both chromosomes and plasmids.

## let's look at cases of identical sequences on chromosomes and plasmids.

both.chr.and.plasmid.cases <- duplicate.proteins %>%
    filter(chromosome_count >= 1 & plasmid_count >= 1) %>%
    arrange(desc(count))

just.chromosome.cases <- duplicate.proteins %>%
    filter(chromosome_count >= 1 & plasmid_count == 0) %>%
    arrange(desc(count))

just.plasmid.cases <- duplicate.proteins %>%
    filter(chromosome_count == 0 & plasmid_count >= 1) %>%
    arrange(desc(count))

################################################################################

## simple HGT analysis.

## get the match between sequences and product annotations.
duplicate.protein.annotations <- duplicate.proteins %>%
    select(product, sequence) %>%
    distinct()

HGT.candidates <- duplicate.proteins %>%
    group_by(sequence) %>%
    summarize(number.of.genomes = n()) %>%
    filter(number.of.genomes > 1) %>%
    arrange(desc(number.of.genomes)) %>%
    left_join(duplicate.protein.annotations)

HGT.candidate.summary <- HGT.candidates %>%
    select(-sequence)

MGE.HGT.candidates <- HGT.candidate.summary %>%
    filter(str_detect(.$product,IS.keywords))

non.MGE.HGT.candidates <- HGT.candidate.summary %>%
    filter(!str_detect(.$product,IS.keywords))

AR.HGT.candidates <- HGT.candidate.summary %>%
    filter(str_detect(.$product,antibiotic.keywords))
