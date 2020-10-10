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
library(ggrepel)
library(data.table)

## annotate source sequence as plasmid or chromosome.
genome.database <- read.csv("../results/AR-gene-duplication/chromosome-plasmid-table.csv")

raw.gbk.annotation <- as_tibble(read.csv("../results/AR-gene-duplication/gbk-annotation-table.csv")) %>% ## get species name annotation from genome.database.
    left_join(genome.database)

length(unique(genome.database$Annotation_Accession))

## This uses the manually-annotated data that I have so far.
## IMPORTANT TODO: update manual annotation version,
## with a script that generates this file.
gbk.annotation <- as_tibble(read.csv("../data/manually-curated-gbk-annotation-table.csv")) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Manual_Annotation = replace_na(Manual_Annotation,"Unannotated")) %>%
    ## get species name annotation from genome.database
    left_join(genome.database)

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
IS.keywords <- "IS|transposon|Transposase|transposase|hypothetical protein|Phage|phage|integrase|Integrase|tail|intron|Mobile|mobile|antitoxin|toxin|capsid|plasmid|Plasmid"

## now look at a few antibiotic-specific annotations.
antibiotic.keywords <- "lactamase|chloramphenicol|quinolone|antibiotic resistance|tetracycline|VanZ"

## Potential TODO: Use the same regular expressions used by Zeevi et al. (2019).
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
## call garbage collector to free up memory.
gc()
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
## Figure 1 A & B: Diagram of the analysis workflow, made in Inkscape/Illustrator.
###########################################################################
## Supplementary Table 1. Isolates with antibiotic resistance genes.

## First column: count the number of isolates in each category.
isolate.totals <- gbk.annotation %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    group_by(Manual_Annotation) %>%
    summarize(total_isolates = n()) %>%
    arrange(desc(total_isolates))

## Second column: count the number of isolates with duplications in each category.
isolates.with.duplicate.genes <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
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
raw.TableS1 <- isolate.totals %>% left_join(isolates.with.duplicate.genes) %>%
    left_join(AR.category.counts) %>%
    mutate(isolates_with_duplicated_AR_genes = replace_na(isolates_with_duplicated_AR_genes,0)) %>%
    arrange(desc(isolates_with_duplicated_AR_genes))

calc.expected.isolates.with.AR.genes <- function(raw.TableS1) {
    total.isolates.with.duplicated.genes <- sum(raw.TableS1$isolates_with_duplicate_genes)
    total.isolates.with.duplicated.AR.genes <- sum(raw.TableS1$isolates_with_duplicated_AR_genes)
    Table <- raw.TableS1 %>%
        mutate(expected_isolates_with_duplicated_AR_genes = total.isolates.with.duplicated.AR.genes * isolates_with_duplicate_genes/total.isolates.with.duplicated.genes)
    return(Table)
}

calc.isolate.AR.gene.enrichment.pvals <- function(raw.TableS1) {
    
    total.isolates.with.duplicated.genes <- sum(raw.TableS1$isolates_with_duplicate_genes)
    total.isolates.with.duplicated.AR.genes <- sum(raw.TableS1$isolates_with_duplicated_AR_genes)

    Table <- raw.TableS1 %>%
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
TableS1 <- calc.expected.isolates.with.AR.genes(raw.TableS1) %>%
    ## Add a fifth column: p-values for deviation from
    ## expected number of duplicated AR genes, using binomial test,
    ## correcting for multiple tests.
    calc.isolate.AR.gene.enrichment.pvals()

## write Table 1 to file.
write.csv(x=TableS1,file="../results/AR-gene-duplication/TableS1.csv")
####################################################
## Figure 2: visual representation of the result in Table S1.

## helper function to wrap some common code that massages data that
## will go into making figures.
table.df.to.figure.df <- function(table.df) {
    table.df %>%
        ## remove Unannotated data.
        filter(Manual_Annotation != "Unannotated") %>%
        mutate(Manual_Annotation = factor(
                   Manual_Annotation,
                   levels = rev(c("Human-host","Livestock","Animal-host",
                                  "Anthropogenic-environment", "Food", "Freshwater",
                                  "Agriculture", "Sediment", "Soil", "Plant-host",
                                  "Marine","Terrestrial", "Fungal-host"))))
}

Fig2A <- ggplot(TableS1, aes(x=log2(isolates_with_duplicate_genes),
                              y=log2(isolates_with_duplicated_AR_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(isolates_with_duplicate_genes),
                     y = log2(isolates_with_duplicated_AR_genes),
                     xend = log2(isolates_with_duplicate_genes),
                     yend = log2(expected_isolates_with_duplicated_AR_genes)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_isolates_with_duplicated_AR_genes)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](Isolates~with~duplicated~genes))) +
    ylab(expression(log[2](Isolates~with~duplicated~ARGs)))    

###########################################################################
## Positive control 1: Make a version of Table S1, examining the distribution
## of AR genes that have NOT duplicated.

## Second column:
## count the number of isolates with singleton AR genes in each category.

AR.singleton.category.counts <- singleton.proteins %>%
    ## remove Unannotated data.
    filter(Manual_Annotation != "Unannotated") %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(Annotation_Accession, Manual_Annotation) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(isolates_with_singleton_AR_genes = n()) %>%
    arrange(desc(isolates_with_singleton_AR_genes))
gc() ## free memory.

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

##############################################################
## Figure 1C: main figure, showing enrichment of AR duplicates
## in human hosts and livestock.

make.Fig1C.df <- function(TableS1, TableS3, ControlTable1, ControlTable3) {
    ## join duplicate and singleton tables to make Fig 1C.

    ## have to remove p-values from the two tables, because
    ## the column names are the same, but the values are different
    ## (because these are two different tests).
    no.pval.TableS1 <- select(TableS1, -corrected.pval)
    no.pval.ControlTable1 <- select(ControlTable1, -corrected.pval)
    
    Fig1C.df <- no.pval.TableS1 %>%
        full_join(no.pval.ControlTable1) %>%
        full_join(TableS3) %>%
        full_join(ControlTable3)
    
    total_isolates.sum <- sum(Fig1C.df$total_isolates)
    isolates_with_duplicate_genes.sum <- sum(Fig1C.df$isolates_with_duplicate_genes)
    isolates_with_singleton_AR_genes.sum <- sum(Fig1C.df$isolates_with_singleton_AR_genes)
    isolates_with_duplicated_AR_genes.sum <- sum(Fig1C.df$isolates_with_duplicated_AR_genes)

    Fig1C.df <- Fig1C.df %>%
        ## calculate y-coordinates for line for duplicate genes.
        mutate(yvals.for.isolates_with_duplicate_genes.line = total_isolates * isolates_with_duplicate_genes.sum/total_isolates.sum) %>%
        ## calculate y-coordinates for line for singleton AR genes.
        mutate(yvals.for.isolates_with_singleton_AR_genes.line = total_isolates * isolates_with_singleton_AR_genes.sum/total_isolates.sum) %>%
        ## calculate y-coordinates for line for duplicate AR genes.
        mutate(yvals.for.isolates_with_duplicated_AR_genes.line = total_isolates * isolates_with_duplicated_AR_genes.sum/total_isolates.sum) %>%
        ## calculate the percentage of genes on plasmids for symbol size.
        mutate(plasmid_duplicate_percent = plasmid_duplicate_genes/(plasmid_duplicate_genes+chromosomal_duplicate_genes+0.1)) %>%
        mutate(plasmid_AR_singleton_percent = plasmid_AR_singleton_genes/(plasmid_AR_singleton_genes+chromosomal_AR_singleton_genes+0.1)) %>%
        mutate(plasmid_AR_duplicate_percent = plasmid_AR_duplicate_genes/(plasmid_AR_duplicate_genes+chromosomal_AR_duplicate_genes+0.1))

    return(Fig1C.df)
}

Fig1C.df <- make.Fig1C.df(TableS1, TableS3, ControlTable1, ControlTable3)

make.Fig1C <- function(Fig1C.df) {

    total_isolates.sum <- sum(Fig1C.df$total_isolates)
    isolates_with_duplicated_AR_genes.sum <- sum(Fig1C.df$isolates_with_duplicated_AR_genes)
    isolates_with_singleton_AR_genes.sum <- sum(Fig1C.df$isolates_with_singleton_AR_genes)
    isolates_with_duplicate_genes.sum <- sum(Fig1C.df$isolates_with_duplicate_genes)

    Fig1C.color.palette <- scales::viridis_pal()(3)

    Fig1C <- ggplot(Fig1C.df, aes(x=total_isolates,
                                  y=isolates_with_duplicated_AR_genes,
                                  label=Manual_Annotation)) +
        theme_classic() +
        geom_point(aes(size = plasmid_AR_duplicate_percent * 0.5),
                   color=Fig1C.color.palette[1], alpha=0.2) +
        geom_point(aes(y=isolates_with_singleton_AR_genes,
                       size=plasmid_AR_singleton_percent * 0.5),
                   color=Fig1C.color.palette[2],alpha=0.2) +
        geom_point(aes(y=isolates_with_duplicate_genes,
                       size=plasmid_duplicate_percent * 0.5),color="gray",alpha=0.2) +
        geom_line(aes(y=yvals.for.isolates_with_duplicated_AR_genes.line),
                  color=Fig1C.color.palette[1]) +
        geom_line(aes(y=yvals.for.isolates_with_singleton_AR_genes.line),
                  color=Fig1C.color.palette[2]) +
        geom_line(aes(y=yvals.for.isolates_with_duplicate_genes.line),
                  color="gray") +
        geom_text_repel(size=2.5) +
        scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        xlab("Total Isolates") +
        ylab("Isolates in given class") +
        annotate("text", x = 90, y = 2.2, label = "Isolates with duplicated ARGs",
                 angle = 32.4, color = Fig1C.color.palette[1],size=3) +
        annotate("text", x = 90, y = 18, label = "Isolates with singleton ARGs",
                 angle = 32.4, color = Fig1C.color.palette[2],size=3) +
        annotate("text", x = 90, y = 28, label = "Isolates with duplicated genes",
                 angle = 32.4, color = "gray",size=3) +
        guides(size=FALSE) 
        
    return(Fig1C)
}

Fig1C <- make.Fig1C(Fig1C.df)
ggsave(Fig1C,file="../results/AR-gene-duplication/Fig1C.pdf")
######################################################################################
## Fig 2B : visual representation of the result in Control Table 1.

Fig2B <- ggplot(ControlTable1, aes(x=log2(total_isolates),
                              y=log2(isolates_with_singleton_AR_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(total_isolates),
                     y = log2(isolates_with_singleton_AR_genes),
                     xend = log2(total_isolates),
                     yend = log2(expected_isolates_with_singleton_AR_genes)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_isolates_with_singleton_AR_genes)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](Total~isolates))) +
    ylab(expression(log[2](Isolates~with~singleton~ARGs)))    

####################################################################
## Supplementary Table S2. Enrichment/deletion analysis of AR genes using duplicated genes,
## rather than number of isolates as in Table 1.

## First column: the number of duplicated genes in each category.
duplicate.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    group_by(Manual_Annotation) %>%
    summarize(duplicate_genes = sum(count)) %>%
    arrange(desc(duplicate_genes))

## Second column: the number of duplicated MGE genes.
duplicate.MGE.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(MGE_duplicates = sum(count)) %>%
    arrange(desc(MGE_duplicates))

## Third column: the number of duplicated AR genes.
duplicate.AR.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(AR_duplicates = sum(count)) %>%
    arrange(desc(AR_duplicates))

## Fourth column: the expected number of duplicated AR genes in each category.
calc.expected.AR.duplicates <- function(raw.TableS2) {
    total.duplicated.genes <- sum(raw.TableS2$duplicate_genes)
    total.AR.duplicates <- sum(raw.TableS2$AR_duplicates)
    Table <- raw.TableS2 %>%
        mutate(expected_AR_duplicates = total.AR.duplicates * duplicate_genes/total.duplicated.genes)
    return(Table)
}

## Fifth column: p-value for enrichment/depletion of duplicated AR genes in each category.
calc.AR.duplicate.enrichment.pvals <- function(raw.TableS2) {

    total.duplicated.genes <- sum(raw.TableS2$duplicate_genes)
    total.AR.duplicates <- sum(raw.TableS2$AR_duplicates)

    Table <- raw.TableS2 %>%
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

TableS2 <- duplicate.genes.count %>% ## first column
    left_join(duplicate.MGE.genes.count) %>% ## second column
    left_join(duplicate.AR.genes.count) %>% ## third column
    mutate(AR_duplicates = replace_na(AR_duplicates, 0)) %>%
    arrange(desc(AR_duplicates)) %>%
    calc.expected.AR.duplicates() %>% ## fourth column
    calc.AR.duplicate.enrichment.pvals() ## fifth column

## write Table 2 to file.
write.csv(x=TableS2,file="../results/AR-gene-duplication/TableS2.csv")

############################################################
## S2 Figure and Figure 2C: visualization of Table S2 results.

S2Fig <- ggplot(TableS2, aes(x=log2(duplicate_genes),
                              y=log2(MGE_duplicates),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_point(color='red') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](Duplicated~genes))) +
    ylab(expression(log[2](Duplicated~MGE~genes)))    

ggsave("../results/AR-gene-duplication/S2Fig.pdf",S2Fig)

Fig2C <- ggplot(TableS2, aes(x=log2(duplicate_genes),
                              y=log2(AR_duplicates),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(duplicate_genes),
                     y = log2(AR_duplicates),
                     xend = log2(duplicate_genes),
                     yend = log2(expected_AR_duplicates)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_AR_duplicates)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](Duplicated~genes))) +
    ylab(expression(log[2](Duplicated~ARGs)))    

############

## Table S2B: the number of types of duplicate genes in each category,
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

TableS2B <- duplicated.gene.type.count %>% 
    left_join(duplicated.MGE.type.count) %>%
    left_join(duplicated.ARG.type.count) %>%
    mutate(duplicate_ARG_types = replace_na(duplicate_ARG_types, 0)) %>%
    mutate(mean.ARG.duplicate.num = replace_na(mean.ARG.duplicate.num, 0)) %>%
    arrange(desc(duplicate_ARG_types))

## write Table S2B to file.
write.csv(x=TableS2B, file="../results/AR-gene-duplication/TableS2B.csv")

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
## S4 and S5 Figures: visualization of control table 2B.

S4Fig <- ggplot(ControlTable2B, aes(x=log2(singleton_genes),
                              y=log2(MGE_singletons),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(singleton_genes),
                     y = log2(MGE_singletons),
                     xend = log2(singleton_genes),
                     yend = log2(MGE_singletons)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='black') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](Singleton~genes))) +
    ylab(expression(log[2](Singleton~MGE~genes)))    

S4Fig <- plot_grid(S4FigA,S4FigB,labels=c('A','B'),ncol=1)
ggsave("../results/AR-gene-duplication/S4Fig.pdf",S4Fig)


Fig2D <- ggplot(ControlTable2B, aes(x=log2(singleton_genes),
                              y=log2(AR_singletons),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(singleton_genes),
                     y = log2(AR_singletons),
                     xend = log2(singleton_genes),
                     yend = log2(expected_AR_singletons)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_AR_singletons)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](Singleton~genes))) +
    ylab(expression(log[2](Singleton~ARGs)))    

### NOW *FINALLY* MAKE FIGURE 2.
Fig2 <- plot_grid(Fig2A,Fig2B,Fig2C,Fig2D,labels=c('A','B','C','D'),nrow=2)
ggsave("../results/AR-gene-duplication/Fig2.pdf",Fig2)

################################################################################

## Table S3. Show number of duplicated genes on chromosomes, and number of
## duplicate genes on plasmids, for each category, for duplicated genes
## and duplicated AR genes. This table of raw data goes into the text. Then,
## sum over all categories for a 2x2 contingency table and report the result of a
## Fisher's exact test for asssociation between duplicated AR genes and plasmids.

## Column 1
duplicate.chromosome.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    group_by(Manual_Annotation) %>%
    summarize(chromosomal_duplicate_genes = sum(chromosome_count))

## Column 2
duplicate.plasmid.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    group_by(Manual_Annotation) %>%
    summarize(plasmid_duplicate_genes = sum(plasmid_count))

## Column 3
duplicate.AR.chromosome.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(chromosomal_AR_duplicate_genes = sum(chromosome_count))

## Column 4
duplicate.AR.plasmid.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(plasmid_AR_duplicate_genes = sum(plasmid_count))

TableS3 <- duplicate.chromosome.genes.count %>%
    left_join(duplicate.plasmid.genes.count) %>%
    left_join(duplicate.AR.chromosome.genes.count) %>%
    mutate(chromosomal_AR_duplicate_genes=replace_na(chromosomal_AR_duplicate_genes, 0)) %>%
    left_join(duplicate.AR.plasmid.genes.count) %>%
    mutate(plasmid_AR_duplicate_genes = replace_na(plasmid_AR_duplicate_genes, 0)) %>%
    arrange(desc(plasmid_AR_duplicate_genes))

## write Table S3 to file.
write.csv(x=TableS3,file="../results/AR-gene-duplication/TableS3.csv")


## get values for Fisher's exact test.
total.chr.AR.duplicates <- sum(TableS3$chromosomal_AR_duplicate_genes)
total.plasmid.AR.duplicates <- sum(TableS3$plasmid_AR_duplicate_genes)

total.chr.duplicates <- sum(TableS3$chromosomal_duplicate_genes)
total.plasmid.duplicates <- sum(TableS3$plasmid_duplicate_genes)

total.nonAR.chr.duplicates <- total.chr.duplicates - total.chr.AR.duplicates
total.nonAR.plasmid.duplicates <- total.plasmid.duplicates - total.plasmid.AR.duplicates

TableS4.contingency.table <- matrix(c(total.chr.AR.duplicates,
                                     total.plasmid.AR.duplicates,
                                     total.nonAR.chr.duplicates,
                                     total.nonAR.plasmid.duplicates),nrow=2)
## label the rows and columns of the contingency table.
rownames(TableS4.contingency.table) <- c("chromosome","plasmid")
colnames(TableS4.contingency.table) <- c("AR duplicate genes","non-AR duplicate genes")

## write Table S4 contingency table to file.
write.csv(x=TableS4.contingency.table,file="../results/AR-gene-duplication/Table4.csv")

## p < 1e-197
fisher.test(TableS4.contingency.table)
fisher.test(TableS4.contingency.table)$p.value
################################################################################
## Figure 3. Visualization of the data in Table S3.

## SEE CODE FOR SUPPLEMENTARY FIGURE S6 FOR PANELS B AND D!!!
## TODO: REFACTOR CODE TO IMPROVE ORGANIZATION!

## Panel 3A.
## add a column for expected number of duplicate ARGs on chromosomes,
## based on distribution of duplicate genes across environments.

## Fifth column: the expected number of duplicate ARGs on chromosomes in each category.
calc.expected.AR.chromosome.duplicates <- function(TableS3) {
    total.chromosome.duplicate.genes <- sum(TableS3$chromosomal_duplicate_genes)
    total.AR.chromosome.duplicate.genes <- sum(TableS3$chromosomal_AR_duplicate_genes)
    Fig3A.df <- TableS3 %>%
        mutate(expected_chromosomal_AR_duplicates = total.AR.chromosome.duplicate.genes * chromosomal_duplicate_genes/total.chromosome.duplicate.genes)
    return(Fig3A.df)
}

Fig3A.data <- calc.expected.AR.chromosome.duplicates(TableS3)

## enrichment of AR duplications on chromosomes.
Fig3A <- ggplot(Fig3A.data, aes(x=log2(chromosomal_duplicate_genes),
                              y=log2(chromosomal_AR_duplicate_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(chromosomal_duplicate_genes),
                     y = log2(chromosomal_AR_duplicate_genes),
                     xend = log2(chromosomal_duplicate_genes),
                     yend = log2(expected_chromosomal_AR_duplicates)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_chromosomal_AR_duplicates)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](chromosomal~~duplicate~genes))) +
    ylab(expression(log[2](chromosomal~duplicate~ARGs)))    

## Panel 3C.

## add a column for expected number of duplicate ARGs on plasmids, based on distribution
## of duplicate genes on plasmids across environments.

## Fifth column: the expected number of duplicate ARGs on plasmids in each category.
calc.expected.AR.plasmid.duplicates <- function(TableS3) {
    total.plasmid.duplicate.genes <- sum(TableS3$plasmid_duplicate_genes)
    total.AR.plasmid.duplicate.genes <- sum(TableS3$plasmid_AR_duplicate_genes)
    Fig3C.df <- TableS3 %>%
        mutate(expected_plasmid_AR_duplicates = total.AR.plasmid.duplicate.genes * plasmid_duplicate_genes/total.plasmid.duplicate.genes)
    return(Fig3C.df)
}

Fig3C.data <- calc.expected.AR.plasmid.duplicates(TableS3)

## enrichment of AR duplications on chromosomes.
Fig3C <- ggplot(Fig3C.data, aes(x=log2(plasmid_duplicate_genes),
                              y=log2(plasmid_AR_duplicate_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(plasmid_duplicate_genes),
                     y = log2(plasmid_AR_duplicate_genes),
                     xend = log2(plasmid_duplicate_genes),
                     yend = log2(expected_plasmid_AR_duplicates)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_plasmid_AR_duplicates)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](plasmid~~duplicate~genes))) +
    ylab(expression(log[2](plasmid~duplicate~ARGs)))

##########################
## old Figure 3. Visualization of the data in Table S3.
## Keep this visualization to help with debugging.
## TODO: consider keeping this figure but putting it into the supplement.

oldFig3A <- ggplot(Fig3.data, aes(x=Manual_Annotation,y=chromosomal_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,90000) +
    ggtitle("Chromosomal duplicate genes")

oldFig3B <- ggplot(Fig3.data, aes(x=Manual_Annotation,y=plasmid_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,90000) +
    ggtitle("Plasmid duplicate genes")

oldFig3C <- ggplot(Fig3.data, aes(x=Manual_Annotation,y=chromosomal_AR_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,500) +
    ggtitle("Chromosomal duplicate ARGs")

oldFig3D <- ggplot(Fig3.data, aes(x=Manual_Annotation,y=plasmid_AR_duplicate_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,500) +
    ggtitle("Plasmid duplicate ARGs")

oldFig3 <- plot_grid(oldFig3A, oldFig3B, oldFig3C, oldFig3D,labels=c('A','B','C','D'),ncol=2)
ggsave("../results/AR-gene-duplication/Fig3-old.pdf",oldFig3, height=4, width=12)
################################################################################
## Positive control 3: look at distribution of singleton AR genes on
## chromosomes and plasmids, to compare with Tables 3 and 4.

## Column 1
singleton.chromosome.genes.count <- singleton.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    group_by(Manual_Annotation) %>%
    summarize(chromosomal_singleton_genes = sum(chromosome_count))
gc() ## free memory when dealing with singleton.proteins.

## Column 2
singleton.plasmid.genes.count <- singleton.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    group_by(Manual_Annotation) %>%
    summarize(plasmid_singleton_genes = sum(plasmid_count))
gc() ## free memory when dealing with singleton.proteins.

## Column 3
singleton.AR.chromosome.genes.count <- singleton.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(chromosomal_AR_singleton_genes = sum(chromosome_count))
gc() ## free memory when dealing with singleton.proteins.

## Column 4
singleton.AR.plasmid.genes.count <- singleton.proteins %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated") %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(plasmid_AR_singleton_genes = sum(plasmid_count))
gc() ## free memory when dealing with singleton.proteins.

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
## S6 Figure: visualization of result in ControlTable3

S6Fig.data <- ControlTable3 %>%
    ## remove Unannotated isolates.
    filter(Manual_Annotation != "Unannotated")


S6FigA <- ggplot(S6Fig.data, aes(x=Manual_Annotation,y=chromosomal_singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,7000000) +
    ggtitle("Chromosomal singleton genes")

S6FigB <- ggplot(S6Fig.data, aes(x=Manual_Annotation,y=plasmid_singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation") + ylim(0,7000000) +
    ggtitle("Plasmid singleton genes")

S6FigC <- ggplot(S6Fig.data, aes(x=Manual_Annotation,y=chromosomal_AR_singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,9000) +
    ggtitle("Chromosomal singleton ARGs")

S6FigD <- ggplot(S6Fig.data, aes(x=Manual_Annotation,y=plasmid_AR_singleton_genes)) +
    geom_bar(stat="identity") +
    theme_classic() + coord_flip() + ylab("Count") +
    xlab("Isolate annotation")  + ylim(0,9000) +
    ggtitle("Plasmid singleton ARGs")

S6Fig <- plot_grid(S6FigA, S6FigB, S6FigC, S6FigD,labels=c('A','B','C','D'),ncol=2)
ggsave("../results/AR-gene-duplication/S6Fig.pdf",S6Fig)
################################################################################
## Finishing Figure 3.

## TODO: REFACTOR THIS CODE so that its organization makes sense.

## Use S6Fig.data to make panels 3B and 3D to finish making Figure 3. 

## Panel 3B.
## add a column for expected number of singleton ARGs on chromosomes, based on distribution
## of singleton genes across environments.

## Fifth column: the expected number of singleton ARGs on chromosomes in each category.
calc.expected.AR.chromosome.singletons <- function(S6Fig.data) {
    total.chromosome.singleton.genes <- sum(S6Fig.data$chromosomal_singleton_genes)
    total.AR.chromosome.singleton.genes <- sum(S6Fig.data$chromosomal_AR_singleton_genes)
    Fig3B.df <- S6Fig.data %>%
        mutate(expected_chromosomal_AR_singletons = total.AR.chromosome.singleton.genes * chromosomal_singleton_genes/total.chromosome.singleton.genes)
    return(Fig3B.df)
}

Fig3B.data <- calc.expected.AR.chromosome.singletons(S6Fig.data)

## enrichment of AR singletons on chromosomes.
Fig3B <- ggplot(Fig3B.data, aes(x=log2(chromosomal_singleton_genes),
                              y=log2(chromosomal_AR_singleton_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(chromosomal_singleton_genes),
                     y = log2(chromosomal_AR_singleton_genes),
                     xend = log2(chromosomal_singleton_genes),
                     yend = log2(expected_chromosomal_AR_singletons)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_chromosomal_AR_singletons)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](chromosomal~~singleton~genes))) +
    ylab(expression(log[2](chromosomal~singleton~ARGs)))    

## Panel 3D.

## add a column for expected number of singleton ARGs on plasmids, based on distribution
## of singleton genes on plasmids across environments.

## Fifth column: the expected number of singleton ARGs on plasmids in each category.
calc.expected.AR.plasmid.singletons <- function(S6Fig.data) {
    total.plasmid.singleton.genes <- sum(S6Fig.data$plasmid_singleton_genes)
    total.AR.plasmid.singleton.genes <- sum(S6Fig.data$plasmid_AR_singleton_genes)
    Fig3D.df <- S6Fig.data %>%
        mutate(expected_plasmid_AR_singletons = total.AR.plasmid.singleton.genes * plasmid_singleton_genes/total.plasmid.singleton.genes)
    return(Fig3D.df)
}

Fig3D.data <- calc.expected.AR.plasmid.singletons(S6Fig.data)

## enrichment of AR singletons on plasmids.
Fig3D <- ggplot(Fig3D.data, aes(x=log2(plasmid_singleton_genes),
                              y=log2(plasmid_AR_singleton_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(plasmid_singleton_genes),
                     y = log2(plasmid_AR_singleton_genes),
                     xend = log2(plasmid_singleton_genes),
                     yend = log2(expected_plasmid_AR_singletons)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_plasmid_AR_singletons)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](plasmid~singleton~genes))) +
    ylab(expression(log[2](plasmid~singleton~ARGs)))

### NOW *FINALLY* MAKE FIGURE 3.
Fig3 <- plot_grid(Fig3A,Fig3B,Fig3C,Fig3D,labels=c('A','B','C','D'),nrow=2)
ggsave("../results/AR-gene-duplication/Fig3.pdf",Fig3)

#######################################################################
## Figure 4.

## Where Figure 3 does vertical comparisons between panels in old Figure 3
## and Figure S6, Figure 4 does horizontal comparisons.

## This compares the distribution of genes on plasmids per category to the
## distribution of genes on chromosomes. This answers, for instance,
## where any categories are enriched with plasmids.

## Fig. 4A.
## add a column for expected number of duplicate genes on plasmids, based on distribution
## of duplicate genes on chromosomes across environments.

## Fifth column: the expected number of duplicate genes on plasmids in each category.
calc.expected.duplicate.genes.on.plasmids <- function(TableS3) {
    total.chromosome.duplicate.genes <- sum(TableS3$chromosomal_duplicate_genes)
    total.plasmid.duplicate.genes <- sum(TableS3$plasmid_duplicate_genes)
    Fig4A.df <- TableS3 %>%
        ## use exp(log(x)) trick to avoid numerical overflow.
        mutate(expected_duplicates_on_plasmids = exp(log(total.plasmid.duplicate.genes) + log(chromosomal_duplicate_genes) - log(total.chromosome.duplicate.genes)))
    return(Fig4A.df)
}

Fig4A.data <- calc.expected.duplicate.genes.on.plasmids(TableS3)

## enrichment of duplicates on plasmids compared to duplicates on chromosomes.
Fig4A <- ggplot(Fig4A.data, aes(x=log2(chromosomal_duplicate_genes),
                              y=log2(plasmid_duplicate_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(chromosomal_duplicate_genes),
                     y = log2(plasmid_duplicate_genes),
                     xend = log2(chromosomal_duplicate_genes),
                     yend = log2(expected_duplicates_on_plasmids)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_duplicates_on_plasmids)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](chromosomal~~duplicate~genes))) +
    ylab(expression(log[2](plasmid~duplicate~genes)))    


## Fig. 4B.
## add a column for expected number of singleton genes on plasmids, based on distribution
## of singleton genes on chromosomes across environments.

## Fifth column: the expected number of singleton genes on plasmids in each category.
calc.expected.singleton.genes.on.plasmids <- function(S6Fig.data) {
    total.chromosome.singleton.genes <- sum(S6Fig.data$chromosomal_singleton_genes)
    total.plasmid.singleton.genes <- sum(S6Fig.data$plasmid_singleton_genes)
    Fig4B.df <- S6Fig.data %>%
        ## use exp(log(x)) trick to avoid numerical overflow.
        mutate(expected_singletons_on_plasmids = exp(log(total.plasmid.singleton.genes) + log(chromosomal_singleton_genes) - log(total.chromosome.singleton.genes)))
    return(Fig4B.df)
}

Fig4B.data <- calc.expected.singleton.genes.on.plasmids(S6Fig.data)

## enrichment of singletons on plasmids compared to singletons on chromosomes.
Fig4B <- ggplot(Fig4B.data, aes(x=log2(chromosomal_singleton_genes),
                              y=log2(plasmid_singleton_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(chromosomal_singleton_genes),
                     y = log2(plasmid_singleton_genes),
                     xend = log2(chromosomal_singleton_genes),
                     yend = log2(expected_singletons_on_plasmids)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_singletons_on_plasmids)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](chromosomal~~singleton~genes))) +
    ylab(expression(log[2](plasmid~singleton~genes)))    

## Panel 4C.
## add a column for expected number of duplicate ARGs on plasmids, based on distribution
## of duplicate ARGs on chromosomes across environments.

## Fifth column: the expected number of duplicate ARGs on plasmids, based on
## distribution of duplicate ARGs on chromosomes in each category.
calc.expected.AR.plasmid.duplicates.compared.to.chromosomes <- function(TableS3) {
    total.AR.plasmid.duplicate.genes <- sum(TableS3$plasmid_AR_duplicate_genes)
    total.AR.chromosome.duplicate.genes <- sum(TableS3$chromosomal_AR_duplicate_genes)
    Fig4C.df <- TableS3 %>%
        mutate(expected_plasmid_AR_duplicates = total.AR.plasmid.duplicate.genes * chromosomal_AR_duplicate_genes/total.AR.chromosome.duplicate.genes)
    return(Fig4C.df)
}

Fig4C.data <- calc.expected.AR.plasmid.duplicates.compared.to.chromosomes(TableS3)

## enrichment of AR duplications on plasmids compared to chromosomes.
Fig4C <- ggplot(Fig4C.data, aes(x=log2(chromosomal_AR_duplicate_genes),
                              y=log2(plasmid_AR_duplicate_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(chromosomal_AR_duplicate_genes),
                     y = log2(plasmid_AR_duplicate_genes),
                     xend = log2(chromosomal_AR_duplicate_genes),
                     yend = log2(expected_plasmid_AR_duplicates)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_plasmid_AR_duplicates)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](chromosomal~~duplicate~ARGs))) +
    ylab(expression(log[2](plasmid~duplicate~ARGs)))    


## Panel 4D.
## add a column for expected number of singleton ARGs on plasmids, based on distribution
## of singleton ARGs on chromosomes across environments.

## Fifth column: the expected number of singleton ARGs on plasmids, based on
## distribution of singleton ARGs on chromosomes in each category.
calc.expected.AR.plasmid.singletons.compared.to.chromosomes <- function(S6Fig.data) {
    total.AR.plasmid.singleton.genes <- sum(S6Fig.data$plasmid_AR_singleton_genes)
    total.AR.chromosome.singleton.genes <- sum(S6Fig.data$chromosomal_AR_singleton_genes)
    Fig4D.df <- S6Fig.data %>%
        mutate(expected_plasmid_AR_singletons = total.AR.plasmid.singleton.genes * chromosomal_AR_singleton_genes/total.AR.chromosome.singleton.genes)
    return(Fig4D.df)
}

Fig4D.data <- calc.expected.AR.plasmid.singletons.compared.to.chromosomes(S6Fig.data)

## enrichment of AR singletons on plasmids compared to chromosomes.
Fig4D <- ggplot(Fig4D.data, aes(x=log2(chromosomal_AR_singleton_genes),
                              y=log2(plasmid_AR_singleton_genes),
                              label=Manual_Annotation)) +
    theme_classic() +
    geom_segment(aes(x = log2(chromosomal_AR_singleton_genes),
                     y = log2(plasmid_AR_singleton_genes),
                     xend = log2(chromosomal_AR_singleton_genes),
                     yend = log2(expected_plasmid_AR_singletons)),
                 color="light gray",size=0.2,linetype='dashed') +
    geom_point(color='red') +
    geom_line(aes(y=log2(expected_plasmid_AR_singletons)),color='yellow') +
    geom_text_repel(size=2.5) +
    xlab(expression(log[2](chromosomal~~singleton~ARGs))) +
    ylab(expression(log[2](plasmid~singleton~ARGs)))    


Fig4 <- plot_grid(Fig4A,Fig4B,Fig4C,Fig4D,labels=c('A','B','C','D'),nrow=2)
ggsave("../results/AR-gene-duplication/Fig4.pdf",Fig4)


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
    table.df.to.figure.df()

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

################################################################################

## do a simple analysis of the frequency of different annotations overall,
## and in different ecological annotation categories.

total.annotation.freq.table <- duplicate.proteins %>%
    group_by(product) %>%
    summarize(annotation.count = n()) %>%
    filter(annotation.count > 1) %>%
    arrange(desc(annotation.count))

## remove MGE associations.
no.MGE.annotation.freq.table <- duplicate.proteins %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(product) %>%
    summarize(annotation.count = n()) %>%
    filter(annotation.count > 1) %>%
    arrange(desc(annotation.count))

## do similar analyses based on sequence identity.

total.seq.freq.table <- duplicate.proteins %>%
    group_by(sequence, product) %>%
    summarize(seq.count = n()) %>%
    filter(seq.count > 1) %>%
    arrange(desc(seq.count))

## remove MGE associations.
no.MGE.seq.freq.table <- duplicate.proteins %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(sequence, product) %>%
    summarize(seq.count = n()) %>%
    filter(seq.count > 1) %>%
    arrange(desc(seq.count))

## let's look at the annotations and sequences of high frequency duplicated
## proteins in each environment, after removing MGEs.

## this function filters duplicate proteins by manual annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
make.annotation.freq.table <- function(manual.annot.string) {
    duplicate.proteins %>%
    filter(Manual_Annotation == manual.annot.string) %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(product) %>%
    summarize(annotation.count = n()) %>%
    filter(annotation.count > 1) %>%
    arrange(desc(annotation.count))
}

animal.host.annotation.freq.table <- make.annotation.freq.table("Animal-host")
anthropogenic.annotation.freq.table <- make.annotation.freq.table("Anthropogenic-environment")
human.host.annotation.freq.table <- make.annotation.freq.table("Human-host")
food.annotation.freq.table <- make.annotation.freq.table("Food")
agriculture.annotation.freq.table <- make.annotation.freq.table("Agriculture")
unannotated.annotation.freq.table <- make.annotation.freq.table("Unannotated")
marine.annotation.freq.table <- make.annotation.freq.table("Marine")
sediment.annotation.freq.table <- make.annotation.freq.table("Sediment")
livestock.annotation.freq.table <- make.annotation.freq.table("Livestock")
soil.annotation.freq.table <- make.annotation.freq.table("Soil")
freshwater.annotation.freq.table <- make.annotation.freq.table("Freshwater")
terrestrial.annotation.freq.table <- make.annotation.freq.table("Terrestrial")
plant.host.annotation.freq.table <- make.annotation.freq.table("Plant-host")
fungal.host.annotation.freq.table <- make.annotation.freq.table("Fungal-host")

## this function filters duplicate proteins by manual annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
make.seq.freq.table <- function(manual.annot.string) {
    duplicate.proteins %>%
    filter(Manual_Annotation == manual.annot.string) %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(sequence,product) %>%
    summarize(seq.count = n()) %>%
    filter(seq.count > 1) %>%
    arrange(desc(seq.count))
}

animal.host.seq.freq.table <- make.seq.freq.table("Animal-host")
anthropogenic.seq.freq.table <- make.seq.freq.table("Anthropogenic-environment")
human.host.seq.freq.table <- make.seq.freq.table("Human-host")
food.seq.freq.table <- make.seq.freq.table("Food")
agriculture.seq.freq.table <- make.seq.freq.table("Agriculture")
unannotated.seq.freq.table <- make.seq.freq.table("Unannotated")
marine.seq.freq.table <- make.seq.freq.table("Marine")
sediment.seq.freq.table <- make.seq.freq.table("Sediment")
livestock.seq.freq.table <- make.seq.freq.table("Livestock")
soil.seq.freq.table <- make.seq.freq.table("Soil")
freshwater.seq.freq.table <- make.seq.freq.table("Freshwater")
terrestrial.seq.freq.table <- make.seq.freq.table("Terrestrial")
plant.host.seq.freq.table <- make.seq.freq.table("Plant-host")
fungal.host.seq.freq.table <- make.seq.freq.table("Fungal-host")

##########################################

## what happens if we restrict to annotations and sequences which are sometimes
## found on plasmids?

on.plas.annotation.freq.table <- duplicate.proteins %>%
    filter(plasmid_count >= 1) %>%
    group_by(product) %>%
    summarize(annotation.count = n()) %>%
    filter(annotation.count > 1) %>%
    arrange(desc(annotation.count))

## remove MGE associations.
on.plas.no.MGE.annotation.freq.table <- duplicate.proteins %>%
    filter(plasmid_count >= 1) %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(product) %>%
    summarize(annotation.count = n()) %>%
    filter(annotation.count > 1) %>%
    arrange(desc(annotation.count))

## do similar analyses based on sequence identity.

on.plas.seq.freq.table <- duplicate.proteins %>%
    filter(plasmid_count >= 1) %>%
    group_by(sequence, product) %>%
    summarize(seq.count = n()) %>%
    filter(seq.count > 1) %>%
    arrange(desc(seq.count))

## remove MGE associations.
on.plas.no.MGE.seq.freq.table <- duplicate.proteins %>%
    filter(plasmid_count >= 1) %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(sequence, product) %>%
    summarize(seq.count = n()) %>%
    filter(seq.count > 1) %>%
    arrange(desc(seq.count))

####################
make.plas.annotation.freq.table <- function(manual.annot.string) {
    duplicate.proteins %>%
        filter(plasmid_count >= 1) %>%
        filter(Manual_Annotation == manual.annot.string) %>%
        filter(!str_detect(.$product,IS.keywords)) %>%
        group_by(product) %>%
        summarize(annotation.count = n()) %>%
        filter(annotation.count > 1) %>%
        arrange(desc(annotation.count))
}

animal.host.plas.annotation.freq.table <- make.plas.annotation.freq.table("Animal-host")
anthropogenic.plas.annotation.freq.table <- make.plas.annotation.freq.table("Anthropogenic-environment")
human.host.plas.annotation.freq.table <- make.plas.annotation.freq.table("Human-host")
food.plas.annotation.freq.table <- make.plas.annotation.freq.table("Food")
agriculture.plas.annotation.freq.table <- make.plas.annotation.freq.table("Agriculture")
unannotated.plas.annotation.freq.table <- make.plas.annotation.freq.table("Unannotated")
marine.plas.annotation.freq.table <- make.plas.annotation.freq.table("Marine")
sediment.plas.annotation.freq.table <- make.plas.annotation.freq.table("Sediment")
livestock.plas.annotation.freq.table <- make.plas.annotation.freq.table("Livestock")
soil.plas.annotation.freq.table <- make.plas.annotation.freq.table("Soil")
freshwater.plas.annotation.freq.table <- make.plas.annotation.freq.table("Freshwater")
terrestrial.plas.annotation.freq.table <- make.plas.annotation.freq.table("Terrestrial")
plant.host.plas.annotation.freq.table <- make.plas.annotation.freq.table("Plant-host")
fungal.host.plas.annotation.freq.table <- make.plas.annotation.freq.table("Fungal-host")

## this function filters duplicate proteins by manual annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
make.on.plas.seq.freq.table <- function(manual.annot.string) {
    duplicate.proteins %>%
        filter(plasmid_count >= 1) %>%
        filter(Manual_Annotation == manual.annot.string) %>%
        filter(!str_detect(.$product,IS.keywords)) %>%
        group_by(sequence,product) %>%
        summarize(seq.count = n()) %>%
        filter(seq.count > 1) %>%
        arrange(desc(seq.count))
}

animal.host.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Animal-host")
anthropogenic.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Anthropogenic-environment")
human.host.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Human-host")
food.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Food")
agriculture.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Agriculture")
unannotated.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Unannotated")
marine.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Marine")
sediment.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Sediment")
livestock.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Livestock")
soil.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Soil")
freshwater.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Freshwater")
terrestrial.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Terrestrial")
plant.host.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Plant-host")
fungal.host.on.plas.seq.freq.table <- make.on.plas.seq.freq.table("Fungal-host")



##########################################
## NOTES AND IDEAS
## signal decomposition algorithms: non-negative matrix factorization,
## ICA, etc.
## represent strains by duplicated genes. Then factorize
## those strains into non-negative combinations of sequences/annotations
## and their counts.
## for a classifier, weight sequences/annotations that are most informative
## of particular environments.

## if I filter out MGEs, then I probably get more insight into WHAT functions
## are being selected.
## BUT MGEs probably carry significant information about niches!
## assume that gene flow networks are more tightly connected within a niche
## in comparison to between niches (since higher probability of interaction).
## This idea goes back at least to Smillie et al. 2011 in Nature from Eric Alm's group.
