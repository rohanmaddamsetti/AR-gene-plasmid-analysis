## Aim1-analysis.R by Rohan Maddamsetti.

## analyse the distribution of AR genes on chromosomes versus plasmids in
## fully-sequenced genomes and plasmids in the NCBI Nucleotide database.

## CRITICAL ANALYSIS TODO: Make sure numbers in genome.database,
## gbk.annotation, and all.proteins, duplicate.proteins, and
## singleton.proteins, in terms of number of isolates in each
## category, are COMPLETELY consistent with each other.

## CRITICAL ANALYSIS TODO FOR A FOLLOW UP PAPER:
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
    ## get species name annotation from genome.database
    left_join(genome.database)

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
## PSNIH2 is missing, but the others are in genome.database
## and AR.results. PSNIH2 is not in the original
## prokaryotes.txt genome reports file.
################################################################################

## Simple analysis of recent protein duplications and HGT between
## chromosome and plasmids.

## remove all genes with the following keywords in the "product" annotation
IS.keywords <- "IS|transposon|Transposase|transposase|hypothetical protein|Phage|phage|integrase|Integrase|tail|intron|Mobile|mobile|antitoxin|toxin|capsid|plasmid|Plasmid|conjug"

## remove all genes that match Elongation Factor Tu (2 copies in most bacteria).
EFTu.keywords <- "Tu | Tu|-Tu"

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

##########################################
## CRITICAL BUG TO FIX:
## there are 7,910 isolates that have annotated proteins in their Genbank
## annotation. 634 are missing protein (CDS) annotation in their genome.
## for now, I have restricted these data by filtering on Annotation_Accession--
###########################################################################
## Analysis for Figure 1C.

## Statistical analysis for isolates with duplicated ARGs
## goes into Supplementary Table S1.

## return the first column for several tables.
## shows the number of isolates in each category.
make.isolate.totals.col <- function(gbk.annotation) {
    isolate.totals <- gbk.annotation %>%
        ## remove Unannotated isolates as we're working with gbk.annotation.
        filter(Annotation != "Unannotated") %>%
        filter(Annotation != "blank") %>%
        group_by(Annotation) %>%
        summarize(total_isolates = n()) %>%
        arrange(desc(total_isolates))
    return(isolate.totals)
}

make.TableS1 <- function(gbk.annotation, duplicate.genes) {

    ## count the number of isolates with duplicated ARGs in each category.
    AR.category.counts <- duplicate.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_duplicated_ARGs = n()) %>%
        arrange(desc(isolates_with_duplicated_ARGs))
    
    ## join columns to make Table 1 with raw data.
    raw.TableS1 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(AR.category.counts) %>%
        mutate(isolates_with_duplicated_ARGs =
                   replace_na(isolates_with_duplicated_ARGs,0)) %>%
        arrange(desc(isolates_with_duplicated_ARGs))
    
    ## For consistency with Figure 1, use the distribution of sampled isolates as the null
    ## distribution.
    calc.expected.isolates.with.ARGs <- function(raw.TableS1) {
        sum.total.isolates <- sum(raw.TableS1$total_isolates)
        total.isolates.with.duplicated.ARGs <- sum(raw.TableS1$isolates_with_duplicated_ARGs)
        Table <- raw.TableS1 %>%
            mutate(expected_isolates_with_duplicated_ARGs = total.isolates.with.duplicated.ARGs * total_isolates/sum.total.isolates)
        return(Table)
    }
    
    calc.isolate.AR.gene.enrichment.pvals <- function(raw.TableS1) {
        sum.total.isolates <- sum(raw.TableS1$total_isolates)
        total.isolates.with.duplicated.ARGs <- sum(raw.TableS1$isolates_with_duplicated_ARGs)
        Table <- raw.TableS1 %>%
            rowwise() %>%
            mutate(binom.test.pval = binom.test
            (
                x = isolates_with_duplicated_ARGs,
                n = total.isolates.with.duplicated.ARGs,
                p = total_isolates/sum.total.isolates
            )$p.value) %>%
            mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
            ## use Benjamini-Hochberg p-value correction.
            select(-binom.test.pval) ## drop original p-value after the correction.
        return(Table)
    }

    ## Add a third column: expected number of isolates with duplicated ARGs,
    ## based on the percentage of isolates with duplicated genes.
    TableS1 <- raw.TableS1 %>% calc.expected.isolates.with.ARGs() %>%
        ## Add a fourth column: p-values for deviation from
        ## expected number of duplicated ARGs, using binomial test,
        ## correcting for multiple tests.
        calc.isolate.AR.gene.enrichment.pvals()

    return(TableS1)
}

TableS1 <- make.TableS1(gbk.annotation, duplicate.genes)
## write Supplementary Table S1 to file.
write.csv(x=TableS1, file="../results/TableS1.csv")

######################
## Control: does the number of isolates with duplicate genes
## follow the sampling distribution of isolates?

## Most follow the expected distribution.
## however, isolates from animal-hosts are signficantly depleted
## in duplicate genes: FDR-corrected p < 10^-5
## while isolates from livestock and agriculture are both significantly enriched
## in duplicate genes: FDR-corrected p < 0.004 for both.

run.duplicate.gene.control <- function(gbk.annotation, duplicate.proteins) {
    
    ## count the number of isolates with duplications in each category.
    isolates.with.duplicate.genes <- duplicate.proteins %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_duplicate_genes = n()) %>%
        arrange(desc(isolates_with_duplicate_genes))
    
    ControlData <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(isolates.with.duplicate.genes) 
    
    ## For consistency with Figure 1, use the distribution of sampled isolates as the null
    ## distribution.

    sum.total.isolates <- sum(ControlData$total_isolates)
    total.isolates.with.duplicate.genes <- sum(ControlData$isolates_with_duplicate_genes)

    Table <- ControlData %>%
        mutate(expected_isolates_with_duplicate_genes = total.isolates.with.duplicate.genes * total_isolates/sum.total.isolates) %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test
        (
            x = isolates_with_duplicate_genes,
            n = total.isolates.with.duplicate.genes,
            p = total_isolates/sum.total.isolates
        )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        ## use Benjamini-Hochberg p-value correction.
        select(-binom.test.pval) ## drop original p-value after the correction.
    return(Table)
}

duplicate.gene.control.df <- run.duplicate.gene.control(gbk.annotation, duplicate.proteins)

######################
## Control: does the distribution of ARG singletons
## (i.e. genes that have NOT duplicated) follow the distribution
## of sampled isolates?

## No categories are enriched with singleton AR genes,
## as most isolates have a gene that matches an antibiotic keyword.
## Animal-host isolates are depleted (perhaps due to aphid bacteria isolates?)

run.singleton.ARG.control <- function(gbk.annotation, singleton.proteins) {

## count the number of isolates with singleton AR genes in each category.
    AR.singleton.category.counts <- singleton.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_singleton_ARGs = n()) %>%
        arrange(desc(isolates_with_singleton_ARGs))
    gc() ## free memory.

    calc.expected.isolates.with.singleton.ARGs <- function(raw.Table) {
        summed.isolates <- sum(raw.Table$total_isolates)
        total.isolates.with.singleton.ARGs <- sum(raw.Table$isolates_with_singleton_ARGs)
        Table <- raw.Table %>%
            mutate(expected_isolates_with_singleton_ARGs = total.isolates.with.singleton.ARGs * total_isolates/summed.isolates)
        return(Table)
    }

    calc.isolate.singleton.ARG.enrichment.pvals <- function(raw.Table) {
    
        summed.isolates <- sum(raw.Table$total_isolates)
        total.isolates.with.singleton.ARGs <- sum(raw.Table$isolates_with_singleton_ARGs)

        Table <- raw.Table %>%
            rowwise() %>%
            mutate(binom.test.pval = binom.test(
                       x = isolates_with_singleton_ARGs,
                       n = total.isolates.with.singleton.ARGs,
                       p = total_isolates/summed.isolates
                   )$p.value) %>%
            mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
            ## use Benjamini-Hochberg p-value correction.
            select(-binom.test.pval) ## drop original p-value after the correction.
        return(Table)
    }
    
    Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(AR.singleton.category.counts) %>%
        mutate(isolates_with_singleton_ARGs=replace_na(isolates_with_singleton_ARGs,0)) %>%
        arrange(desc(isolates_with_singleton_ARGs)) %>%
        calc.expected.isolates.with.singleton.ARGs() %>%
        calc.isolate.singleton.ARG.enrichment.pvals()
    
    return(Table)
}

## This data frame will be used for Figure 1C.
ControlTable1 <- run.singleton.ARG.control(gbk.annotation, singleton.proteins)
## write Control 1 to file.
write.csv(x=ControlTable1, file="../results/ControlTable1.csv")

gc() ## free memory after dealing with singleton data.
###################################
###################################
## Tables to plot percentage of genes on plasmids, for Figure 1C.

## Table 2. Show number of duplicated genes on chromosomes, and number of
## duplicate genes on plasmids, for each category, for duplicated genes
## and duplicated AR genes. This table of raw data goes into the text. Then,
## sum over all categories for a 2x2 contingency table and report the result of a
## Fisher's exact test for asssociation between duplicated AR genes and plasmids.

make.Table2 <- function(duplicate.proteins) {
    ## Column 1
    duplicate.chromosome.genes.count <- duplicate.proteins %>%
        group_by(Annotation) %>%
        summarize(chromosomal_duplicate_genes = sum(chromosome_count))
    
    ## Column 2
    duplicate.plasmid.genes.count <- duplicate.proteins %>%
        group_by(Annotation) %>%
        summarize(plasmid_duplicate_genes = sum(plasmid_count))
    
    ## Column 3
    duplicate.chromosome.ARGs.count <- duplicate.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        group_by(Annotation) %>%
        summarize(chromosomal_duplicate_ARGs = sum(chromosome_count))
    
    ## Column 4
    duplicate.plasmid.ARGs.count <- duplicate.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        group_by(Annotation) %>%
        summarize(plasmid_duplicate_ARGs = sum(plasmid_count))
    
    Table <- duplicate.chromosome.genes.count %>%
        left_join(duplicate.plasmid.genes.count) %>%
        left_join(duplicate.chromosome.ARGs.count) %>%
        mutate(chromosomal_duplicate_ARGs =
                   replace_na(chromosomal_duplicate_ARGs, 0)) %>%
        left_join(duplicate.plasmid.ARGs.count) %>%
        mutate(plasmid_duplicate_ARGs = replace_na(plasmid_duplicate_ARGs, 0)) %>%
        arrange(desc(plasmid_duplicate_ARGs))
    
    return(Table)
}

Table2 <- make.Table2(duplicate.proteins)
## write Table 2 to file.
write.csv(x=Table2,file="../results/Table2.csv")

################
## Analysis of Table 2: Duplicate ARGs are associated with plasmids.

plasmid.chromosome.duplicate.ARG.contingency.test <- function(Table2) {
    ## get values for Fisher's exact test.
    total.chr.AR.duplicates <- sum(Table2$chromosomal_duplicate_ARGs)
    total.plasmid.AR.duplicates <- sum(Table2$plasmid_duplicate_ARGs)

    total.chr.duplicates <- sum(Table2$chromosomal_duplicate_genes)
    total.plasmid.duplicates <- sum(Table2$plasmid_duplicate_genes)

    total.nonAR.chr.duplicates <- total.chr.duplicates - total.chr.AR.duplicates
    total.nonAR.plasmid.duplicates <- total.plasmid.duplicates - total.plasmid.AR.duplicates

    contingency.table <- matrix(c(total.chr.AR.duplicates,
                                  total.plasmid.AR.duplicates,
                                  total.nonAR.chr.duplicates,
                                  total.nonAR.plasmid.duplicates),nrow=2)
    ## label the rows and columns of the contingency table.
    rownames(contingency.table) <- c("chromosome","plasmid")
    colnames(contingency.table) <- c("AR duplicate genes","non-AR duplicate genes")

    ## p < 1e-300
    print(fisher.test(contingency.table))
    print(fisher.test(contingency.table)$p.value)
    return(contingency.table)
}

plasmid.chromosome.duplicate.ARG.contingency.test(Table2)

####################
## Control Table 2: look at distribution of singleton AR genes on
## chromosomes and plasmids.

## This control shows singleton AR genes are highly enriched on plasmids,
## based on a comparison with the distribution of singleton genes overall.
## Therefore AR genes are generally associated with plasmids, regardless of
## status of being a duplication or not.

## This does NOT invalidate the main result of this analysis, that duplicate AR
## genes are more enriched on plasmids in comparison to the distribution of
## duplicate genes overall.

## Notably, the majority of duplicate ARGs are on plasmids, while the
## majority of singleton ARGs are on chromosomes.

make.ControlTable2 <- function(singleton.proteins) {
    ## Column 1
    singleton.chromosome.genes.count <- singleton.proteins %>%
        group_by(Annotation) %>%
        summarize(chromosomal_singleton_genes = sum(chromosome_count))
    gc() ## free memory when dealing with singleton.proteins.

    ## Column 2
    singleton.plasmid.genes.count <- singleton.proteins %>%
        group_by(Annotation) %>%
        summarize(plasmid_singleton_genes = sum(plasmid_count))
    gc() ## free memory when dealing with singleton.proteins.

    ## Column 3
    singleton.chromosome.ARGs.count <- singleton.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        group_by(Annotation) %>%
        summarize(chromosomal_singleton_ARGs = sum(chromosome_count))
    gc() ## free memory when dealing with singleton.proteins.
    
    ## Column 4
    singleton.plasmid.ARGs.count <- singleton.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        group_by(Annotation) %>%
        summarize(plasmid_singleton_ARGs = sum(plasmid_count))
    gc() ## free memory when dealing with singleton.proteins.
    
    ControlTable2 <- singleton.chromosome.genes.count %>%
        left_join(singleton.plasmid.genes.count) %>%
        left_join(singleton.chromosome.ARGs.count) %>%
        mutate(chromosomal_singleton_ARGs =
                    replace_na(chromosomal_singleton_ARGs, 0)) %>%
        left_join(singleton.plasmid.ARGs.count) %>%
        mutate(plasmid_singleton_ARGs = replace_na(plasmid_singleton_ARGs, 0)) %>%
        arrange(desc(plasmid_singleton_ARGs))
}

ControlTable2 <- make.ControlTable2(singleton.proteins)
gc()

plasmid.chromosome.singleton.ARG.contingency.test <- function(ControlTable2) {
    ## get values for Fisher's exact test.
    total.chr.AR.singletons <- sum(ControlTable2$chromosomal_singleton_ARGs)
    total.plasmid.AR.singletons <- sum(ControlTable2$plasmid_singleton_ARGs)

    total.chr.singletons <- sum(ControlTable2$chromosomal_singleton_genes)
    total.plasmid.singletons <- sum(ControlTable2$plasmid_singleton_genes)

    total.nonAR.chr.singletons <- total.chr.singletons - total.chr.AR.singletons
    total.nonAR.plasmid.singletons <- total.plasmid.singletons - total.plasmid.AR.singletons

    contingency.table <- matrix(c(total.chr.AR.singletons,
                                  total.plasmid.AR.singletons,
                                  total.nonAR.chr.singletons,
                                  total.nonAR.plasmid.singletons),nrow=2)
    ## label the rows and columns of the contingency table.
    rownames(contingency.table) <- c("chromosome","plasmid")
    colnames(contingency.table) <- c("AR singleton genes","non-AR singleton genes")
    
    print(fisher.test(contingency.table))
    print(fisher.test(contingency.table)$p.value)

    return(contingency.table)
}

plasmid.chromosome.singleton.ARG.contingency.test(ControlTable2)

################################################################################
## Supplementary Figure S1:
## Histogram visualization of ARGs on plasmids and chromosomes.

make.S1Fig <- function(ControlTable2,Table2) {

    S1Fig.data <- full_join(ControlTable2, Table2) %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(c("Human-host","Livestock","Animal-host",
                                  "Anthropogenic-environment", "Food", "Freshwater",
                                  "Agriculture", "Sediment", "Soil", "Plant-host",
                                  "Marine","Terrestrial", "Fungal-host"))))

    format.S1Fig.panel <- function(S1.panel, xlim.max, my.title) {
        full.S1.panel <- S1.panel +
            geom_bar(stat="identity") +
            theme_classic() + xlab("Count") +
            ylab("Isolate annotation") +
            xlim(0, xlim.max) +
            ggtitle(my.title)
        return(full.S1.panel)
    }
    
    S1FigA <- ggplot(S1Fig.data, aes(y=Annotation, x=chromosomal_singleton_genes)) %>%
        format.S1Fig.panel(41000000, "Chromosomal singleton genes")
    
    S1FigB <- ggplot(S1Fig.data, aes(y=Annotation,x=plasmid_singleton_genes)) %>%
        format.S1Fig.panel(41000000, "Plasmid singleton genes")
    
    S1FigC <- ggplot(S1Fig.data, aes(y=Annotation,x=chromosomal_singleton_ARGs)) %>%
        format.S1Fig.panel(50000, "Chromosomal singleton ARGs")
        
    S1FigD <- ggplot(S1Fig.data, aes(y=Annotation,x=plasmid_singleton_ARGs)) %>%
        format.S1Fig.panel(50000, "Plasmid singleton ARGs")
    
    S1FigE <- ggplot(S1Fig.data, aes(y=Annotation,x=chromosomal_duplicate_genes)) %>%
        format.S1Fig.panel(550000, "Chromosomal duplicate genes")
    
    S1FigF <- ggplot(S1Fig.data, aes(y=Annotation,x=plasmid_duplicate_genes)) %>%
        format.S1Fig.panel(550000,"Plasmid duplicate genes")
    
    S1FigG <- ggplot(S1Fig.data, aes(y=Annotation,x=chromosomal_duplicate_ARGs)) %>%
        format.S1Fig.panel(4000, "Chromosomal duplicate ARGs")
    
    S1FigH <- ggplot(S1Fig.data, aes(y=Annotation,x=plasmid_duplicate_ARGs)) %>%
        format.S1Fig.panel(4000, "Plasmid duplicate ARGs")
    
    S1Fig <- plot_grid(S1FigA, S1FigB, S1FigC, S1FigD,
                       S1FigE, S1FigF, S1FigG, S1FigH,
                       labels=c('A','B','C','D','E','F','G','H'), ncol=2)
    return(S1Fig)
}

S1Fig <- make.S1Fig(ControlTable2,Table2)
ggsave("../results/S1Fig.pdf",S1Fig,width=9,height=11)

#######################
## Figure 2: stacked bar plot version of S1Fig.

## Make Figure 2A:
## Stacked bar chart of multi-copy proteins
## on chromosomes and plasmids.

categorize.as.MGE.ARG.or.other <- function(product) {
    if (is.na(product))
        return("Other function")
    else if (str_detect(product, antibiotic.keywords))
        return("ARG")
    else if (str_detect(product, EFTu.keywords))
        return("EF-Tu")
    else if (str_detect(product, IS.keywords))
        return("MGE")
    else
        return("Other function")
}

Fig2A.data <- duplicate.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation, Category) %>%
    summarize(Plasmid = sum(plasmid_count), Chromosome = sum(chromosome_count)) %>%
    pivot_longer(cols = c("Plasmid", "Chromosome"),
                 names_to = "Episome",
                 values_to = "Count") %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                              "Anthropogenic-environment", "Food", "Freshwater",
                              "Agriculture", "Sediment", "Soil", "Plant-host",
                              "Marine","Terrestrial", "Fungal-host"))))

Fig2B.data <- singleton.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation, Category) %>%
    summarize(Plasmid = sum(plasmid_count), Chromosome = sum(chromosome_count)) %>%
    pivot_longer(cols = c("Plasmid", "Chromosome"),
                 names_to = "Episome",
                 values_to = "Count") %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(c("Human-host","Livestock","Animal-host",
                              "Anthropogenic-environment", "Food", "Freshwater",
                              "Agriculture", "Sediment", "Soil", "Plant-host",
                              "Marine","Terrestrial", "Fungal-host"))))


Fig2A <- ggplot(Fig2A.data, aes(y = Annotation, x = Count, fill = Category)) +
    geom_bar(stat="identity") +
    facet_wrap(.~Episome) +
    theme_classic() +
    ggtitle("Distribution of multi-copy proteins") +
    scale_x_continuous(labels=fancy_scientific)

Fig2.legend <- cowplot::get_legend(Fig2A)

## now remove the legend from Fig. 2A.
Fig2A <- Fig2A + guides(fill = FALSE)

Fig2B <- ggplot(Fig2B.data, aes(y = Annotation, x = Count, fill = Category)) +
    geom_bar(stat="identity") +
    facet_wrap(.~Episome) +
    theme_classic() +
    ggtitle("Distribution of single-copy proteins") +
    guides(fill = FALSE) +
    scale_x_continuous(labels=fancy_scientific)

rm(Fig2A.data) ## to save memory.
rm(Fig2B.data) ## to save memory.
gc() ## run garbage collection.

Fig2 <- plot_grid(Fig2A, Fig2B, labels = c("A", "B"), ncol = 1)
ggsave("../results/Fig2.pdf", Fig2)

##################################################################################
## Figure 1 A & B: Diagram of the analysis workflow, made in Inkscape/Illustrator.
##################################################################################
## Figure 1C: main figure, showing enrichment of AR duplicates
## in human hosts and livestock.



make.Fig1C.df <- function(TableS1, ControlTable1, Table2, ControlTable2, duplicate.gene.control.df) {
    ## join duplicate and singleton tables to make Fig 1C.

    ## have to remove p-values from the two tables, because
    ## the column names are the same, but the values are different
    ## (because these are two different tests).
    no.pval.TableS1 <- select(TableS1, -corrected.pval)
    no.pval.ControlTable1 <- select(ControlTable1, -corrected.pval)
    
    Fig1C.df <- no.pval.TableS1 %>%
        full_join(no.pval.ControlTable1) %>%
        full_join(Table2) %>%
        full_join(ControlTable2) %>%
        full_join(duplicate.gene.control.df)
    
    total_isolates.sum <- sum(Fig1C.df$total_isolates)
    isolates_with_duplicate_genes.sum <- sum(Fig1C.df$isolates_with_duplicate_genes)
    isolates_with_singleton_ARGs.sum <- sum(Fig1C.df$isolates_with_singleton_ARGs)
    isolates_with_duplicated_ARGs.sum <- sum(Fig1C.df$isolates_with_duplicated_ARGs)
    
    Fig1C.df <- Fig1C.df %>%
    ## calculate y-coordinates for line for duplicate genes.
        mutate(yvals.for.isolates_with_duplicate_genes.line = total_isolates * isolates_with_duplicate_genes.sum/total_isolates.sum) %>%
        ## calculate y-coordinates for line for singleton ARGs.
        mutate(yvals.for.isolates_with_singleton_ARGs.line = total_isolates * isolates_with_singleton_ARGs.sum/total_isolates.sum) %>%
        ## calculate y-coordinates for line for duplicate ARGs.
        mutate(yvals.for.isolates_with_duplicated_ARGs.line = total_isolates * isolates_with_duplicated_ARGs.sum/total_isolates.sum) %>%
        ## calculate the percentage of genes on plasmids for symbol size.
        ## add a pseudocount of 0.1 denominator to avoid division by zero.
        mutate(plasmid_duplicate_percent = plasmid_duplicate_genes/(plasmid_duplicate_genes+chromosomal_duplicate_genes + 0.1)) %>%
        mutate(plasmid_AR_singleton_percent = plasmid_singleton_ARGs/(plasmid_singleton_ARGs+chromosomal_singleton_ARGs + 0.1)) %>%
        mutate(plasmid_AR_duplicate_percent = plasmid_duplicate_ARGs/(plasmid_duplicate_ARGs + chromosomal_duplicate_ARGs + 0.1)) %>%
        ## add a column for label colors in the Figure.
        mutate(annotation_label_color = ifelse(isolates_with_duplicated_ARGs >
                                               expected_isolates_with_duplicated_ARGs,
                                               "red", "black"))
    return(Fig1C.df)
}

Fig1C.df <- make.Fig1C.df(TableS1, ControlTable1, Table2, ControlTable2,duplicate.gene.control.df)

make.Fig1C <- function(Fig1C.df) {

    ## This is for adding a scale for the percent of ARGs on plasmids to Fig1C.
    plasmid_legend_df <- data.frame(plasmid_legend_percent = c(0, 0.25, 0.5, 0.75, 1),
                                    total_isolates = c(10000,10000,10000,10000,10000),
                                    y_pos = c(10^-0.5, 10^-0.25, 10^0, 10^0.25, 10^0.5),
                                    Annotation = c("0%", "25%", "50%", "75%", "100%"))

    
    total_isolates.sum <- sum(Fig1C.df$total_isolates)
    isolates_with_duplicated_ARGs.sum <- sum(Fig1C.df$isolates_with_duplicated_ARGs)
    isolates_with_singleton_ARGs.sum <- sum(Fig1C.df$isolates_with_singleton_ARGs)
    isolates_with_duplicate_genes.sum <- sum(Fig1C.df$isolates_with_duplicate_genes)
    
    Fig1C.color.palette <- scales::viridis_pal()(3)

    Fig1C <- ggplot(Fig1C.df, aes(x=total_isolates,
                                  y=isolates_with_duplicated_ARGs,
                                  label=Annotation)) +
        theme_classic() +
        geom_point(aes(size = plasmid_AR_duplicate_percent * 0.5),
                   color=Fig1C.color.palette[1], alpha=0.2) +
        geom_point(aes(y=isolates_with_singleton_ARGs,
                       size=plasmid_AR_singleton_percent * 0.5),
                   color=Fig1C.color.palette[2],alpha=0.2) +
        geom_point(aes(y=isolates_with_duplicate_genes,
                       size=plasmid_duplicate_percent * 0.5),color="gray",alpha=0.2) +
        geom_line(aes(y=yvals.for.isolates_with_duplicated_ARGs.line),
                  color=Fig1C.color.palette[1]) +
        geom_line(aes(y=yvals.for.isolates_with_singleton_ARGs.line),
                  color=Fig1C.color.palette[2]) +
        geom_line(aes(y=yvals.for.isolates_with_duplicate_genes.line),
                  color="gray") +
        ##        geom_text_repel(size=3) +
        geom_text_repel(size=3,aes(color=annotation_label_color)) +
        scale_color_identity() +
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
        annotate("text", x = 65, y = 1.6, label = "Isolates with duplicated ARGs",
                 angle = 34, color = Fig1C.color.palette[1],size=3) +
        annotate("text", x = 65, y = 13.5, label = "Isolates with singleton ARGs",
                 angle = 34, color = Fig1C.color.palette[2],size=3) +
        annotate("text", x = 65, y = 22, label = "Isolates with duplicated genes",
                 angle = 34, color = "gray",size=3) +
        guides(size=FALSE) +
        ## Add the scale for the percent of the given gene class on plasmids.
        geom_point(data = plasmid_legend_df,
                   aes(y=y_pos,
                       size=plasmid_legend_percent * 0.5),
                   color="black",alpha=0.2) +
        annotate("text", x = 10000, y = 10^0.75, size = 3, label = "Percent on\nplasmids") +
        annotate("text", x = 10000, y = 10^0.5, size = 2, label = "100") +
        annotate("text", x = 10000, y = 10^0.25, size = 2, label = "75") +
        annotate("text", x = 10000, y = 10^0, size = 2, label = "50") +
        annotate("text", x = 10000, y = 10^-0.25, size = 2, label = "25")
    
    return(Fig1C)
}

Fig1C <- make.Fig1C(Fig1C.df)
ggsave(Fig1C,file="../results/Fig1C.pdf",width=6,height=6)


make.grant.Fig1C <- function(Fig1C.df) {

    ## simplify the labels.
    grant.Fig1C.df <- Fig1C.df %>%
        mutate(Annotation = sapply(Annotation, function(x)
            str_replace_all(x, c("Human-host" = "Human",
                                 "Anthropogenic-environment" = "Cities",
                                 "Animal-host" = "Animal",
                                 "Plant-host" = "Plant",
                                 "Fungal-host" = "Fungi"))))
            
        
    
    total_isolates.sum <- sum(grant.Fig1C.df$total_isolates)
    isolates_with_duplicated_ARGs.sum <- sum(grant.Fig1C.df$isolates_with_duplicated_ARGs)
    isolates_with_singleton_ARGs.sum <- sum(grant.Fig1C.df$isolates_with_singleton_ARGs)
    isolates_with_duplicate_genes.sum <- sum(grant.Fig1C.df$isolates_with_duplicate_genes)
    
    Fig1C.color.palette <- scales::viridis_pal()(3)

    Fig1C <- ggplot(grant.Fig1C.df, aes(x=total_isolates,
                                  y=isolates_with_duplicated_ARGs,
                                  label=Annotation)) +
        theme_classic() +
        geom_point(aes(size = plasmid_AR_duplicate_percent * 0.5),
                   color=Fig1C.color.palette[1], alpha=0.2) +
        geom_point(aes(y=isolates_with_singleton_ARGs,
                       size=plasmid_AR_singleton_percent * 0.5),
                   color=Fig1C.color.palette[2],alpha=0.2) +
        geom_point(aes(y=isolates_with_duplicate_genes,
                       size=plasmid_duplicate_percent * 0.5),color="gray",alpha=0.2) +
        geom_line(aes(y=yvals.for.isolates_with_duplicated_ARGs.line),
                  color=Fig1C.color.palette[1]) +
        geom_line(aes(y=yvals.for.isolates_with_singleton_ARGs.line),
                  color=Fig1C.color.palette[2]) +
        geom_line(aes(y=yvals.for.isolates_with_duplicate_genes.line),
                  color="gray") +
        geom_text_repel(size=3) +
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
        annotate("text", x = 85, y = 3, label = "duplicated ARGs",
                 angle = 30, color = Fig1C.color.palette[1],size=3) +
        annotate("text", x = 85, y = 14.5, label = "singleton ARGs",
                 angle = 30, color = Fig1C.color.palette[2],size=3) +
        annotate("text", x = 85, y = 30, label = "duplicated genes",
                 angle = 30, color = "gray",size=3) +
        guides(size=FALSE) 
        
    return(Fig1C)
}

grant.Fig1C <- make.grant.Fig1C(Fig1C.df)
ggsave(grant.Fig1C,file="../results/grant_Fig1C.pdf",width=3,height=2)


################################################################################
## Make a figure similar to Fig1C, using total number of genes as the baseline.

## We're going to make one dataframe with lots of columns for plotting
## Figure S2.

## Then the particular statistics and supplementary tables are going to be
## calculated on relevant subsets of this dataframe (subsetting columns, not rows).

## The number of duplicated AR genes.
duplicate.AR.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Annotation != "Unannotated") %>%
    filter(Annotation != "blank") %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Annotation) %>%
    summarize(AR_duplicates = sum(count)) %>%
    arrange(desc(AR_duplicates))

## The number of duplicated MGE genes.
duplicate.MGE.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Annotation != "Unannotated") %>%
    filter(Annotation != "blank") %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Annotation) %>%
    summarize(MGE_duplicates = sum(count)) %>%
    arrange(desc(MGE_duplicates))

## The number of duplicated genes in each category.
duplicate.genes.count <- duplicate.proteins %>%
    ## remove Unannotated isolates.
    filter(Annotation != "Unannotated") %>%
    filter(Annotation != "blank") %>%
    group_by(Annotation) %>%
    summarize(duplicate_genes = sum(count)) %>%
    arrange(desc(duplicate_genes))

## The number of singleton genes in each category.
singleton.genes.count <- singleton.proteins %>%
    group_by(Annotation) %>%
    summarize(singleton_genes = sum(count)) %>%
    arrange(desc(singleton_genes))
gc() ## free memory

## The number of singleton AR genes.
singleton.AR.genes.count <- singleton.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Annotation) %>%
    summarize(AR_singletons = sum(count)) %>%
    arrange(desc(AR_singletons))
gc() ## free memory

## The number of singleton MGE genes.
singleton.MGE.genes.count <- singleton.proteins %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Annotation) %>%
    summarize(MGE_singletons = sum(count)) %>%
    arrange(desc(MGE_singletons))
gc() ## free memory

## Sum up the totals for duplicate genes and singleton genes.
## This is the baseline column for statistical tests.
total.genes.count <- duplicate.genes.count %>%
    full_join(singleton.genes.count) %>%
    mutate(total_genes = duplicate_genes + singleton_genes) %>%
    select(Annotation, total_genes)

big.gene.analysis.df <- duplicate.AR.genes.count %>%
    full_join(duplicate.MGE.genes.count) %>%
    full_join(duplicate.genes.count) %>%
    full_join(singleton.AR.genes.count) %>%
    full_join(singleton.MGE.genes.count) %>%
    full_join(singleton.genes.count) %>%
    full_join(total.genes.count) %>%
    mutate(AR_duplicates = replace_na(AR_duplicates, 0)) %>%
    arrange(desc(AR_duplicates))

## make Figure S2. The big figure!

make.FigS2.df <- function(big.gene.analysis.df) {
        
    total_genes.sum <- sum(big.gene.analysis.df$total_genes)
    duplicate_ARGs.sum <- sum(big.gene.analysis.df$AR_duplicates)
    duplicate_genes.sum <- sum(big.gene.analysis.df$duplicate_genes)
    AR_singletons.sum <- sum(big.gene.analysis.df$AR_singletons)
    duplicate_MGEs.sum <- sum(big.gene.analysis.df$MGE_duplicates)
    MGE_singletons.sum <- sum(big.gene.analysis.df$MGE_singletons)
    singleton_genes.sum <- sum(big.gene.analysis.df$singleton_genes)
            
    FigS2.df <- big.gene.analysis.df %>%
        ## calculate y-coordinates for line for duplicate genes.
        mutate(yvals.for.duplicate_genes.line = exp(log(total_genes) + log(duplicate_genes.sum) - log(total_genes.sum))) %>%
        ## calculate y-coordinates for line for singleton ARGs.
        mutate(yvals.for.singleton_ARGs.line = exp(log(total_genes) + log(AR_singletons.sum) - log(total_genes.sum))) %>%
        ## calculate y-coordinates for line for singleton genes.
        mutate(yvals.for.singleton_genes.line = exp(log(total_genes) + log(singleton_genes.sum) - log(total_genes.sum))) %>%
        ## calculate y-coordinates for line for duplicate ARGs.
        mutate(yvals.for.duplicated_ARGs.line = exp(log(total_genes) + log(duplicate_ARGs.sum) - log(total_genes.sum))) %>%
    ## calculate y-coordinates for line for singleton MGEs.
        mutate(yvals.for.MGE_singletons.line = exp(log(total_genes) + log(MGE_singletons.sum) - log(total_genes.sum))) %>%
    ## calculate y-coordinates for line for duplicate MGEs.
        mutate(yvals.for.MGE_duplicates.line = exp(log(total_genes) + log(duplicate_MGEs.sum) - log(total_genes.sum)))
    return(FigS2.df)
}

make.FigS2 <- function(FigS2.df) {
    
    FigS2.color.palette <- scales::viridis_pal()(4)

    FigS2 <- ggplot(FigS2.df, aes(x = total_genes,
                                  y = AR_duplicates,
                                  label = Annotation)) +
        theme_classic() +
        scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        xlab("Total genes") +
        ylab("Genes from isolates in given class") +
        guides(size = FALSE) +
        
        geom_point(color = FigS2.color.palette[1], alpha = 0.2) +
        geom_point(aes(y = AR_singletons),
                   color = FigS2.color.palette[2], alpha = 0.2) +
        geom_point(aes(y = MGE_singletons),
                   color = FigS2.color.palette[3], alpha=0.2) +
        geom_point(aes(y = MGE_duplicates),
                   color = FigS2.color.palette[4], alpha = 0.2) +
        geom_point(aes(y = duplicate_genes), color="black", alpha = 0.2) +
        geom_point(aes(y = singleton_genes), color="gray", alpha = 0.2) +
        
        geom_line(aes(y = yvals.for.duplicated_ARGs.line),
                  color = FigS2.color.palette[1]) +
        geom_line(aes(y = yvals.for.singleton_ARGs.line),
                  color = FigS2.color.palette[2]) +
        
        geom_line(aes(y = yvals.for.MGE_singletons.line),
                  color = FigS2.color.palette[3]) +
        geom_line(aes(y = yvals.for.MGE_duplicates.line),
                  color = FigS2.color.palette[4]) +
        geom_line(aes(y = yvals.for.duplicate_genes.line),
                  color = "black") +
        geom_line(aes(y = yvals.for.singleton_genes.line),
                  color = "gray") +
        
        geom_text_repel(size = 2.5) +
        ##geom_text_repel(aes(y = AR_singletons), size = 2.5) +
        ##geom_text_repel(aes(y = duplicate_genes), size = 2.5) +
               
        annotate("text", x = 500000, y = 40, label = "Duplicated ARGs",
                 angle = 25, color = FigS2.color.palette[1],size=3) +
        annotate("text", x = 500000, y = 600 , label = "Singleton ARGs",
                 angle = 25, color = FigS2.color.palette[2],size=3) +
        
        annotate("text", x = 500000, y = 100000, label = "MGE singletons",
                 angle = 25, color = FigS2.color.palette[3],size=3) +
        annotate("text", x = 500000, y = 15000 , label = "MGE duplicates",
                 angle = 25, color = FigS2.color.palette[4],size=3) +

        annotate("text", x = 500000, y = 12000, label = "Duplicated genes",
                 angle = 25, color = "black",size=3) +
        annotate("text", x = 500000, y = 700000, label = "Singleton genes",
                 angle = 25, color = "gray",size=3)

    return(FigS2)
}

FigS2.df <- make.FigS2.df(big.gene.analysis.df)
FigS2 <- make.FigS2(FigS2.df)
ggsave(FigS2,file="../results/FigS2.pdf",width=6,height=6)
#############################

## Use the data in Figure S1 to make Figure 3.
## The point of this figure is to show that ARGs are generally
## enriched on plasmids, especially multi-copy ARGs.

## helper functions to calculate statistics of enrichment on
## chromosomes and plasmids, separately.

calc.chromosomal.ARG.duplicate.enrichment.pvals <- function(raw.Table,
                                                            pval.threshold = 0.01) {

    chromosomal_duplicates.sum <- sum(raw.Table$chromosomal_duplicate_genes)
    chromosomal_ARG_duplicates.sum <- sum(raw.Table$chromosomal_duplicate_ARGs)

    Table <- raw.Table %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = chromosomal_duplicate_ARGs,
                   n = chromosomal_ARG_duplicates.sum,
                   p = chromosomal_duplicate_genes/chromosomal_duplicates.sum
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        mutate(expected_chromosomal_duplicate_ARGs = exp(log(chromosomal_ARG_duplicates.sum) + log(chromosomal_duplicate_genes) - log(chromosomal_duplicates.sum))) %>%
        ## use Benjamini-Hochberg p-value correction,
        ## and only keep relevant columns.
        select(Annotation,
               chromosomal_duplicate_genes,
               chromosomal_duplicate_ARGs,
               expected_chromosomal_duplicate_ARGs,
               corrected.pval) %>%
        ## add a column for label colors in the Figure.
        mutate(annotation_label_color = ifelse(corrected.pval < pval.threshold &&
                                               (chromosomal_duplicate_ARGs >
                                                expected_chromosomal_duplicate_ARGs),
                                               "red", "black"))
    return(Table)
}

calc.plasmid.ARG.duplicate.enrichment.pvals <- function(raw.Table,
                                                        pval.threshold = 0.01) {

    plasmid_duplicates.sum <- sum(raw.Table$plasmid_duplicate_genes)
    plasmid_ARG_duplicates.sum <- sum(raw.Table$plasmid_duplicate_ARGs)

    Table <- raw.Table %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = plasmid_duplicate_ARGs,
                   n = plasmid_ARG_duplicates.sum,
                   p = plasmid_duplicate_genes/plasmid_duplicates.sum
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        mutate(expected_plasmid_duplicate_ARGs = exp(log(plasmid_ARG_duplicates.sum) + log(plasmid_duplicate_genes) - log(plasmid_duplicates.sum))) %>%
        ## use Benjamini-Hochberg p-value correction,
        ## and only keep relevant columns.
        select(Annotation,
               plasmid_duplicate_genes,
               plasmid_duplicate_ARGs,
               expected_plasmid_duplicate_ARGs,
               corrected.pval) %>%
        ## add a column for label colors in the Figure.
        mutate(annotation_label_color = ifelse(corrected.pval < pval.threshold &&
                                               (plasmid_duplicate_ARGs >
                                                expected_plasmid_duplicate_ARGs),
                                               "red", "black"))    
    return(Table)
}

calc.chromosomal.ARG.singleton.enrichment.pvals <- function(raw.Table,
                                                            pval.threshold = 0.01) {

    chromosomal_singleton.sum <- sum(raw.Table$chromosomal_singleton_genes)
    chromosomal_ARG_singleton.sum <- sum(raw.Table$chromosomal_singleton_ARGs)

    Table <- raw.Table %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = chromosomal_singleton_ARGs,
                   n = chromosomal_ARG_singleton.sum,
                   p = chromosomal_singleton_genes/chromosomal_singleton.sum
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        mutate(expected_chromosomal_singleton_ARGs = exp(log(chromosomal_ARG_singleton.sum) + log(chromosomal_singleton_genes) - log(chromosomal_singleton.sum))) %>%
        ## use Benjamini-Hochberg p-value correction,
        ## and only keep relevant columns.
        select(Annotation,
               chromosomal_singleton_genes,
               chromosomal_singleton_ARGs,
               expected_chromosomal_singleton_ARGs,
               corrected.pval) %>%
        ## add a column for label colors in the Figure.
        mutate(annotation_label_color = ifelse(corrected.pval < pval.threshold &&
                                               (chromosomal_singleton_ARGs >
                                                expected_chromosomal_singleton_ARGs),
                                               "red", "black"))
    return(Table)
}

calc.plasmid.ARG.singleton.enrichment.pvals <- function(raw.Table,
                                                        pval.threshold = 0.01) {

    plasmid_singleton.sum <- sum(raw.Table$plasmid_singleton_genes)
    plasmid_ARG_singleton.sum <- sum(raw.Table$plasmid_singleton_ARGs)

    Table <- raw.Table %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = plasmid_singleton_ARGs,
                   n = plasmid_ARG_singleton.sum,
                   p = plasmid_singleton_genes/plasmid_singleton.sum
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        mutate(expected_plasmid_singleton_ARGs = exp(log(plasmid_ARG_singleton.sum) + log(plasmid_singleton_genes) - log(plasmid_singleton.sum))) %>%
        ## use Benjamini-Hochberg p-value correction,
        ## and only keep relevant columns.
        select(Annotation,
               plasmid_singleton_genes,
               plasmid_singleton_ARGs,
               expected_plasmid_singleton_ARGs,
               corrected.pval) %>%
        ## add a column for label colors in the Figure.
        mutate(annotation_label_color = ifelse(corrected.pval < pval.threshold &&
                                               (plasmid_singleton_ARGs >
                                                expected_plasmid_singleton_ARGs),
                                               "red", "black"))
    return(Table)
}


make.Fig3.df <- function(ControlTable2, Table2) {

    Fig3.data <- full_join(ControlTable2, Table2) %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(c("Human-host","Livestock","Animal-host",
                                  "Anthropogenic-environment", "Food", "Freshwater",
                                  "Agriculture", "Sediment", "Soil", "Plant-host",
                                  "Marine","Terrestrial", "Fungal-host"))))

    ## This follows the logic in the function make.FigS2.df--
    ## See Figure 1B for the intuitive explanation.
    ## basically, calculate the expected fraction in the given category,
    ## based on the slope, which is the (genes in category)/(total genes in category),
    ## under the null hypothesis that ARGs follow the slope (the sampling distribution).
    
    chromosomal_duplicates.sum <- sum(Fig3.data$chromosomal_duplicate_genes)
    chromosomal_ARG_duplicates.sum <- sum(Fig3.data$chromosomal_duplicate_ARGs)
    plasmid_duplicates.sum <- sum(Fig3.data$plasmid_duplicate_genes)
    plasmid_ARG_duplicates.sum <- sum(Fig3.data$plasmid_duplicate_ARGs)

    chromosomal_singleton.sum <- sum(Fig3.data$chromosomal_singleton_genes)
    chromosomal_ARG_singleton.sum <- sum(Fig3.data$chromosomal_singleton_ARGs)
    plasmid_singleton.sum <- sum(Fig3.data$plasmid_singleton_genes)
    plasmid_ARG_singleton.sum <- sum(Fig3.data$plasmid_singleton_ARGs)

    Fig3.df <- Fig3.data %>%
        ## calculate y-coordinates for chromosomal duplicate ARGs.
        mutate(yvals.for.chromosomal_ARG_duplicates.line = exp(log(chromosomal_ARG_duplicates.sum) + log(chromosomal_duplicate_genes) - log(chromosomal_duplicates.sum))) %>%
        ## calculate y-coordinates for plasmid duplicate ARGs.
        mutate(yvals.for.plasmid_ARG_duplicates.line = exp(log(plasmid_ARG_duplicates.sum) + log(plasmid_duplicate_genes) - log(plasmid_duplicates.sum))) %>%
        ## calculate y-coordinates for chromosomal singleton ARGs.
        mutate(yvals.for.chromosomal_ARG_singletons.line = exp(log(chromosomal_ARG_singleton.sum) + log(chromosomal_singleton_genes) - log(chromosomal_singleton.sum))) %>%
        ## calculate y-coordinates for plasmid singleton ARGs.
        mutate(yvals.for.plasmid_ARG_singletons.line = exp(log(plasmid_ARG_singleton.sum) + log(plasmid_singleton_genes) - log(plasmid_singleton.sum)))
    
    return(Fig3.df)
}


Fig3.df <- make.Fig3.df(ControlTable2, Table2)

## calculate formal statistics, and use to add columns for coloring annotations.
chromosome.ARG.duplicate.pvals <- Fig3.df %>%
    calc.chromosomal.ARG.duplicate.enrichment.pvals()    
plasmid.ARG.duplicate.pvals <- Fig3.df %>%
    calc.plasmid.ARG.duplicate.enrichment.pvals()
chromosome.ARG.singleton.pvals <- Fig3.df %>%
    calc.chromosomal.ARG.singleton.enrichment.pvals()
plasmid.ARG.singleton.pvals <- Fig3.df %>%
    calc.plasmid.ARG.singleton.enrichment.pvals()

Fig3A <- ggplot(Fig3.df, aes(x=chromosomal_duplicate_genes,
                                  y=chromosomal_duplicate_ARGs,
                                  label=Annotation)) +
        theme_classic() +
        geom_point(alpha=0.2) +                  
        geom_line(aes(y=yvals.for.chromosomal_ARG_duplicates.line)) +
    geom_text_repel(data=chromosome.ARG.duplicate.pvals,
                    aes(color=annotation_label_color), size=3) +
    scale_color_identity() +
        scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        xlab("Chromosomal multi-copy proteins") +
        ylab("Chromosomal multi-copy ARGs")

Fig3B <- ggplot(Fig3.df, aes(x=plasmid_duplicate_genes,
                                  y=plasmid_duplicate_ARGs,
                                  label=Annotation)) +
    theme_classic() +
    geom_point(alpha=0.2) +                  
    geom_line(aes(y=yvals.for.plasmid_ARG_duplicates.line)) +
    geom_text_repel(data=plasmid.ARG.duplicate.pvals,
                    aes(color=annotation_label_color), size=3) +
    scale_color_identity() +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    xlab("Plasmid multi-copy proteins") +
    ylab("Plasmid multi-copy ARGs")

Fig3C <- ggplot(Fig3.df, aes(x=chromosomal_singleton_genes,
                             y=chromosomal_singleton_ARGs,
                             label=Annotation)) +
    theme_classic() +
    geom_point(alpha=0.2) +                  
    geom_line(aes(y=yvals.for.chromosomal_ARG_singletons.line)) +
    geom_text_repel(data=chromosome.ARG.singleton.pvals,
                    aes(color=annotation_label_color), size=3) +
    scale_color_identity() +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    xlab("Chromosomal single-copy proteins") +
    ylab("Chromosomal single-copy ARGs")

Fig3D <- ggplot(Fig3.df, aes(x=plasmid_singleton_genes,
                             y=plasmid_singleton_ARGs,
                             label=Annotation)) +
    theme_classic() +
    geom_point(alpha=0.2) +                  
    geom_line(aes(y=yvals.for.plasmid_ARG_singletons.line)) +
    geom_text_repel(data=plasmid.ARG.singleton.pvals,
                    aes(color=annotation_label_color), size=3) +
    scale_color_identity() +    
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    xlab("Plasmid single-copy proteins") +
    ylab("Plasmid single-copy ARGs")

Fig3 <- plot_grid(Fig3A, Fig3B, Fig3C, Fig3D, labels = c("A","B","C","D"), nrow=2)
ggsave("../results/Fig3.pdf", Fig3)



#############################
## Duplication Index Ratio calculations.

duplication.index.df <- big.gene.analysis.df %>%
    ## calculate the ratio of duplicated ARGs to singleton ARGs.
    mutate(DI.ARGs = AR_duplicates/AR_singletons)  %>%
    mutate(DI.all = duplicate_genes/singleton_genes) %>%
    mutate(Annotation = factor(Annotation,
                                      levels=rev(big.gene.analysis.df$Annotation))) %>%
    ## make a bar graph comparing DI for ARGs to DI for all genes.
    pivot_longer(c(DI.ARGs, DI.all), names_to = "DI.type", values_to = "DI")

DI.ARGs.to.all.fig <- ggplot(duplication.index.df,
                       aes(x = Annotation, y = DI,fill=DI.type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_classic() +
    ## make a better legend.
    scale_fill_manual(values=c("gray","black"),
                      name = "",
                      labels=c("All genes", "ARGs")) +
    theme(axis.text.x  = element_text(angle=90,vjust=0.2,hjust=0.95)) +
    ylab("Duplication Index") +
    theme(legend.position="top") +
    xlab("Annotation")

ggsave("../results/DI-index.pdf", DI.ARGs.to.all.fig)

#############################
## Supplementary Table S2.
## Enrichment/deletion analysis of AR genes using total genes,
## rather than number of isolates as in Supplementary Table S1.

## the expected number of duplicated AR genes in each category.
calc.expected.AR.duplicates <- function(raw.Table) {
    total.genes.sum <- sum(raw.Table$total_genes)
    total.AR.duplicates <- sum(raw.Table$AR_duplicates)
    Table <- raw.Table %>%
        mutate(expected_AR_duplicates = total.AR.duplicates * total_genes/total.genes.sum)
    return(Table)
}

## Seventh column: p-value for enrichment/depletion of duplicated AR genes in each category.
calc.AR.duplicate.enrichment.pvals <- function(raw.Table) {

    total.genes.sum <- sum(raw.Table$total_genes)
    total.AR.duplicates <- sum(raw.Table$AR_duplicates)

    Table <- raw.Table %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = AR_duplicates,
                   n = total.AR.duplicates,
                   p = total_genes/total.genes.sum
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        ## use Benjamini-Hochberg p-value correction.
        select(-binom.test.pval) ## drop original p-value after the correction.
    
    return(Table)
}

TableS2 <- big.gene.analysis.df %>%
    select(Annotation, AR_duplicates, total_genes) %>%
    calc.expected.AR.duplicates() %>% 
    calc.AR.duplicate.enrichment.pvals() %>%
    mutate(deviation.from.expected = AR_duplicates - expected_AR_duplicates) %>%
    arrange(desc(deviation.from.expected)) %>%
    select(-deviation.from.expected)


## write Supplementary Table S2 to file.
write.csv(x=TableS2,file="../results/TableS2.csv")

############################################################
## Positive control: Examine the distribution of ARGs that have NOT duplicated.
## This shows that Soil and Agricultural isolates are highly enriched in
## singleton ARGs.

## calculate the expected number of singleton AR genes in each category.
calc.expected.AR.singletons <- function(raw.Table) {
    total.genes.sum <- sum(raw.Table$total_genes)
    total.AR.singletons <- sum(raw.Table$AR_singletons)
    Table <- raw.Table %>%
        mutate(expected_AR_singletons = exp(log(total.AR.singletons) + log(total_genes) - log(total.genes.sum)))
    return(Table)
}

## p-value for enrichment/depletion of singleton AR genes in each category.
calc.AR.singleton.enrichment.pvals <- function(raw.Table) {
    total.genes.sum <- sum(raw.Table$total_genes)
    total.AR.singletons <- sum(raw.Table$AR_singletons)

    Table <- raw.Table %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = AR_singletons,
                   n = total.AR.singletons,
                   p = total_genes/total.genes.sum
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        ## use Benjamini-Hochberg p-value correction.
        select(-binom.test.pval) ## drop original p-value after the correction.
    
    return(Table)
}

TableS3 <- big.gene.analysis.df %>%
    select(Annotation, AR_singletons, total_genes) %>%
    calc.expected.AR.singletons() %>%
    calc.AR.singleton.enrichment.pvals() %>%
    mutate(deviation.from.expected = AR_singletons - expected_AR_singletons) %>%
    arrange(desc(deviation.from.expected)) %>%
    select(-deviation.from.expected)

write.csv(x=TableS3, file="../results/TableS3.csv")
############################################################
## Positive control: Examine the distribution of duplicated genes.

## The distribution of total genes is a terrible null model for the distribution
## of duplicated genes. The number of duplicated genes is not proportional
## to the total number of genes.

## calculate the expected number of duplicated in each category.
calc.expected.duplicated <- function(raw.Table) {
    total.genes.sum <- sum(raw.Table$total_genes)
    total.duplicates <- sum(raw.Table$duplicate_genes)
    Table <- raw.Table %>%
        mutate(expected_duplicates = exp(log(total.duplicates) + log(total_genes) - log(total.genes.sum)))
    return(Table)
}

## p-value for enrichment/depletion of duplicate genes in each category.
calc.duplicate.enrichment.pvals <- function(raw.Table) {
    total.genes.sum <- sum(raw.Table$total_genes)
    total.duplicates <- sum(raw.Table$duplicate_genes)

    Table <- raw.Table %>%
        rowwise() %>%
        mutate(binom.test.pval = binom.test(
                   x = duplicate_genes,
                   n = total.duplicates,
                   p = total_genes/total.genes.sum
               )$p.value) %>%
        mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
        ## use Benjamini-Hochberg p-value correction.
        select(-binom.test.pval) ## drop original p-value after the correction.
    
    return(Table)
}

duplicate.table.test <- big.gene.analysis.df %>%
    select(Annotation, duplicate_genes, total_genes) %>%
    calc.expected.duplicated() %>%
    calc.duplicate.enrichment.pvals() %>%
    mutate(deviation.from.expected = duplicate_genes - expected_duplicates) %>%
    arrange(desc(deviation.from.expected)) %>%
    select(-deviation.from.expected)

######################################################
## Control for Table S2 that Lingchong asked me to make.
## examine the number of types of duplicate genes in each category,
## and the average num.

## the number of duplicate gene types.
duplicated.gene.type.count <- duplicate.proteins %>%
    group_by(Annotation) %>%
    ## each row corresponds to a type of duplicated gene.
    summarize(duplicate_gene_types = n(),
              mean.duplicate.num = sum(count)/n()) %>%
    arrange(desc(duplicate_gene_types))

## the number of duplicated MGE gene types
duplicated.MGE.type.count <- duplicate.proteins %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Annotation) %>%
    ## each row corresponds to a type of duplicated MGE.
    summarize(duplicate_MGE_types = n(),
              mean.MGE.duplicate.num = sum(count)/n()) %>%
    arrange(desc(duplicate_MGE_types))

## the number of duplicated AR gene types.
duplicated.ARG.type.count <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Annotation) %>%
    ## each row corresponds to a type of duplicated ARG.
    summarize(duplicate_ARG_types = n(),
              mean.ARG.duplicate.num = sum(count)/n()) %>%
    arrange(desc(duplicate_ARG_types))

Control.for.TableS2 <- duplicated.gene.type.count %>% 
    left_join(duplicated.MGE.type.count) %>%
    left_join(duplicated.ARG.type.count) %>%
    mutate(duplicate_ARG_types = replace_na(duplicate_ARG_types, 0)) %>%
    mutate(mean.ARG.duplicate.num = replace_na(mean.ARG.duplicate.num, 0)) %>%
    arrange(desc(duplicate_ARG_types))

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
## include genes from Unannotated genomes.

## get the match between sequences and product annotations.
duplicate.protein.annotations <- all.duplicate.proteins %>%
    select(product, sequence) %>%
    distinct() %>%
## PROBLEM: often, an IDENTICAL protein sequence will have multiple 'product'
## annotations! See the solution here.
## https://stackoverflow.com/questions/19944334/extract-rows-for-the-first-occurrence-of-a-variable-in-a-data-frame
group_by(sequence) %>% 
  filter(product == max(product)) %>% 
  slice(1) %>% # takes the first occurrence if there is a tie
  ungroup()
    
HGT.candidates <- all.duplicate.proteins %>%
    group_by(sequence) %>%
    summarize(number.of.genomes = n(),
              total.copies = sum(count),
              chromosome.copies = sum(chromosome_count),
              plasmid.copies = sum(plasmid_count)) %>%
    filter(number.of.genomes > 1) %>%
    left_join(duplicate.protein.annotations) %>%
    arrange(desc(number.of.genomes))

HGT.candidate.summary <- HGT.candidates %>%
    select(-sequence)

MGE.HGT.candidates <- HGT.candidate.summary %>%
    filter(str_detect(.$product,IS.keywords))

non.MGE.HGT.candidates <- HGT.candidate.summary %>%
    filter(!str_detect(.$product,IS.keywords))

AR.HGT.candidates <- HGT.candidate.summary %>%
    filter(str_detect(.$product,antibiotic.keywords))

plasmid.only.HGT.candidates <- HGT.candidate.summary %>%
    filter(chromosome.copies == 0)

non.MGE.plasmid.only.HGT.candidates <- plasmid.only.HGT.candidates %>%
    filter(!str_detect(.$product,IS.keywords))

################################################################################

## Calculate TF-IDF (Term Frequency times Inverse Document Frequency)
## for each ecological category, using protein sequences.
## see "Mining of Massive Datasets"
## and https://en.wikipedia.org/wiki/Tf%E2%80%93idf.

## References for R package tidytext:
## https://www.tidytextmining.com/
## https://www.tidytextmining.com/tfidf.html

## TF-IDF is good at finding proteins/terms that are
## specific to particular ecological annotations.

## let's look at the annotations and sequences of high frequency duplicated
## proteins in each environment, after removing MGEs, and EF-Tu.

## this function filters proteins by manual annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
.make.annotation.freq.table <- function(data.df, manual.annot.string, remove.MGEs = TRUE) {
    df <- data.df %>%
        filter(Annotation == manual.annot.string) %>%
        filter(!str_detect(.$product, EFTu.keywords))

    if (remove.MGEs)
        df <- df %>% filter(!str_detect(.$product,IS.keywords))
    
    df %>% group_by(Annotation, product) %>%
        summarize(annotation.count = n()) %>%
        arrange(desc(annotation.count)) %>%
        ungroup() %>%
        mutate(total.duplicate.proteins = sum(annotation.count))
}

## IMPORTANT: These are the functions that are actually used.
make.dup.annotation.freq.table <- partial(.f = .make.annotation.freq.table, duplicate.proteins)
make.sing.annotation.freq.table <- partial(.f = .make.annotation.freq.table, singleton.proteins)

make.dup.with.MGEs.annotation.freq.table <- partial(.f = .make.annotation.freq.table,
                                                    duplicate.proteins, remove.MGEs = FALSE)

make.sing.with.MGEs.annotation.freq.table <- partial(.f = .make.annotation.freq.table,
                                                     singleton.proteins, remove.MGEs = FALSE)


## let's make a big table of product annotations per annotation category,
## for TF-IDF analysis.

## The analysis here closely follows the text mining example here:
## https://www.tidytextmining.com/tfidf.html

big.dup.prot.annotation.freq.table <- map_dfr(unique(duplicate.proteins$Annotation),
                                              .f = make.dup.annotation.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all annotations in the table.
    mutate(total.annotation.count = sum(annotation.count))

dup.prot.annotation.tf_idf <- big.dup.prot.annotation.freq.table %>%
  bind_tf_idf(product, Annotation, annotation.count) %>%
  arrange(desc(tf_idf))

top.dup.prot.annotation.tf_idf <- dup.prot.annotation.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup()

dup.prot.annotation.tf_idf.plot <- top.dup.prot.annotation.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~Annotation, ncol = 2, scales = "free") +
  labs(x = "tf-idf", y = NULL)

ggsave("../results/duplicate-protein-annotation-TF-IDF.pdf",
       dup.prot.annotation.tf_idf.plot,
       height=21,width=21)


## repeat the analysis, but keep MGE sequences this time.
## this is a valuable comparison, because it shows how MGEs completely dominate the
## functional annotation of multicopy proteins across ecological categories.

big.dup.with.MGEs.prot.annotation.freq.table <- map_dfr(unique(duplicate.proteins$Annotation),
                                              .f = make.dup.with.MGEs.annotation.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all annotations in the table.
    mutate(total.annotation.count = sum(annotation.count))

dup.with.MGEs.prot.annotation.tf_idf <- big.dup.with.MGEs.prot.annotation.freq.table %>%
  bind_tf_idf(product, Annotation, annotation.count) %>%
  arrange(desc(tf_idf))

top.dup.with.MGEs.prot.annotation.tf_idf <- dup.with.MGEs.prot.annotation.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup()

dup.with.MGEs.prot.annotation.tf_idf.plot <- top.dup.with.MGEs.prot.annotation.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~Annotation, ncol = 2, scales = "free") +
  labs(x = "tf-idf", y = NULL)

ggsave("../results/duplicate-with-MGEs-protein-annotation-TF-IDF.pdf",
       dup.with.MGEs.prot.annotation.tf_idf.plot,
       height=21,width=21)



#### Now repeat this analysis, but with singleton proteins.
#### This is a really valuable comparison. It is clear that
#### the annotations of duplicate proteins carry far more ecological information
#### than the annotations of singleton proteins.

big.sing.prot.annotation.freq.table <- map_dfr(unique(singleton.proteins$Annotation),
                                              .f = make.sing.annotation.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all annotations in the table.
    mutate(total.annotation.count = sum(annotation.count))

sing.prot.annotation.tf_idf <- big.sing.prot.annotation.freq.table %>%
  bind_tf_idf(product, Annotation, annotation.count) %>%
  arrange(desc(tf_idf))

top.sing.prot.annotation.tf_idf <- sing.prot.annotation.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup()

sing.prot.annotation.tf_idf.plot <- top.sing.prot.annotation.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~Annotation, ncol = 2, scales = "free") +
  labs(x = "tf-idf", y = NULL)

ggsave("../results/singleton-protein-annotation-TF-IDF.pdf",
       sing.prot.annotation.tf_idf.plot,
       height=21,width=21)

## now repeat the singleton analysis, keeping MGEs.

big.sing.with.MGEs.prot.annotation.freq.table <- map_dfr(unique(singleton.proteins$Annotation),
                                              .f = make.sing.with.MGEs.annotation.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all annotations in the table.
    mutate(total.annotation.count = sum(annotation.count))

sing.with.MGEs.prot.annotation.tf_idf <- big.sing.with.MGEs.prot.annotation.freq.table %>%
  bind_tf_idf(product, Annotation, annotation.count) %>%
  arrange(desc(tf_idf))

top.sing.with.MGEs.prot.annotation.tf_idf <- sing.with.MGEs.prot.annotation.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup()

sing.with.MGEs.prot.annotation.tf_idf.plot <- top.sing.with.MGEs.prot.annotation.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~Annotation, ncol = 2, scales = "free") +
  labs(x = "tf-idf", y = NULL)

ggsave("../results/singleton-with-MGEs-protein-annotation-TF-IDF.pdf",
       sing.with.MGEs.prot.annotation.tf_idf.plot,
       height=21,width=21)


#########
## this function filters duplicate proteins by annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
make.seq.freq.table <- function(annot.string) {
    duplicate.proteins %>%
        filter(Annotation == annot.string) %>%
        filter(!str_detect(.$product,IS.keywords)) %>%
        filter(!str_detect(.$product, EFTu.keywords)) %>%
        group_by(Annotation, sequence) %>%
        summarize(seq.count = n()) %>%
        arrange(desc(seq.count)) %>%
        ungroup() %>%
        mutate(total.annotation.duplicate.seqs = sum(seq.count))
}

## let's make a big table of product annotations per annotation category,
## for TF-IDF analysis.

## The analysis here closely follows the text mining example here:
## https://www.tidytextmining.com/tfidf.html

big.dup.prot.seq.freq.table <- map_dfr(unique(duplicate.proteins$Annotation),
                                              .f = make.seq.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all seqs in the table.
    mutate(total.seqs = sum(seq.count))

## let's make canonical protein annotations (allow only one annotation per sequence).
canonical.dup.prot.seq.annotations <- duplicate.proteins %>%
    select(sequence, product) %>%
    group_by(sequence) %>%
    ## take the first product annotation as the canonical annotation.
    filter(row_number() == 1) %>%
    ungroup()

dup.prot.seq.tf_idf <- big.dup.prot.seq.freq.table %>%
  bind_tf_idf(sequence, Annotation, seq.count) %>%
    arrange(desc(tf_idf)) %>%
    ## now merge the canonical sequence annotations
    left_join(canonical.dup.prot.seq.annotations)

top.dup.prot.seq.tf_idf <- dup.prot.seq.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup()

dup.prot.seq.tf_idf.plot <- top.dup.prot.seq.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
    geom_col(show.legend = TRUE) +
    facet_wrap(~Annotation, ncol = 2, scales = "free") +
    labs(x = "tf-idf", y = NULL)

ggsave("../results/duplicate-protein-seq-TF-IDF.pdf",
       dup.prot.seq.tf_idf.plot,
       height=21,width=21)

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
make.plas.annotation.freq.table <- function(annot.string) {
    duplicate.proteins %>%
        filter(plasmid_count >= 1) %>%
        filter(Annotation == annot.string) %>%
        filter(!str_detect(.$product,IS.keywords)) %>%
        group_by(Annotation, product) %>%
        summarize(annotation.count = n()) %>%
        arrange(desc(annotation.count))
}


big.dup.plas.annotation.freq.table <- map_dfr(unique(duplicate.proteins$Annotation),
                                              .f = make.plas.annotation.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all annotations in the table.
    mutate(total.annotation.count = sum(annotation.count))

dup.plas.annotation.tf_idf <- big.dup.plas.annotation.freq.table %>%
  bind_tf_idf(product, Annotation, annotation.count) %>%
    arrange(desc(tf_idf))

top.dup.plas.annotation.tf_idf <- dup.plas.annotation.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup()

dup.plas.annotation.tf_idf.plot <- top.dup.plas.annotation.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~Annotation, ncol = 2, scales = "free") +
  labs(x = "tf-idf", y = NULL)

ggsave("../results/duplicate-plasmid-protein-annotation-TF-IDF.pdf",
       dup.plas.annotation.tf_idf.plot,
       height=21,width=21)


## this function filters duplicate proteins by annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
make.on.plas.seq.freq.table <- function(annot.string) {
    duplicate.proteins %>%
        filter(plasmid_count >= 1) %>%
        filter(Annotation == annot.string) %>%
        filter(!str_detect(.$product,IS.keywords)) %>%
        group_by(Annotation, sequence) %>%
        summarize(seq.count = n()) %>%
        filter(seq.count > 1) %>%
        arrange(desc(seq.count))
}


big.dup.plas.seq.freq.table <- map_dfr(unique(duplicate.proteins$Annotation),
                                              .f = make.on.plas.seq.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all seqs in the table.
    mutate(total.seqs = sum(seq.count))

dup.plas.seq.tf_idf <- big.dup.plas.seq.freq.table %>%
  bind_tf_idf(sequence, Annotation, seq.count) %>%
    arrange(desc(tf_idf)) %>%
    ## now merge the canonical sequence annotations
    left_join(canonical.dup.prot.seq.annotations)

top.dup.plas.seq.tf_idf <- dup.plas.seq.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup()

dup.plas.seq.tf_idf.plot <- top.dup.plas.seq.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
    geom_col(show.legend = TRUE) +
    facet_wrap(~Annotation, ncol = 2, scales = "free") +
    labs(x = "tf-idf", y = NULL)

ggsave("../results/duplicate-plasmid-seq-TF-IDF.pdf",
       dup.plas.seq.tf_idf.plot,
       height=21,width=21)

## Figure 4. annotations of multi-copy proteins are more informative
## about ecology than the annotations of single-copy proteins.

best.dup.prot.annotation.tf_idf <- dup.prot.annotation.tf_idf %>%
    filter(Annotation %in% c("Agriculture", "Anthropogenic-environment",
                             "Human-host")) %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 5) %>%
    ungroup()

best.sing.prot.annotation.tf_idf <- sing.prot.annotation.tf_idf %>%
    filter(Annotation %in% c("Agriculture", "Anthropogenic-environment",
                             "Human-host")) %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 5) %>%
    ungroup()


Fig4A <- ggplot(best.dup.prot.annotation.tf_idf,
                aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
    geom_col(show.legend = FALSE) +
    labs(x = "tf-idf", y = NULL) +
    theme_classic() +
    xlim(0,0.04) +
    facet_wrap(.~Annotation, ncol=1, scales = "free_y") +
    ggtitle("Most informative multi-copy protein annotations")

Fig4B <- ggplot(best.sing.prot.annotation.tf_idf,
                aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
    geom_col(show.legend = FALSE) +
    labs(x = "tf-idf", y = NULL) +
    theme_classic() +
    xlim(0,0.04) +
    facet_wrap(.~Annotation, ncol=1, scales = "free_y") +
    ggtitle("Most informative single-copy protein annotations")


Fig4 <- plot_grid(Fig4A, Fig4B, ncol = 1, labels = c('A','B'), align = "v")

ggsave("../results/Fig4.pdf", Fig4, width=8, height = 8)


##########################################
## NOTES AND IDEAS

## signal decomposition algorithms: non-negative matrix factorization,
## ICA, etc.
## represent strains by duplicated genes. Then factorize those strains
## into non-negative combinations of sequences/annotations
## and their counts. for a classifier, weight sequences/annotations that
## are most informative of particular environments.

## if I filter out MGEs, then I probably get more insight into WHAT
## functions are being selected.
## BUT MGEs probably carry significant information about niches!
## assume that gene flow networks are more tightly connected within a niche
## in comparison to between niches (since higher probability of interaction).
## This idea goes back at least to Smillie et al. 2011 in Nature from Eric Alm's group.
