## ARG-duplication-analysis.R by Rohan Maddamsetti.

## analyse the distribution of AR genes on chromosomes versus plasmids in
## fully-sequenced genomes and plasmids in the NCBI Nucleotide database.

## This is a rewrite of Aim1-analysis.R, using a different statistical framework.

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

################################################################################
## Regular expressions used in this analysis.

## match MGE genes using the following keywords in the "product" annotation
IS.keywords <- "IS|transposon|Transposase|transposase|Transposable|transposable|hypothetical protein|Phage|phage|integrase|Integrase|tail|intron|Mobile|mobile|antitoxin|toxin|capsid|plasmid|Plasmid|conjug"

## Elongation Factor Tu (2 copies in most bacteria).
## \\b is a word boundary.
## see: https://stackoverflow.com/questions/62430498/detecting-whole-words-using-str-detect-in-r
EFTu.keywords <- "\\bTu | Tu\\b|-Tu\\b"

## antibiotic-specific keywords.
antibiotic.keywords <- "lactamase|chloramphenicol|quinolone|antibiotic resistance|tetracycline|VanZ"

## unknown protein keywords
unknown.protein.keywords <- "unknown|Unknown|hypothetical|Hypothetical|Uncharacterized|Uncharacterised|uncharacterized|uncharacterised|DUF|unknow|putative protein in bacteria|Unassigned|unassigned"

## The regular expressions used by Zeevi et al. (2019).
## These are not used in this analysis, but nice to have on hand.
## Transposon: ‘transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon’
## plasmid: ‘relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation’
## phage: ‘capsid|phage|tail|head|tape measure|antiterminatio’
## other HGT mechanisms: ‘integrase|excision\S*|exo- nuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip’
## antibiotic resistance: ‘azole resistance|antibiotic resistance|TetR|tetracycline resistance|VanZ|betalactam\S*|beta-lactam|antimicrob\S*|lantibio\S*’.


################################################################################
## Set up the key data structures for the analysis:
## gbk.annotation, in particular.

## import the 17GB file containing all proteins, including singletons.
## I can save a ton of memory if I don't import the sequence column,
## and by using the data.table package for import.
all.proteins <- data.table::fread("../results/all-proteins.csv",
                                  drop="sequence")

## annotate source sequences as plasmid or chromosome.
episome.database <- read.csv("../results/chromosome-plasmid-table.csv") %>%
    as_tibble()

gbk.annotation <- read.csv(
    "../results/computationally-annotated-gbk-annotation-table.csv") %>%
    as_tibble() %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## get species name annotation from episome.database.
    left_join(episome.database) %>%
    ## Annotate the genera.
    mutate(Genus = stringr::word(Organism, 1)) %>%
    ## CRITICAL STEP: remove the NCBI_Nucleotide_Accession and SequenceType columns.
    ## This is absolutely critical, otherwise each row is duplicated for every
    ## chromosome and plasmid, breaking the invariant that each row refers to one sequence,
    ## when we add this annotation to duplicate.proteins and singleton.proteins.
    select(-NCBI_Nucleotide_Accession, -SequenceType) %>%
    ## and we have to explicitly remove redundant rows now.
    distinct()

## Some strains in chromosome-and-plasmid-table.csv and
## gbk-annotation-table.csv are missing from
## all-proteins.csv
## These should be the genomes that do not have
## CDS annotated in their GFF annotation.
## list the 1,064 strains missing from the singletons data.
missing.ones <- gbk.annotation %>%
    filter(!(Annotation_Accession %in% all.proteins$Annotation_Accession))
write.csv(missing.ones, file= "../results/strains-without-proteins.csv")

## CRITICAL STEP: remove all genomes that do not have proteins annotated.
gbk.annotation <- anti_join(gbk.annotation, missing.ones) %>%
    ## And now remove all Unannotated genomes, since these are not analyzed
    ## at all in this first paper.
    filter(Annotation != "Unannotated") %>%
    ## and remove any strains (although none should fall in this category)
    ## that were not annotated by annotate-ecological-category.py.
    filter(Annotation != "blank")

## return the first column for several tables.
## shows the number of isolates in each category.
make.isolate.totals.col <- function(gbk.annotation) {
    isolate.totals <- gbk.annotation %>%
        group_by(Annotation) %>%
        summarize(total_isolates = n()) %>%
        arrange(desc(total_isolates))
    return(isolate.totals)
}


## This vector is used for ordering axes in figures and tables.
order.by.total.isolates <- make.isolate.totals.col(gbk.annotation)$Annotation

## and filter episome.database to be consistent with gbk.annotation.
episome.database <- episome.database %>%
    filter(Annotation_Accession %in% gbk.annotation$Annotation_Accession)

## now get the singleton protein by filtering.
singleton.proteins <- all.proteins %>%
    filter(count == 1) %>%
    inner_join(gbk.annotation)

## read in duplicate proteins with sequences, using a separate file.
## I want the sequence column for the duplicate genes,
## but not for the singletons, to save memory.
duplicate.proteins <- read.csv("../results/duplicate-proteins.csv") %>%
    ## now merge with gbk annotation.
    inner_join(gbk.annotation)

## For Teng (and myself), let's make a table of uncategorized duplicate proteins.
unmatched.duplicate.proteins <- duplicate.proteins %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    filter(!str_detect(.$product,unknown.protein.keywords)) %>%
    filter(!str_detect(.$product,EFTu.keywords)) %>%
    filter(!str_detect(.$product,antibiotic.keywords)) %>%
    as_tibble()
write.csv(unmatched.duplicate.proteins,
          file= "../results/non-MGE-non-ARG-duplicate-proteins.csv",
          row.names=FALSE)
unmatched.duplicate.protein.product.annotations <- unique(unmatched.duplicate.proteins$product)

## free up memory by deallocating all.proteins,
rm(all.proteins)
## and running garbage collection.
gc()

########################################################################
cds.counts <- read.csv("../results/protein_db_CDS_counts.csv")

protein.db.metadata <- episome.database %>%
    inner_join(gbk.annotation) %>%
    inner_join(cds.counts)

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

###########################################################################
## Analysis for Figure 1.
## Panel A is a schematic of the analysis pipeline.

## See Wikipedia reference:
## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

## Make Z-distributed confidence intervals for the fraction of isolates with
## duplicated ARGs (panel B),
## the fraction of isolates with single-copy ARGs (panel C),
## the fraction of isolates with duplicated genes (panel D).

## Count data for isolates with duplicated ARGs
## goes into Supplementary Table S1.

calc.isolate.confints <- function(df) {
    df %>%
        ## use the normal approximation for binomial proportion conf.ints
        mutate(se = sqrt(p*(1-p)/total_isolates)) %>%
        ## and the Rule of Three to handle zeros.
        ## See Wikipedia reference:
        ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
        mutate(Left = ifelse(p > 0, p - 1.96*se, 0)) %>%
        mutate(Right = ifelse(p > 0, p + 1.96*se, 3/total_isolates)) %>%
        ## Sort every table by the total number of isolates.
        arrange(desc(total_isolates))
}


make.TableS1 <- function(gbk.annotation, duplicate.proteins) {

    ## count the number of isolates with duplicated ARGs in each category.
    ARG.category.counts <- duplicate.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_duplicated_ARGs = n)
    
    ## join columns to make Table S1.
    TableS1 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(ARG.category.counts) %>%
        mutate(isolates_with_duplicated_ARGs =
                   replace_na(isolates_with_duplicated_ARGs,0)) %>%
        mutate(p = isolates_with_duplicated_ARGs/total_isolates) %>%
        calc.isolate.confints()
    
    return(TableS1)
}


make.confint.figure.panel <- function(Table, order.by.total.isolates, title) {    
    Fig1.panel <- Table %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total.isolates))) %>%
        ggplot(aes(y = Annotation, x = p)) +     
        geom_point(size=2) +
        ylab("Annotation") +
        xlab("Proportion of Isolates") +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=1, size=0.5)
    return(Fig1.panel)
}


## Figure 1B: normal-approximation confidence intervals for the percentage
## of isolates with duplicated ARGs.
TableS1 <- make.TableS1(gbk.annotation, duplicate.proteins)
## write Supplementary Table S1 to file.
write.csv(x=TableS1, file="../results/TableS1.csv")
Fig1B <- make.confint.figure.panel(TableS1, order.by.total.isolates, "Duplicated ARGs")

######################
## Table S2. Control: does the distribution of ARG singletons
## (i.e. genes that have NOT duplicated) follow the distribution
## of sampled isolates?

## No categories are enriched with ARG singletons,
## as most isolates have a gene that matches an antibiotic keyword.
## Animal-host isolates are depleted (perhaps due to aphid bacteria isolates?)

make.TableS2 <- function(gbk.annotation, singleton.proteins) {

## count the number of isolates with singleton AR genes in each category.
    ARG.category.counts <- singleton.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_singleton_ARGs = n()) %>%
        arrange(desc(isolates_with_singleton_ARGs))
    gc() ## free memory.
    
   
    ## join columns to make Table S2.
    TableS2 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(ARG.category.counts) %>%
        mutate(isolates_with_singleton_ARGs =
                   replace_na(isolates_with_singleton_ARGs,0)) %>%
        mutate(p = isolates_with_singleton_ARGs/total_isolates) %>%
        calc.isolate.confints()
    return(TableS2)
}

## This data frame will be used for Figure 1C.
TableS2 <- make.TableS2(gbk.annotation, singleton.proteins)
## write TableS2 to file.
write.csv(x=TableS2, file="../results/TableS2.csv")
Fig1C <- make.confint.figure.panel(TableS2, order.by.total.isolates, "Single-copy ARGs")
gc() ## free memory after dealing with singleton data.

###############################################################################
## Table S3. Control: does the number of isolates with duplicate genes
## follow the sampling distribution of isolates?

## Most follow the expected distribution.
## however, isolates from animal-hosts are signficantly depleted
## in duplicate genes: FDR-corrected p = 0.0000314
## while isolates from anthropogenic environments are weakly enriched
## in multi-copy genes: FDR-corrected p = 0.0212.

make.TableS3 <- function(gbk.annotation, duplicate.proteins) {
    ## count the number of isolates with duplicated genes in each category.
    category.counts <- duplicate.proteins %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_duplicated_genes = n()) %>%
        arrange(desc(isolates_with_duplicated_genes))
    
    ## join columns to make Table S3.
    TableS3 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_duplicated_genes =
                   replace_na(isolates_with_duplicated_genes, 0)) %>%
        mutate(p = isolates_with_duplicated_genes/total_isolates) %>%
        calc.isolate.confints()
    return(TableS3)
}


TableS3 <- make.TableS3(gbk.annotation, duplicate.proteins)
## write TableS3 to file.
write.csv(x=TableS3, file="../results/TableS3.csv")

Fig1D <- make.confint.figure.panel(TableS3, order.by.total.isolates, "All Duplicated Genes")

Fig1BCD <- plot_grid(Fig1B, Fig1C, Fig1D, labels=c('B','C','D'),ncol=1)
ggsave("../results/Fig1BCD.pdf", width = 5)
######################################################################
## Control for Taxonomy (simpler control than for phylogeny)

## let's look at the taxonomic distribution of strains with duplicated ARGs.
duplicated.ARG.seq.genera.summary <- duplicate.proteins %>%
    tibble() %>% filter(str_detect(product,antibiotic.keywords)) %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%
    group_by(Genus) %>%
    summarize(duplicated.ARG.count = n()) %>%
    arrange(desc(duplicated.ARG.count))

duplicated.genera.seq.summary <- duplicate.proteins %>%
    tibble() %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%
    group_by(Genus) %>%
    summarize(duplicated.seq.count = n()) %>%
    arrange(desc(duplicated.seq.count))

all.genera.isolate.summary <- gbk.annotation %>%
    filter(Annotation != "Unannotated") %>%
    filter(Annotation != "blank") %>%
    tibble() %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%
    group_by(Genus) %>%
    summarize(genome.count = n()) %>%
    arrange(desc(genome.count))

duplicated.genera.isolate.summary <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(Annotation_Accession, Organism, Strain, Annotation) %>%
    distinct() %>%
    tibble() %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%
    group_by(Genus) %>%
    summarize(duplicated.ARG.genome.count = n()) %>%
    arrange(desc(duplicated.ARG.genome.count))

## Let's mark the most commonly sampled genera, and recalculate a version
## of TableS1, and the associated figure.
top.ARG.genera <- c("Klebsiella", "Escherichia",
                    "Acinetobacter", "Salmonella",
                    "Enterobacter", "Staphylococcus",
                    "Citrobacter", "Pseudomonas",
                    "Enterococcus","Bacillus")

genera.isolate.comparison.df <- full_join(
    duplicated.genera.isolate.summary,
    all.genera.isolate.summary) %>%
    ## turn NAs to zeros.
    replace(is.na(.), 0) %>%
    mutate(percent.genomes.with.dup.ARGs = duplicated.ARG.genome.count/genome.count) %>%
    mutate(in.top.ARG.genera = ifelse(Genus %in% top.ARG.genera, TRUE, FALSE)) %>%
    arrange(desc(percent.genomes.with.dup.ARGs))


## Hmmm... sampling biases could be problematic,
## since antibiotic-resistant bacteria are more likely to be sequenced.
S1FigA <- genera.isolate.comparison.df %>%
    ggplot(aes(x=sqrt(genome.count),
               y = sqrt(duplicated.ARG.genome.count),
               label = Genus,
               color = in.top.ARG.genera)) +
    theme_classic() + geom_jitter() + geom_text_repel() +
    scale_color_manual(values=c("black", "red")) +
    guides(color=FALSE)

## 4,559 isolates are in the top ARG genera.
top.ARG.genera.isolates <- gbk.annotation %>%
    filter(Genus %in% top.ARG.genera)

## 3,780 isolates are in the remaining genera.
filtered.gbk.annotation <- gbk.annotation %>%
    filter(!(Genus %in% top.ARG.genera))

filtered.duplicate.proteins <- duplicate.proteins %>%
    filter(!(Genus %in% top.ARG.genera))

filtered.TableS1 <- make.TableS1(filtered.gbk.annotation, filtered.duplicate.proteins)
S1FigB <- make.confint.figure.panel(filtered.TableS1, order.by.total.isolates, "Duplicated ARGs after filtering top genera")

S1Fig <- plot_grid(S1FigA, S1FigB, labels = c("A", "B"), nrow = 2)
ggsave("../results/S1Fig.pdf", S1Fig)

###########################################################################
## Figure S2: Distribution of proportions of duplicated ARGs per genome:
## For each Annotation category:
## panel A: the % of genes in each genome that are duplicated ARGs.
## panel B: the % of genes in each genome that are duplicated MGEs.
## panel C: the % of genes in each genome that are duplicated.

## calculate the number of duplicate genes in each genome.
duplicated.genes.per.genome <- duplicate.proteins %>%
    group_by(Annotation_Accession) %>%
    summarize(gene.duplicates.per.genome = sum(count))
    
duplicated.MGEs.per.genome <- duplicate.proteins %>%
    filter(str_detect(.$product, IS.keywords)) %>%
    group_by(Annotation_Accession) %>%
    summarize(MGE.duplicates.per.genome = sum(count))

duplicated.ARGs.per.genome <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Annotation_Accession) %>%
    summarize(ARG.duplicates.per.genome = sum(count))
    
boxplot.df <- protein.db.metadata %>%
    ## sum the CDS in all chromosomes and plasmids per genome.
    group_by(Annotation_Accession, Annotation) %>%
    summarize(total_genes = sum(CDS_count)) %>%
    ## calculate percentage for all duplicated.genes.
    left_join(duplicated.genes.per.genome) %>%
    mutate(gene.duplicates.per.genome =
               replace_na(gene.duplicates.per.genome, 0)) %>%
    mutate(proportion.of.duplicated.genes = gene.duplicates.per.genome/total_genes) %>%
    ## calculate percentage for duplicated MGE-associated genes.
    left_join(duplicated.MGEs.per.genome) %>%
    mutate(MGE.duplicates.per.genome =
               replace_na(MGE.duplicates.per.genome, 0)) %>%
    mutate(proportion.of.duplicated.MGEs = MGE.duplicates.per.genome/total_genes) %>%
    ## calculate percentage for duplicated ARGs.
    left_join(duplicated.ARGs.per.genome) %>%
    mutate(ARG.duplicates.per.genome =
               replace_na(ARG.duplicates.per.genome, 0)) %>%
    mutate(proportion.of.duplicated.ARGs = ARG.duplicates.per.genome/total_genes)
   
S2FigA <- boxplot.df %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates))) %>%
    ggplot(aes(y = Annotation, x = proportion.of.duplicated.ARGs)) +
    geom_jitter(size=0.2) +
    ylab("Annotation") +
    xlab("Proportion of duplicated ARGs per genome") +
    theme_classic() +
    ggtitle("Duplicated ARGs")

S2FigB <- boxplot.df %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates))) %>%
    ggplot(aes(y = Annotation, x = proportion.of.duplicated.MGEs)) +
    geom_jitter(size=0.2) +
    ylab("Annotation") +
    xlab("Proportion of duplicated MGE genes per genome") +
    theme_classic() +
    ggtitle("Duplicated MGE genes")

S2FigC <- boxplot.df %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates))) %>%
    ggplot(aes(y = Annotation, x = proportion.of.duplicated.genes)) +
    geom_jitter(size=0.2) +
    ylab("Annotation") +
    xlab("Proportion of duplicated genes per genome") +
    theme_classic() +
    ggtitle("Duplicated genes")

S2Fig <- plot_grid(S2FigA, S2FigB, S2FigC, labels = c('A','B','C'), nrow = 3)
ggsave("../results/S2Fig.pdf", S2Fig)

######################################################################
## Figure 2: Visualization of ARGs on plasmids and chromosomes.

categorize.as.MGE.ARG.or.other <- function(product) {
    if (is.na(product))
        return("Other function")
    else if (str_detect(product, antibiotic.keywords))
        return("ARG")
    else if (str_detect(product, EFTu.keywords))
        return("EF-Tu")
    else if (str_detect(product, "\\bribosomal protein\\b"))
        return("ribosomal protein")
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
               levels = rev(order.by.total.isolates)))

Fig2B.data <- singleton.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation, Category) %>%
    summarize(Plasmid = sum(plasmid_count), Chromosome = sum(chromosome_count)) %>%
    pivot_longer(cols = c("Plasmid", "Chromosome"),
                 names_to = "Episome",
                 values_to = "Count") %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

Fig2A <- ggplot(Fig2A.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    facet_wrap(.~Episome) +
    theme_classic() +
    ggtitle("Distribution of duplicated genes") +
    xlab("Proportion of genes") +
    theme(legend.position="bottom")

Fig2legend <- get_legend(Fig2A)
Fig2A <- Fig2A + guides(fill = FALSE)

Fig2B <- ggplot(Fig2B.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    facet_wrap(.~Episome) +
    theme_classic() +
    ggtitle("Distribution of single-copy genes") +
    xlab("Proportion of genes") +
    guides(fill = FALSE)

stackedbar.Fig2A <- ggplot(Fig2A.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity") +
    facet_wrap(.~Episome, scales = "free") +
    theme_classic() +
    ggtitle("Distribution of duplicated genes") +
    scale_x_continuous(labels=fancy_scientific)

stackedbar.Fig2B <- ggplot(Fig2B.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity") +
    facet_wrap(.~Episome, scales = "free") +
    theme_classic() +
    ggtitle("Distribution of single-copy genes") +
    guides(fill = FALSE) +
    scale_x_continuous(labels=fancy_scientific)


Fig2 <- plot_grid(Fig2A, Fig2B, Fig2legend, labels = c("A", "B"),
                  ncol = 1, rel_heights = c(1,1,0.25))
ggsave("../results/Fig2.pdf", Fig2, height = 7, width = 6)

## This visualization is also useful.
stackedbar.Fig2 <- plot_grid(stackedbar.Fig2A,
                             stackedbar.Fig2B, labels = c("A", "B"), ncol = 1)
ggsave("../results/stackedbar-Fig2.pdf", stackedbar.Fig2)

#########################
## S3 Figure.
## Analysis of duplicate pairs found just on chromosome, just on plasmid, or
## on both chromosomes and plasmids.

## let's look at cases of identical sequences on chromosomes and plasmids.

both.chr.and.plasmid.cases <- duplicate.proteins %>%
    filter(chromosome_count >= 1 & plasmid_count >= 1) %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    arrange(desc(count)) %>%
    tibble()

just.chromosome.cases <- duplicate.proteins %>%
    filter(chromosome_count >= 1 & plasmid_count == 0) %>%
    arrange(desc(count)) %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    tibble()

just.plasmid.cases <- duplicate.proteins %>%
    filter(chromosome_count == 0 & plasmid_count >= 1) %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    arrange(desc(count)) %>%
    tibble()

both.chr.and.plasmid.summary <- both.chr.and.plasmid.cases %>%
    group_by(Annotation, Category) %>%
    summarize(Count = sum(count)) %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

just.chromosome.summary <- just.chromosome.cases %>%
    group_by(Annotation, Category) %>%
    summarize(Count = sum(count)) %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

just.plasmid.summary <- just.plasmid.cases %>%
    group_by(Annotation, Category) %>%
    summarize(Count = sum(count)) %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

S3FigA <- ggplot(both.chr.and.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("Both chromosome and plasmid") +
    theme(legend.position="bottom")

S3Fig.legend <- get_legend(S3FigA)
S3FigA <- S3FigA + guides(fill = FALSE)

S3FigB <- ggplot(just.chromosome.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("chromosome only") +
    guides(fill = FALSE)

S3FigC <- ggplot(just.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("plasmid only") +
    guides(fill = FALSE)

S3Fig <- plot_grid(NULL, S3FigA, S3FigB, S3FigC, S3Fig.legend, ncol = 1,
                   labels = c("Genomic distribution of duplicated genes", "A","B","C"),
                   rel_heights=c(0.2,1,1,1,0.25))
ggsave("../results/S3Fig.pdf", S3Fig)

##################################################################################

## Let's take a closer look at duplicated EF-Tu sequences,
## and examine ribosomal proteins too.

Fig2A.summary <- Fig2A.data %>% group_by(Episome, Category) %>%
    summarize(count = sum(Count))

## Remember: singleton is based on 100% sequence identity.
## so a pair of duplicate gene with a single amino acid difference
## is counted as a pair of singletons.

Fig2B.summary <- Fig2B.data %>% group_by(Episome, Category) %>%
    summarize(count = sum(Count))

duplicated.EF.Tu <- duplicate.proteins %>%
    tibble() %>% filter(str_detect(product,EFTu.keywords))
singleton.EF.Tu <- singleton.proteins %>%
    tibble() %>% filter(str_detect(product,EFTu.keywords))

## find average number of EF-Tu sequences per genome. 1.52 per genome.
num.genomes <- nrow(gbk.annotation)
(sum(filter(Fig2A.summary,Category=="EF-Tu")$count) +
 sum(filter(Fig2B.summary,Category=="EF-Tu")$count))/num.genomes

duplicated.EF.Tu.genera.summary <- duplicated.EF.Tu %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%
    group_by(Genus, Annotation) %>%
    summarize(duplicated.EF.Tu.count = n()) %>%
    arrange(desc(duplicated.EF.Tu.count))

singleton.EF.Tu.genera.summary <- singleton.EF.Tu %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%    
    group_by(Genus, Annotation) %>%
    summarize(singleton.EF.Tu.count = n()) %>%
    arrange(desc(singleton.EF.Tu.count))

EF.Tu.genera.summary <- full_join(
    duplicated.EF.Tu.genera.summary,
    singleton.EF.Tu.genera.summary) %>%
    ## turn NAs to zeros.
    replace(is.na(.), 0)

duplicate.ribosomal.proteins <- duplicate.proteins %>%
    tibble() %>% filter(str_detect(product,"\\bribosomal protein\\b"))

duplicate.ribosomal.protein.genera.summary <- duplicate.ribosomal.proteins %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%
    group_by(Genus, Annotation) %>%
    summarize(duplicated.ribosomal.proteins.count = n()) %>%
    arrange(desc(duplicated.ribosomal.proteins.count))

duplicate.ribosomal.protein.genera.isolate.summary <- duplicate.ribosomal.proteins %>%
    ## next two lines is to count isolates rather than genes
    select(Annotation_Accession, Organism, Strain, Annotation) %>%
    distinct() %>%
    tibble() %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%
    group_by(Genus) %>%
    summarize(duplicated.ribosomal.proteins.isolates = n()) %>%
    arrange(desc(duplicated.ribosomal.proteins.isolates))

######################################################################
## Table S4. Show number of duplicated genes on chromosomes, and number of
## duplicate genes on plasmids, for each category, for duplicated genes
## and duplicated AR genes. This table of raw data goes into the text. Then,
## sum over all categories for a 2x2 contingency table and report the result of a
## Fisher's exact test for asssociation between duplicated AR genes and plasmids.

make.TableS4 <- function(duplicate.proteins) {
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

TableS4 <- make.TableS4(duplicate.proteins)
## write Table S4 to file.
write.csv(x=TableS4,file="../results/TableS4.csv")

################
## Analysis of Table S4: Duplicate ARGs are associated with plasmids.

plasmid.chromosome.duplicate.ARG.contingency.test <- function(TableS4) {
    ## get values for Fisher's exact test.
    total.chr.AR.duplicates <- sum(TableS4$chromosomal_duplicate_ARGs)
    total.plasmid.AR.duplicates <- sum(TableS4$plasmid_duplicate_ARGs)

    total.chr.duplicates <- sum(TableS4$chromosomal_duplicate_genes)
    total.plasmid.duplicates <- sum(TableS4$plasmid_duplicate_genes)

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

plasmid.chromosome.duplicate.ARG.contingency.test(TableS4)

####################
## Table S5: look at distribution of single-copy ARGs on
## chromosomes and plasmids.

## This control shows that single-copy ARGs are highly enriched on plasmids,
## based on a comparison with the distribution of singleton genes overall.
## Therefore AR genes are generally associated with plasmids, regardless of
## status of being a duplication or not.

## This does NOT invalidate the main result of this analysis, that duplicate AR
## genes are more enriched on plasmids in comparison to the distribution of
## duplicate genes overall.

## Notably, the majority of duplicate ARGs are on plasmids, while the
## majority of singleton ARGs are on chromosomes.

make.TableS5 <- function(singleton.proteins) {
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
    
    TableS5 <- singleton.chromosome.genes.count %>%
        left_join(singleton.plasmid.genes.count) %>%
        left_join(singleton.chromosome.ARGs.count) %>%
        mutate(chromosomal_singleton_ARGs =
                    replace_na(chromosomal_singleton_ARGs, 0)) %>%
        left_join(singleton.plasmid.ARGs.count) %>%
        mutate(plasmid_singleton_ARGs = replace_na(plasmid_singleton_ARGs, 0)) %>%
        arrange(desc(plasmid_singleton_ARGs))
}

TableS5 <- make.TableS5(singleton.proteins)
## write Table S5 to file.
write.csv(x=TableS5,file="../results/TableS5.csv")
gc()


plasmid.chromosome.singleton.ARG.contingency.test <- function(TableS5) {
    ## get values for Fisher's exact test.
    total.chr.AR.singletons <- sum(TableS5$chromosomal_singleton_ARGs)
    total.plasmid.AR.singletons <- sum(TableS5$plasmid_singleton_ARGs)

    total.chr.singletons <- sum(TableS5$chromosomal_singleton_genes)
    total.plasmid.singletons <- sum(TableS5$plasmid_singleton_genes)

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

plasmid.chromosome.singleton.ARG.contingency.test(TableS5)

################################################################################
## Use the data in Tables S4 and S5 to make Figure 3.
## The point of this figure is to show that the distribution of
## duplicated ARGs is not predicted by the distribution of single-copy ARGs
## in the ecological categories.

make.Fig3.df <- function(TableS1, TableS4, TableS5) {
    order.by.total_isolates <- TableS1$Annotation
    
    Fig3.df <- full_join(TableS4, TableS5) %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = order.by.total_isolates)) %>%
        mutate(total_duplicate_genes =
                   plasmid_duplicate_genes + chromosomal_duplicate_genes) %>%
        mutate(total_singleton_genes =
                   plasmid_singleton_genes + chromosomal_singleton_genes) %>%
        mutate(total_duplicate_ARGs =
                   plasmid_duplicate_ARGs + chromosomal_duplicate_ARGs) %>%
        mutate(total_singleton_ARGs = plasmid_singleton_ARGs + chromosomal_singleton_ARGs) %>%
        mutate(total_chromosomal_genes = chromosomal_duplicate_genes + chromosomal_singleton_genes) %>%
        mutate(total_plasmid_genes = plasmid_duplicate_genes + plasmid_singleton_genes) %>%
        mutate(total_genes = total_duplicate_genes + total_singleton_genes)
        
    return(Fig3.df)
}


Fig3.df <- make.Fig3.df(TableS1, TableS4, TableS5)

## Figure 3.
## Plot point estimates for the fraction of chromosomal genes that are
## the fraction of genes that are duplicated ARGs (panel A),
## the fraction of genes that are single-copy ARGs (panel B).
## duplicated ARGs (panel C),
## the fraction of plasmid genes that are duplicated ARGs (panel D),
## the fraction of chromosomal genes that are single-copy ARGs (panel E),
## the fraction of plasmid genes that are single-copy ARGs (panel F),

## Fig3A
Fig3A.df <- Fig3.df %>%
    mutate(p = total_duplicate_ARGs/(total_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_genes)) %>%
    ## and the Rule of Three to handle zeros.
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = ifelse(p > 0, p - 1.96*se, 0)) %>%
    mutate(Right = ifelse(p > 0, p + 1.96*se, 3/total_genes)) %>%
    select(Annotation, total_duplicate_ARGs, total_genes,
           p, Left, Right)

Fig3B.df <- Fig3.df %>%
    mutate(p = total_singleton_ARGs/(total_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_genes)) %>%
    ## and the Rule of Three to handle zeros.
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = ifelse(p > 0, p - 1.96*se, 0)) %>%
    mutate(Right = ifelse(p > 0, p + 1.96*se, 3/total_genes)) %>%
    select(Annotation, total_singleton_ARGs, total_genes,
           p, Left, Right)


## Fig3C: the fraction of chromosomal genes that are duplicated ARGs.
Fig3C.df <- Fig3.df %>%
    mutate(p = chromosomal_duplicate_ARGs/(total_chromosomal_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_chromosomal_genes)) %>%
    ## and the Rule of Three to handle zeros.
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = ifelse(p > 0, p - 1.96*se, 0)) %>%
    mutate(Right = ifelse(p > 0, p + 1.96*se, 3/total_chromosomal_genes)) %>%
    select(Annotation, chromosomal_duplicate_ARGs, total_chromosomal_genes,
           p, Left, Right)


Fig3D.df <- Fig3.df %>%
    mutate(p = plasmid_duplicate_ARGs/(total_plasmid_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_plasmid_genes)) %>%
    ## and the Rule of Three to handle zeros.
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = ifelse(p > 0, p - 1.96*se, 0)) %>%
    mutate(Right = ifelse(p > 0, p + 1.96*se, 3/total_plasmid_genes)) %>%
    select(Annotation, plasmid_duplicate_ARGs, total_plasmid_genes,
           p, Left, Right)

Fig3E.df <- Fig3.df %>%
    mutate(p = chromosomal_singleton_ARGs/(total_chromosomal_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_chromosomal_genes)) %>%
    ## and the Rule of Three to handle zeros.
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = ifelse(p > 0, p - 1.96*se, 0)) %>%
    mutate(Right = ifelse(p > 0, p + 1.96*se, 3/total_chromosomal_genes)) %>%
    select(Annotation, chromosomal_singleton_ARGs, total_chromosomal_genes,
           p, Left, Right)


Fig3F.df <- Fig3.df %>%
    mutate(p = plasmid_singleton_ARGs/(total_plasmid_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_plasmid_genes)) %>%
    ## and the Rule of Three to handle zeros.
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = ifelse(p > 0, p - 1.96*se, 0)) %>%
    mutate(Right = ifelse(p > 0, p + 1.96*se, 3/total_plasmid_genes)) %>%
    select(Annotation, plasmid_singleton_ARGs, total_plasmid_genes,
           p, Left, Right)



make.Fig3.panel <- function(Table, order.by.total.isolates, title, xlabel) {
    Fig3.panel <- Table %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total.isolates))) %>%
        ggplot(aes(y = Annotation, x = p)) +     
        geom_point(size=2) +
        ylab("Annotation") +
        xlab(xlabel) +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=1, size=0.5)
    return(Fig3.panel)
}


Fig3A <- make.Fig3.panel(Fig3A.df, order.by.total.isolates,
                         "Duplicated ARGs",
                         "Proportion of all genes")
Fig3B <- make.Fig3.panel(Fig3B.df, order.by.total.isolates,
                         "Single-copy ARGs",
                         "Proportion of all genes")
Fig3C <- make.Fig3.panel(Fig3C.df, order.by.total.isolates,
                         "Duplicated ARGs on the chromosome",
                         "Proportion of chromosomal genes")
Fig3D <- make.Fig3.panel(Fig3D.df, order.by.total.isolates,
                         "Duplicated ARGs on plasmids",
                         "Proportion of plasmid genes")
Fig3E <- make.Fig3.panel(Fig3E.df, order.by.total.isolates,
                         "Single-copy ARGs on the chromosome",
                         "Proportion of chromosomal genes")
Fig3F <- make.Fig3.panel(Fig3F.df, order.by.total.isolates,
                         "Single-copy ARGs on plasmids",
                         "Proportion of plasmid genes")

Fig3 <- plot_grid(Fig3A, Fig3B, Fig3C, Fig3D, Fig3E, Fig3F,
                  labels = c('A','B','C','D','E','F'), ncol=2)
ggsave("../results/Fig3.pdf", Fig3)

################################################################################
## Figure 4: 
## The observed ecological distribution of duplicate genes is driven by either
## selection, HGT, or associations with MGEs.

## In the absence of selection, HGT, or association with MGEs,
## the distribution of non-MGE duplicated genes should be a random sample of
## non-MGE singletons.

## Null hypothesis: ratio of duplicated ARGs to all duplicated genes
## should be proportional to the number of singleton ARGs out of all
## singleton genes.

## Deviation from the null hypothesis indicates selection, HGT, or linkage with
## MGEs, thus causing enrichment.

## Panel A is a schematic figure in Illustrator to show the rationale.

duplicated.ARGs.per.category <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Annotation) %>%
    summarize(ARG.duplicates = sum(count))

duplicated.genes.per.category <- duplicate.proteins %>%
    group_by(Annotation) %>%
    summarize(gene.duplicates = sum(count))

duplicated.MGEs.per.category <- duplicate.proteins %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Annotation) %>%
    summarize(MGE.duplicates = sum(count))

singleton.ARGs.per.category <- singleton.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Annotation) %>%
    summarize(singleton.ARGs = sum(count))

singleton.genes.per.category <- singleton.proteins %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(Annotation) %>%
    summarize(singleton.genes = sum(count))

singleton.MGEs.per.category <- singleton.proteins %>%
    filter(str_detect(.$product,IS.keywords)) %>%
    group_by(Annotation) %>%
    summarize(singleton.MGEs = sum(count))

selection.test.df <- duplicated.ARGs.per.category %>%
    full_join(duplicated.genes.per.category) %>%
    full_join(duplicated.MGEs.per.category) %>%
    full_join(singleton.ARGs.per.category) %>%
    full_join(singleton.genes.per.category) %>%
    full_join(singleton.MGEs.per.category) %>%
    ## turn NAs to zeros.
    replace(is.na(.), 0) %>%
    mutate(p1 = ARG.duplicates / gene.duplicates) %>%
    mutate(p2 = MGE.duplicates / gene.duplicates) %>%
    mutate(q1 = singleton.ARGs / singleton.genes) %>%
    mutate(q2 = singleton.MGEs / singleton.genes) %>%
    mutate(ARG.dup.singleton.ratio = p1/q1) %>%
    mutate(MGE.dup.singleton.ratio = p2/q2) %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

ARG.selection.plot <- selection.test.df %>%
    ggplot(aes(y = Annotation, x = ARG.dup.singleton.ratio)) +
    geom_point() + theme_classic() +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    xlim(0,4) +
    xlab("proportion of ARGs among duplicate genes / proportion of ARGs among single-copy genes")

MGE.selection.plot <- selection.test.df %>%
    ggplot(aes(y = Annotation, x = MGE.dup.singleton.ratio)) +
    geom_point() + theme_classic() +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    xlim(0,4) +
    xlab("proportion of MGE genes among duplicate genes / proportion of MGE genes among single-copy genes")

Fig4BC <- plot_grid(ARG.selection.plot, MGE.selection.plot, labels = c('B','C'), nrow = 2)
ggsave("../results/Fig4BC.pdf", Fig4BC, width=9.5, height=5)

################################################################################
## Figure S4. A deterministic ODE model demonstrates that selection can
## drive the evolution of duplicated ARGs on plasmids.

## The panels of this figure are generated in my Pluto notebook:
## duplication-linear-ODE-model.jl.
################################################################################
## Figure 5: examples that indicate generality of our method.
## let's examine some other functions that we expect to be enriched in some, but
## not all ecological annotations.

## generic version of make.TableS1, for examining classes of genes other than
## antibiotic resistance genes.

make.IsolateEnrichmentTable <- function(gbk.annotation, duplicate.genes, keywords) {
    ## count the number of isolates with duplicated genes of interest in each category.
    category.counts <- duplicate.proteins %>%
        filter(str_detect(.$product, keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_duplicated_function = n)
    
    ## join columns to make the Table.
    Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_duplicated_function =
                   replace_na(isolates_with_duplicated_function,0)) %>%
        mutate(p = isolates_with_duplicated_function/total_isolates) %>%
        calc.isolate.confints()
    return(Table)
}


photosynthesis.table <- make.IsolateEnrichmentTable(gbk.annotation,
                                                    duplicate.genes,
                                                    "photosystem")
write.csv(x=photosynthesis.table, file="../results/TableS6.csv")

N2.fixation.table <- make.IsolateEnrichmentTable(gbk.annotation,
                                                 duplicate.genes,
                                                 "nitrogenase")
write.csv(x=N2.fixation.table, file="../results/TableS7.csv")

toxic.metal.table <- make.IsolateEnrichmentTable(gbk.annotation,
                                                 duplicate.genes,
                                                 "mercury|cadmium|arsen")
write.csv(x=toxic.metal.table, file="../results/TableS8.csv")

heme.table <- make.IsolateEnrichmentTable(gbk.annotation,
                                          duplicate.genes,
                                          "heme")
write.csv(x=heme.table, file="../results/TableS9.csv")


Fig5A <- make.confint.figure.panel(
    photosynthesis.table, order.by.total.isolates,
    "Photosynthesis") +
    ylab("Isolates with duplicated photosystem genes")

Fig5B <- make.confint.figure.panel(
    N2.fixation.table, order.by.total.isolates,
    "Nitrogen fixation") +
    ylab("Isolates with duplicated nitrogenase genes")

Fig5C <- make.confint.figure.panel(
    toxic.metal.table, order.by.total.isolates,
    "Toxic-metal resistance") +
    ylab("Isolates with duplicated toxic-metal resistance genes")

Fig5D <- make.confint.figure.panel(
    heme.table, order.by.total.isolates,
    "Heme degradation") +
    ylab("Isolates with duplicated heme degradation genes")

Fig5 <- plot_grid(Fig5A, Fig5B, Fig5C, Fig5D,
                  labels = c("A","B","C","D"),
                  nrow = 2)

ggsave(Fig5, file = "../results/Fig5.pdf", width = 8.5, height = 8.5)
##########################################################################

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
## proteins in each environment, after removing MGEs, EF-Tu, and unknown,
## hypothetical, or uncharacterized proteins (since many different protein families
## can be described this way).

## This stuff will go into the Supplementary Material.

## this function filters proteins by manual annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
.make.annotation.freq.table <- function(data.df, manual.annot.string, remove.MGEs = TRUE) {

    df <- data.df %>%
        filter(Annotation == manual.annot.string) %>%
        filter(!str_detect(.$product, unknown.protein.keywords)) %>%
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
    slice_max(tf_idf, n = 5) %>%
    ungroup()

dup.prot.annotation.tf_idf.plot <- top.dup.prot.annotation.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~Annotation, ncol = 2, scales = "free") +
  labs(x = "tf-idf", y = NULL)

ggsave("../results/duplicate-protein-annotation-TF-IDF.pdf",
       dup.prot.annotation.tf_idf.plot,
       height=21,width=21)

#### Now repeat this analysis, but with singleton proteins.
#### This is a really valuable comparison.
#### Do annotations of duplicate proteins carry more ecological information
#### than the annotations of singleton proteins?

big.sing.prot.annotation.freq.table <- map_dfr(
    unique(singleton.proteins$Annotation),
    .f = make.sing.annotation.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all annotations in the table.
    mutate(total.annotation.count = sum(annotation.count))

sing.prot.annotation.tf_idf <- big.sing.prot.annotation.freq.table %>%
  bind_tf_idf(product, Annotation, annotation.count) %>%
  arrange(desc(tf_idf))

top.sing.prot.annotation.tf_idf <- sing.prot.annotation.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 100) %>%
    ungroup()

sing.prot.annotation.tf_idf.plot <- top.sing.prot.annotation.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~Annotation, ncol = 2, scales = "free") +
  labs(x = "tf-idf", y = NULL)

ggsave("../results/singleton-protein-annotation-TF-IDF.pdf",
       sing.prot.annotation.tf_idf.plot,
       height=21,width=21)

#################

ranked.dup.prot.annotation.tf_idf <- dup.prot.annotation.tf_idf %>%
    group_by(Annotation) %>%
    mutate(Rank = rank(desc(tf_idf))) %>%
    ungroup()

ranked.sing.prot.annotation.tf_idf <- sing.prot.annotation.tf_idf %>%
    group_by(Annotation) %>%
    mutate(Rank = rank(desc(tf_idf))) %>%
    ungroup()

## let's plot the tf_idf distributions per Annotation.
dup.prot.annotation.tf_idf.dist.plot <- ranked.dup.prot.annotation.tf_idf %>%
    ggplot(aes(x=Rank, y = tf_idf)) +
    geom_line() +
    facet_wrap(.~Annotation)

sing.prot.annotation.tf_idf.dist.plot <- ranked.sing.prot.annotation.tf_idf %>%
    ggplot(aes(x=Rank, y = tf_idf)) +
    geom_line() +
    facet_wrap(.~Annotation)

## now make eCDF plots.
dup.prot.annotation.tf_idf.cdf.plot <- ranked.dup.prot.annotation.tf_idf %>%
    ggplot(aes(x = tf_idf)) +
    stat_ecdf(geom = "step") +
    facet_wrap(.~Annotation)

sing.prot.annotation.tf_idf.cdf.plot <- ranked.sing.prot.annotation.tf_idf %>%
    ggplot(aes(x = tf_idf)) +
    stat_ecdf(geom = "step") +
    facet_wrap(.~Annotation)



retrieve.genomes.by.term <- function(tf_idf.df, annotated.proteins) {
    ## map each query term to genome accessions and annotations.
    retrieved.genomes <- tf_idf.df$product %>%
        map_dfr(.f = ~filter(annotated.proteins, product == .x)) %>%
        select(Annotation_Accession, host, isolation_source,
               Annotation, Organism, Strain) %>%
                distinct()
}


calc.precision <- function(retrieved.genomes, true.Annotation) {
    ## precision := (# of relevant & retrieved genomes)/(# of retrieved genomes).
    relevant.retrieved.genomes <- retrieved.genomes %>%
        filter(Annotation == true.Annotation)
    precision <- nrow(relevant.retrieved.genomes)/nrow(retrieved.genomes)
    return(precision)
}


calc.recall <- function(gbk.annotation, retrieved.genomes, true.Annotation) {
    ## recall := (# of retrieved & relevant genomes)/(# of relevant genomes).
    retrieved.relevant.genomes <- retrieved.genomes %>%
        filter(Annotation == true.Annotation)
    relevant.genomes <- filter(gbk.annotation, Annotation == true.Annotation)
    recall <- nrow(retrieved.relevant.genomes)/nrow(relevant.genomes)
    return(recall)
}


make.precision.recall.df <- function(gbk.annotation, annotated.proteins, tf_idf.df) {
    ## for each Annotation category, get the terms of interest,
    ## calculate precision for the category,
    ## calculate recall for the category,
    ## and return a dataframe of the results.
   
    precision.recall.helper.func <- function(true.Annotation) {
        queries <- filter(tf_idf.df, Annotation == true.Annotation)
        retrieved.genomes <- retrieve.genomes.by.term(queries, annotated.proteins)
        
        precision <- calc.precision(retrieved.genomes, true.Annotation)
        recall <- calc.recall(gbk.annotation, retrieved.genomes, true.Annotation)
        ret.df <- data.frame(Annotation = true.Annotation,
                             Precision = precision,
                             Recall = recall)
        return(ret.df)
    }


    all.annotations.vec <- unique(tf_idf.df$Annotation)

    precision.recall.df <- map_dfr(
        .x = all.annotations.vec,
        .f = precision.recall.helper.func)

    return(precision.recall.df)
}

dup.tf_idf_precision.recall.df <- make.precision.recall.df(
    gbk.annotation,
    duplicate.proteins,
    top.dup.prot.annotation.tf_idf) %>%
    mutate(class = "Multicopy")

sing.tf_idf_precision.recall.df <- make.precision.recall.df(
    gbk.annotation,
    singleton.proteins,
    top.sing.prot.annotation.tf_idf) %>%
    mutate(class = "Single-copy")

combined.tf_idf_precision.recall.df <- rbind(
    dup.tf_idf_precision.recall.df,
    sing.tf_idf_precision.recall.df) %>%
    mutate(F1_score = 2*(Precision*Recall)/(Precision+Recall))

precision.plot <- combined.tf_idf_precision.recall.df %>%
    ggplot(aes(x=Annotation, y = Precision, color = class)) +
    geom_point() +
    theme_classic() + ggtitle("Precision")

recall.plot <- combined.tf_idf_precision.recall.df %>%
    ggplot(aes(x=Annotation, y = Recall, color = class)) +
    geom_point() +
    theme_classic() + ggtitle("Recall")

F1_score.plot <- combined.tf_idf_precision.recall.df %>%
    ggplot(aes(x=Annotation, y = F1_score, color = class)) +
    geom_point() +
    theme_classic() + ggtitle("F1 score")

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
## Figure S5? the annotations of duplicated proteins are informative about
## ecology.

best.dup.prot.annotation.tf_idf <- dup.prot.annotation.tf_idf %>%
    filter(Annotation %in% c("Agriculture", "Anthropogenic-environment",
                             "Human-host")) %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 5) %>%
    ungroup()

S5Fig <- ggplot(best.dup.prot.annotation.tf_idf,
                aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
    geom_col(show.legend = FALSE) +
    labs(x = "tf-idf", y = NULL) +
    theme_classic() +
    xlim(0,0.04) +
    facet_wrap(.~Annotation, ncol=1, scales = "free_y") +
    ggtitle("Most informative multi-copy protein annotations")

ggsave("../results/S5Fig.pdf", S5Fig, width=8, height = 8)

