## Aim1-analysis.R by Rohan Maddamsetti.

## analyse the distribution of AR genes on chromosomes versus plasmids in
## fully-sequenced genomes and plasmids in the NCBI Nucleotide database.

## TODO: CAREFULLY re-number Supplementary Tables based on the manuscript.
## work backwards to avoid interchanging data structures between pieces of code.

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
## Analysis for Figure 1C.

## Statistical analysis for isolates with duplicated ARGs
## goes into Supplementary Table S1.

## return the first column for several tables.
## shows the number of isolates in each category.
make.isolate.totals.col <- function(gbk.annotation) {
    isolate.totals <- gbk.annotation %>%
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

## Use the total isolate order in TableS1 for organizing the pie chart figures.
order.by.total_isolates.vec <- TableS1$Annotation

###############################################################################
## Control: does the number of isolates with duplicate genes
## follow the sampling distribution of isolates?

## Most follow the expected distribution.
## however, isolates from animal-hosts are signficantly depleted
## in duplicate genes: FDR-corrected p = 0.0000314
## while isolates from anthropogenic environments are weakly enriched
## in multi-copy genes: FDR-corrected p = 0.0212.

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
               levels = rev(order.by.total_isolates.vec)))

Fig2B.data <- singleton.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation, Category) %>%
    summarize(Plasmid = sum(plasmid_count), Chromosome = sum(chromosome_count)) %>%
    pivot_longer(cols = c("Plasmid", "Chromosome"),
                 names_to = "Episome",
                 values_to = "Count") %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total_isolates.vec)))

Fig2A <- ggplot(Fig2A.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) + coord_polar() +
    facet_wrap(.~Episome) +
    theme_classic() +
    ggtitle("Distribution of multi-copy proteins") +
    xlab("Proportion of genes") +
    guides(fill = FALSE)

Fig2B <- ggplot(Fig2B.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) + coord_polar() +
    facet_wrap(.~Episome) +
    theme_classic() +
    ggtitle("Distribution of single-copy proteins") +
    xlab("Proportion of genes") +
    guides(fill = FALSE)

stackedbar.Fig2A <- ggplot(Fig2A.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity") +
    facet_wrap(.~Episome, scales = "free") +
    theme_classic() +
    ggtitle("Distribution of multi-copy proteins") +
    scale_x_continuous(labels=fancy_scientific)

stackedbar.Fig2B <- ggplot(Fig2B.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity") +
    facet_wrap(.~Episome, scales = "free") +
    theme_classic() +
    ggtitle("Distribution of single-copy proteins") +
    guides(fill = FALSE) +
    scale_x_continuous(labels=fancy_scientific)


Fig2AB <- plot_grid(Fig2A, Fig2B, labels = c("A", "B"), ncol = 1)
ggsave("../results/Fig2AB.pdf", Fig2AB, height = 9, width = 9)

## This visualization is also useful.
stackedbar.Fig2AB <- plot_grid(stackedbar.Fig2A,
                             stackedbar.Fig2B, labels = c("A", "B"), ncol = 1)
ggsave("../results/stackedbar-Fig2AB.pdf", stackedbar.Fig2AB, height = 9, width = 9)

#########################
## Figure 2C. TODO: FLESH OUT THIS SECTION.
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
               levels = rev(order.by.total_isolates.vec)))

just.chromosome.summary <- just.chromosome.cases %>%
    group_by(Annotation, Category) %>%
    summarize(Count = sum(count)) %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total_isolates.vec)))

just.plasmid.summary <- just.plasmid.cases %>%
    group_by(Annotation, Category) %>%
    summarize(Count = sum(count)) %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total_isolates.vec)))

Fig2C_1 <- ggplot(both.chr.and.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) + coord_polar() +
    theme_classic() +
    ggtitle("Both chromosome and plasmid") +
    guides(fill = FALSE) +
    guides(fill = FALSE) +
    theme(axis.line=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()#,
#          plot.background=element_blank()
          )


Fig2C_2 <- ggplot(just.chromosome.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) + coord_polar() +
    theme_classic() +
    ggtitle("chromosome only") +
    guides(fill = FALSE) +
    theme(axis.line=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()#,
#          plot.background=element_blank()
          )

Fig2C_3 <- ggplot(just.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) + coord_polar() +
    theme_classic() +
    ggtitle("plasmid only") +
    guides(fill = FALSE) +
    theme(axis.line=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()#,
          ##          plot.background=element_blank()
          )


Fig2C <- plot_grid(Fig2C_1, Fig2C_2, Fig2C_3, nrow = 1,
                   labels = c("C","Multi-copy proteins"))
ggsave("../results/Fig2C.pdf", Fig2C, height = 4, width = 9)

Fig2 <- plot_grid(Fig2AB, Fig2C, ncol = 1, rel_heights = c(9,4))
ggsave("../results/Fig2.pdf", Fig2, height = 13, width = 9)
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
                                    total_isolates = c(2000, 2000, 2000, 2000, 2000),
                                    y_pos = c(10^-0.55, 10^-0.30, 10^-0.05, 10^0.20, 10^0.45),
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
        annotate("text", x = 27, y = 2.6, label = "Isolates with multi-copy ARGs",
                 angle = 33.2, color = Fig1C.color.palette[1],size=3) +
        annotate("text", x = 27, y = 18.5, label = "Isolates with single-copy ARGs",
                 angle = 33.2, color = Fig1C.color.palette[2],size=3) +
        annotate("text", x = 27, y = 30, label = "Isolates with multi-copy genes",
                 angle = 33.2, color = "gray",size=3) +
        guides(size=FALSE) +
        ## Add the scale for the percent of the given gene class on plasmids.
        geom_point(data = plasmid_legend_df,
                   aes(y=y_pos,
                       size=plasmid_legend_percent * 0.5),
                   color="black",alpha=0.2) +
        annotate("text", x = 2000, y = 10^0.75, size = 3, label = "Percent on\nplasmids") +
        annotate("text", x = 2000, y = 10^0.45, size = 2, label = "100") +
        annotate("text", x = 2000, y = 10^0.20, size = 2, label = "75") +
        annotate("text", x = 2000, y = 10^-0.05, size = 2, label = "50") +
        annotate("text", x = 2000, y = 10^-0.30, size = 2, label = "25")
    
    return(Fig1C)
}

Fig1C <- make.Fig1C(Fig1C.df)
ggsave(Fig1C,file="../results/Fig1C.pdf",width=4.5,height=4.5)

################################################################################
make.big.gene.analysis.df <- function(duplicate.proteins, singleton.proteins) {
    ## We're going to make one dataframe with lots of columns for
    ## statistics and analyses. This is important for Figure 3.

    ## The particular statistics and supplementary tables are going to be
    ## calculated on relevant subsets of this dataframe (subsetting columns, not rows).

    ## The number of duplicated AR genes.
    duplicate.AR.genes.count <- duplicate.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        group_by(Annotation) %>%
        summarize(AR_duplicates = sum(count)) %>%
        arrange(desc(AR_duplicates))
    
    ## The number of duplicated MGE genes.
    duplicate.MGE.genes.count <- duplicate.proteins %>%
        filter(str_detect(.$product,IS.keywords)) %>%
        group_by(Annotation) %>%
        summarize(MGE_duplicates = sum(count)) %>%
        arrange(desc(MGE_duplicates))
    
    ## The number of duplicated genes in each category.
    duplicate.genes.count <- duplicate.proteins %>%
        group_by(Annotation) %>%
        summarize(duplicate_genes = sum(count)) %>%
        arrange(desc(duplicate_genes))
    
    ## The number of singleton genes in each category.
    singleton.genes.count <- singleton.proteins %>%
        group_by(Annotation) %>%
        summarize(singleton_genes = sum(count)) %>%
        arrange(desc(singleton_genes))
    
    ## The number of singleton AR genes.
    singleton.AR.genes.count <- singleton.proteins %>%
        filter(str_detect(.$product,antibiotic.keywords)) %>%
        group_by(Annotation) %>%
        summarize(AR_singletons = sum(count)) %>%
        arrange(desc(AR_singletons))

    ## The number of singleton MGE genes.
    singleton.MGE.genes.count <- singleton.proteins %>%
        filter(str_detect(.$product,IS.keywords)) %>%
        group_by(Annotation) %>%
        summarize(MGE_singletons = sum(count)) %>%
        arrange(desc(MGE_singletons))

    ## Sum up the totals for duplicate genes and singleton genes.
    ## This is the baseline column for statistical tests.
    total.genes.count <- duplicate.genes.count %>%
        full_join(singleton.genes.count) %>%
        mutate(total_genes = duplicate_genes + singleton_genes) %>%
        select(Annotation, total_genes)
    
    ## This data structure is important for later code.
    big.gene.analysis.df <- duplicate.AR.genes.count %>%
        full_join(duplicate.MGE.genes.count) %>%
        full_join(duplicate.genes.count) %>%
        full_join(singleton.AR.genes.count) %>%
        full_join(singleton.MGE.genes.count) %>%
        full_join(singleton.genes.count) %>%
        full_join(total.genes.count) %>%
        mutate(AR_duplicates = replace_na(AR_duplicates, 0)) %>%
        arrange(desc(AR_duplicates))
    
    return(big.gene.analysis.df)
}

big.gene.analysis.df <- make.big.gene.analysis.df(duplicate.proteins, singleton.proteins)
gc() ## free memory
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


make.Fig3.df <- function(ControlTable2, Table2, order.by.total_isolates.vec) {

    Fig3.data <- full_join(ControlTable2, Table2) %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total_isolates.vec)))

    ## See Figure 1B for the intuitive explanation for this code.
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


Fig3.df <- make.Fig3.df(ControlTable2, Table2, order.by.total_isolates.vec)

make.Fig3 <- function(Fig3.df) {

    ## calculate formal statistics, and use to add columns for coloring annotations.
    chromosome.ARG.duplicate.pvals <- Fig3.df %>%
        calc.chromosomal.ARG.duplicate.enrichment.pvals()    
    plasmid.ARG.duplicate.pvals <- Fig3.df %>%
        calc.plasmid.ARG.duplicate.enrichment.pvals()
    chromosome.ARG.singleton.pvals <- Fig3.df %>%
        calc.chromosomal.ARG.singleton.enrichment.pvals()
    plasmid.ARG.singleton.pvals <- Fig3.df %>%
        calc.plasmid.ARG.singleton.enrichment.pvals()
    
    Fig3.color.palette <- scales::viridis_pal()(3)
    
    Fig3A <- ggplot(Fig3.df, aes(x=chromosomal_duplicate_genes,
                                 y=chromosomal_duplicate_ARGs,
                                 label=Annotation)) +
        theme_classic() +
        geom_point(alpha=0.2, color = Fig3.color.palette[1]) +                  
        geom_line(aes(y=yvals.for.chromosomal_ARG_duplicates.line),
                  color = Fig3.color.palette[1]) +
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
        xlab("Chromosomal multi-copy genes") +
        ylab("Chromosomal multi-copy ARGs")
    
    Fig3B <- ggplot(Fig3.df, aes(x=plasmid_duplicate_genes,
                                 y=plasmid_duplicate_ARGs,
                                 label=Annotation)) +
        theme_classic() +
        geom_point(alpha=0.2, color = Fig3.color.palette[1]) +                  
        geom_line(aes(y=yvals.for.plasmid_ARG_duplicates.line),
                  color = Fig3.color.palette[1]) +
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
        xlab("Plasmid multi-copy genes") +
        ylab("Plasmid multi-copy ARGs")
    
    Fig3C <- ggplot(Fig3.df, aes(x=chromosomal_singleton_genes,
                                 y=chromosomal_singleton_ARGs,
                                 label=Annotation)) +
        theme_classic() +
        geom_point(alpha=0.2, color = Fig3.color.palette[2]) +                  
        geom_line(aes(y=yvals.for.chromosomal_ARG_singletons.line),
                  color = Fig3.color.palette[2]) +
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
        xlab("Chromosomal single-copy genes") +
        ylab("Chromosomal single-copy ARGs")
    
    Fig3D <- ggplot(Fig3.df, aes(x=plasmid_singleton_genes,
                                 y=plasmid_singleton_ARGs,
                                 label=Annotation)) +
        theme_classic() +
        geom_point(alpha=0.2, color = Fig3.color.palette[2]) +                  
        geom_line(aes(y=yvals.for.plasmid_ARG_singletons.line),
                  color = Fig3.color.palette[2]) +
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
        xlab("Plasmid single-copy genes") +
        ylab("Plasmid single-copy ARGs")
    Fig3 <- plot_grid(Fig3A, Fig3B, Fig3C, Fig3D, labels = c("A","B","C","D"), nrow=2)
    return(Fig3)
}

Fig3 <- make.Fig3(Fig3.df)
ggsave("../results/Fig3.pdf", Fig3)

## formal statistics for the results in Figure 3 (the colored points are enriched).
chromosome.ARG.duplicate.pvals <- Fig3.df %>%
    calc.chromosomal.ARG.duplicate.enrichment.pvals()    
plasmid.ARG.duplicate.pvals <- Fig3.df %>%
    calc.plasmid.ARG.duplicate.enrichment.pvals()
chromosome.ARG.singleton.pvals <- Fig3.df %>%
    calc.chromosomal.ARG.singleton.enrichment.pvals()
plasmid.ARG.singleton.pvals <- Fig3.df %>%
    calc.plasmid.ARG.singleton.enrichment.pvals()


#############################
## Supplementary Table S2. TODO: change to supplementary Table S4??
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

## I am playing around here.

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

genera.isolate.comparison.df <- full_join(
    duplicated.genera.isolate.summary,
    all.genera.isolate.summary) %>%
    ## turn NAs to zeros.
    replace(is.na(.), 0) %>%
    mutate(percent.genomes.with.dup.ARGs = duplicated.ARG.genome.count/genome.count) %>%
    arrange(desc(percent.genomes.with.dup.ARGs))
    
genera.duplicated.ARG.comparison.plot <- genera.isolate.comparison.df %>%
    ggplot(aes(x=log2(1+genome.count),
               y = log2(1+duplicated.ARG.genome.count),
               label = Genus)) +
    theme_classic() + geom_jitter() + geom_text_repel()

## Hmmm... sampling biases could be problematic.
## AKA antibiotic-resistance bacteria might be most likely to be sequenced.

unannotated.duplicate.proteins <- all.duplicate.proteins %>%
    filter(Annotation == "Unannotated")

unannotated.duplicated.ARG.isolate.genera.summary <- unannotated.duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(Annotation_Accession, Organism, Strain, Annotation) %>%
    distinct() %>%
    tibble() %>%
    mutate(Genus = stringr::word(Organism, 1)) %>%
    group_by(Genus) %>%
    summarize(duplicated.ARG.genome.count = n()) %>%
    arrange(desc(duplicated.ARG.genome.count))

################################################################################
## Figure 4: examples that indicate generality of our method.
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
        group_by(Annotation) %>%
        summarize(isolates_with_duplicated_function = n()) %>%
        arrange(desc(isolates_with_duplicated_function))
    
    ## join columns to make Table 1 with raw data.
    raw.Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_duplicated_function =
                   replace_na(isolates_with_duplicated_function,0)) %>%
        arrange(desc(isolates_with_duplicated_function))
    
    ## For consistency with Figure 1, use the distribution of sampled isolates as the null
    ## distribution.
    calc.expected.isolates.with.function <- function(raw.Table) {
        sum.total.isolates <- sum(raw.Table$total_isolates)
        total.isolates.with.duplicated.function <- sum(
            raw.Table$isolates_with_duplicated_function)
        Table <- raw.Table %>%
            mutate(expected_isolates_with_duplicated_function = total.isolates.with.duplicated.function * total_isolates/sum.total.isolates)
        return(Table)
    }
    
    calc.isolate.function.gene.enrichment.pvals <- function(raw.Table) {
        sum.total.isolates <- sum(raw.Table$total_isolates)
        total.isolates.with.duplicated.function <- sum(raw.Table$isolates_with_duplicated_function)
        Table <- raw.Table %>%
            rowwise() %>%
            mutate(binom.test.pval = binom.test
            (
                x = isolates_with_duplicated_function,
                n = total.isolates.with.duplicated.function,
                p = total_isolates/sum.total.isolates
            )$p.value) %>%
            mutate(corrected.pval = p.adjust(binom.test.pval,"BH")) %>%
            ## use Benjamini-Hochberg p-value correction.
            select(-binom.test.pval) ## drop original p-value after the correction.
        return(Table)
    }

    ## Add a third column: expected number of isolates with duplicated function genes,
    ## based on the percentage of isolates with duplicated genes.
    Table <- raw.Table %>% calc.expected.isolates.with.function() %>%
        ## Add a fourth column: p-values for deviation from
        ## expected number of duplicated function genes, using binomial test,
        ## correcting for multiple tests.
        calc.isolate.function.gene.enrichment.pvals()

    return(Table)
}


make.Fig4.panel.df <- function(Fig4.panel.table, pval.threshold = 0.05) {

    total_isolates.sum <- sum(Fig4.panel.table$total_isolates)
    isolates_with_duplicated_function.sum <- sum(Fig4.panel.table$isolates_with_duplicated_function)
    
    Fig4.panel.df <- Fig4.panel.table %>%
        ## calculate y-coordinates for line for duplicate ARGs.
        mutate(yvals.for.isolates_with_duplicated_function.line = total_isolates * isolates_with_duplicated_function.sum/total_isolates.sum) %>%
        ## add a column for label colors in the Figure.
        mutate(annotation_label_color = ifelse(
        (isolates_with_duplicated_function > expected_isolates_with_duplicated_function)
        && (corrected.pval < pval.threshold),
        "red", "black"))
    
    return(Fig4.panel.df)

}

make.Fig4.panel <- function(Fig4.panel.df) {

    Fig4.panel <- ggplot(Fig4.panel.df, aes(x=total_isolates,
                                 y=isolates_with_duplicated_function,
                                 label=Annotation)) +
        theme_classic() +
        geom_point(alpha=0.2) +                  
        geom_line(aes(y=yvals.for.isolates_with_duplicated_function.line)) +
        geom_text_repel(aes(color=annotation_label_color), size=3) +
        scale_color_identity() +
        scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        xlab("Total isolates")

    return(Fig4.panel)
}

photosynthesis.table <- make.IsolateEnrichmentTable(gbk.annotation,
                                                    duplicate.genes,
                                                    "photosystem")
write.csv(x=photosynthesis.table, file="../results/TableS4.csv")

N2.fixation.table <- make.IsolateEnrichmentTable(gbk.annotation,
                                                 duplicate.genes,
                                                 "nitrogenase")
write.csv(x=N2.fixation.table, file="../results/TableS5.csv")

toxic.metal.table <- make.IsolateEnrichmentTable(gbk.annotation,
                                                 duplicate.genes,
                                                 "mercury|cadmium|arsen")
write.csv(x=toxic.metal.table, file="../results/TableS6.csv")

heme.table <- make.IsolateEnrichmentTable(gbk.annotation,
                                          duplicate.genes,
                                          "heme")
write.csv(x=heme.table, file="../results/TableS7.csv")


Fig4A.df <- make.Fig4.panel.df(photosynthesis.table)
Fig4A <- make.Fig4.panel(Fig4A.df) +
    ylab("Isolates with multi-copy photosystem genes") +
    ggtitle("Photosynthesis")

Fig4B.df <- make.Fig4.panel.df(N2.fixation.table)
Fig4B <- make.Fig4.panel(Fig4B.df) +
    ylab("Isolates with multi-copy nitrogenase genes") +
    ggtitle("Nitrogen fixation")

Fig4C.df <- make.Fig4.panel.df(toxic.metal.table)
Fig4C <- make.Fig4.panel(Fig4C.df) +
    ylab("Isolates with multi-copy toxic-metal resistance genes") +
    ggtitle("Toxic-metal resistance")

Fig4D.df <- make.Fig4.panel.df(heme.table)
Fig4D <- make.Fig4.panel(Fig4D.df) +
    ylab("Isolates with multi-copy heme degradation genes") +
    ggtitle("Heme degradation")

Fig4 <- plot_grid(Fig4A, Fig4B, Fig4C, Fig4D,
                  labels = c("A","B","C","D"),
                  nrow = 2)

ggsave(Fig4, file = "../results/Fig4.pdf", width = 8.5, height = 8.5)
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

## Figure S1. annotations of multi-copy proteins are informative about
## ecology.

best.dup.prot.annotation.tf_idf <- dup.prot.annotation.tf_idf %>%
    filter(Annotation %in% c("Agriculture", "Anthropogenic-environment",
                             "Human-host")) %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 5) %>%
    ungroup()

S1Fig <- ggplot(best.dup.prot.annotation.tf_idf,
                aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
    geom_col(show.legend = FALSE) +
    labs(x = "tf-idf", y = NULL) +
    theme_classic() +
    xlim(0,0.04) +
    facet_wrap(.~Annotation, ncol=1, scales = "free_y") +
    ggtitle("Most informative multi-copy protein annotations")

ggsave("../results/S1Fig.pdf", S1Fig, width=8, height = 8)
