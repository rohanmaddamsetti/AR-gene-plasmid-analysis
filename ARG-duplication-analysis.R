## ARG-duplication-analysis.R by Rohan Maddamsetti.
## analyse the distribution of antibiotic resistance genes (ARGs)
## on chromosomes versus plasmids in  fully-sequenced genomes and plasmids
## in the NCBI Nucleotide database.
## This is a rewrite of Aim1-analysis.R, using a different statistical framework.

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)
library(tidytext) ## for text mining with R.
library(forcats)


fancy_scientific <- function(x) {
    ## function for plotting better axis labels.
    ## see solution here for nice scientific notation on axes.
    ## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

################################################################################
## Regular expressions used in this analysis.

## unknown protein keywords.
unknown.protein.keywords <- "unknown|Unknown|hypothetical|Hypothetical|Uncharacterized|Uncharacterised|uncharacterized|uncharacterised|DUF|unknow|putative protein in bacteria|Unassigned|unassigned"

## NOTE: some hypothetical proteins are "ISXX family insertion sequence hypothetical protein"
## so filter out those cases, when counting unknown proteins.


## match MGE genes using the following keywords in the "product" annotation
IS.keywords <- "IS|transposon|Transposase|transposase|Transposable|transposable|virus|Phage|phage|integrase|Integrase|baseplate|tail|intron|Mobile|mobile|antitoxin|toxin|capsid|plasmid|Plasmid|conjug|Tra"

MGE.or.unknown.protein.keywords <- paste(IS.keywords,unknown.protein.keywords,sep="|")

## Elongation Factor Tu (2 copies in most bacteria).
## \\b is a word boundary.
## see: https://stackoverflow.com/questions/62430498/detecting-whole-words-using-str-detect-in-r
EFTu.keywords <- "\\bTu | Tu\\b|-Tu\\b"

## antibiotic-specific keywords.
chloramphenicol.keywords <- "chloramphenicol|Chloramphenicol"
tetracycline.keywords <- "tetracycline|Tetracycline"
MLS.keywords <- "macrolide|lincosamide|streptogramin"
multidrug.keywords <- "multidrug"
beta.lactam.keywords <- "lactamase"
glycopeptide.keywords <- "glycopeptide resistance|VanZ"
polypeptide.keywords <- "bacitracin|polymyxin B"
diaminopyrimidine.keywords <- "trimethoprim-resistant"
sulfonamide.keywords <- "sulfonamide-resistant"
quinolone.keywords <- "quinolone|Quinolone|oxacin"
aminoglycoside.keywords <- "aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin"
macrolide.keywords <- "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythomycin"

antibiotic.keywords <- "chloramphenicol|Chloramphenicol|tetracycline|Tetracycline|macrolide|lincosamide|streptogramin|multidrug|lactamase|glycopeptide resistance|VanZ|bacitracin|polymyxin B|trimethoprim-resistant|sulfonamide-resistant|quinolone|Quinolone|oxacin|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythomycin|antibiotic resistance"

antibiotic.or.IS.keywords <- paste(IS.keywords,antibiotic.keywords,sep="|")

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

## import the 35GB file containing all proteins, including singletons.
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
## list the 2,848 strains missing from the singletons data.
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
## Figure 1 is a schematic of the analysis pipeline.
###########################################################################
## Analysis for Figure 2ABC.

## See Wikipedia reference:
## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

## Make Z-distributed confidence intervals for the fraction of isolates with
## duplicated ARGs (panel A),
## the fraction of isolates with single-copy ARGs (panel B),
## the fraction of isolates with duplicated genes (panel C).

## Count data for isolates with duplicated ARGs
## goes into Supplementary Table S1.

calc.isolate.confints <- function(df) {
    df %>%
        ## use the normal approximation for binomial proportion conf.ints
        mutate(se = sqrt(p*(1-p)/total_isolates)) %>%
        ## See Wikipedia reference:
        ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
        mutate(Left = p - 1.96*se) %>%
        mutate(Right = p + 1.96*se) %>%
        ## truncate confidence limits to interval [0,1].
        rowwise() %>% mutate(Left = max(0, Left)) %>%
        rowwise() %>% mutate(Right = min(1, Right)) %>%
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


## generic version of make.TableS1, for examining classes of genes other than
## antibiotic resistance genes.
make.IsolateEnrichmentTable <- function(gbk.annotation, duplicate.proteins, keywords) {
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


make.IsolateEnrichmentControlTable <- function(gbk.annotation, singleton.proteins, keywords) {
    ## count the number of isolates with singleton genes of interest in each category.
    category.counts <- singleton.proteins %>%
        filter(str_detect(.$product, keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_singleton_function = n)
    
    ## join columns to make the Table.
    Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_singleton_function =
                   replace_na(isolates_with_singleton_function,0)) %>%
        mutate(p = isolates_with_singleton_function/total_isolates) %>%
        calc.isolate.confints()
    return(Table)
}


make.confint.figure.panel <- function(Table, order.by.total.isolates, title,
                                      no.category.label = FALSE) {    
    Fig.panel <- Table %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total.isolates))) %>%
        ggplot(aes(y = Annotation, x = p)) +
        geom_point(size=1) +
        ylab("") +
        xlab("Proportion of Isolates") +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=0.2, size=0.2)
    
    if (no.category.label)
        Fig.panel <- Fig.panel +
            theme(axis.text.y=element_blank())
    
    return(Fig.panel)
}


## Figure 2A: normal-approximation confidence intervals for the percentage
## of isolates with duplicated ARGs.
TableS1 <- make.TableS1(gbk.annotation, duplicate.proteins)
## write Supplementary Table S1 to file.
write.csv(x=TableS1, file="../results/TableS1.csv")

Fig2A <- make.confint.figure.panel(TableS1, order.by.total.isolates, "Duplicated ARGs")

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

## This data frame will be used for Figure 2B.
TableS2 <- make.TableS2(gbk.annotation, singleton.proteins)
## write TableS2 to file.
write.csv(x=TableS2, file="../results/TableS2.csv")

Fig2B <- make.confint.figure.panel(TableS2, order.by.total.isolates,
                                   "Single-copy ARGs", no.category.label=TRUE)
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

Fig2C <- make.confint.figure.panel(TableS3, order.by.total.isolates,
                                   "All Duplicated Genes", no.category.label=TRUE)

Fig2ABC.title <- title_theme <- ggdraw() +
    draw_label("Isolate-level analysis",fontface="bold")

Fig2ABC <- plot_grid(Fig2A, Fig2B, Fig2C, labels=c('A','B','C'),
                     rel_widths = c(1.5,1,1), nrow=1)

Fig2ABC.with.title <- plot_grid(Fig2ABC.title, Fig2ABC, ncol = 1, rel_heights = c(0.2, 1))

## the rest of Figure 2 -- Fig2DEFGHI -- requires Tables S4 and S5.
## See below.

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

## Let's mark the genera with the most ARGs, and recalculate a version
## of TableS1, and the associated figure.
top.ARG.genera <- c("Klebsiella", "Escherichia",
                    "Acinetobacter", "Salmonella",
                    "Staphylococcus", "Enterobacter",
                    "Pseudomonas", "Proteus", "Citrobacter")

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
    theme_classic() + geom_jitter() + geom_text_repel(fontface = "italic") +
    scale_color_manual(values=c("black", "red")) +
    guides(color=FALSE) +
    xlab("sqrt(Number of isolates)") +
    ylab("sqrt(Number of isolates with duplicated ARGs)")

## 7,195 isolates are in the top ARG genera.
top.ARG.genera.isolates <- gbk.annotation %>%
    filter(Genus %in% top.ARG.genera)

## 11,509 isolates are in the remaining genera.
filtered.gbk.annotation <- gbk.annotation %>%
    filter(!(Genus %in% top.ARG.genera))

filtered.duplicate.proteins <- duplicate.proteins %>%
    filter(!(Genus %in% top.ARG.genera))

filtered.TableS1 <- make.TableS1(filtered.gbk.annotation, filtered.duplicate.proteins)
S1FigB <- make.confint.figure.panel(filtered.TableS1, order.by.total.isolates, "Duplicated ARGs after filtering top genera")

## Let's try an alternative strategy: downsample the data such that only one
## sample for each organism is allowed.

single.organism.gbk.annotation <- gbk.annotation %>%
    group_by(Organism) %>%
    filter(row_number() == 1) %>% ## take the first one in the group.
    ungroup()

single.organism.duplicate.proteins <- duplicate.proteins %>%
    filter(Annotation_Accession %in% single.organism.gbk.annotation$Annotation_Accession)

single.organism.singleton.proteins <- singleton.proteins %>%
    filter(Annotation_Accession %in% single.organism.gbk.annotation$Annotation_Accession)

single.organism.TableS1 <- make.TableS1(
    single.organism.gbk.annotation, single.organism.duplicate.proteins)
S1FigC <- make.confint.figure.panel(
    single.organism.TableS1, order.by.total.isolates, "Duplicated ARGs after downsampling species")

S1Fig <- plot_grid(S1FigA, S1FigB, S1FigC, labels = c("A", "B", "C"), nrow = 3, rel_heights = c(2,1,1))
ggsave("../results/S1Fig.pdf", S1Fig, height=10)

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
    ylab("") +
    xlab("Proportion of duplicated ARGs per genome") +
    theme_classic() +
    ggtitle("Duplicated ARGs")

S2FigB <- boxplot.df %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates))) %>%
    ggplot(aes(y = Annotation, x = proportion.of.duplicated.MGEs)) +
    geom_jitter(size=0.2) +
    ylab("") +
    xlab("Proportion of duplicated MGE genes per genome") +
    theme_classic() +
    ggtitle("Duplicated MGE genes")

S2FigC <- boxplot.df %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates))) %>%
    ggplot(aes(y = Annotation, x = proportion.of.duplicated.genes)) +
    geom_jitter(size=0.2) +
    ylab("") +
    xlab("Proportion of duplicated genes per genome") +
    theme_classic() +
    ggtitle("Duplicated genes")

S2Fig <- plot_grid(S2FigA, S2FigB, S2FigC, labels = c('A','B','C'), nrow = 3)
ggsave("../results/S2Fig.pdf", S2Fig)

################################################################################
## Figures S3 and S4. Proportion of isolates with duplicated or single-copy ARGs
## for 12 different antibiotic classes.

########################################
## Figure S3: Duplicated ARGs.


chloramphenicol.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    chloramphenicol.keywords)

tetracycline.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    tetracycline.keywords)

MLS.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    MLS.keywords)

multidrug.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    multidrug.keywords)

beta.lactam.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    beta.lactam.keywords)

glycopeptide.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    glycopeptide.keywords)

polypeptide.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    polypeptide.keywords)

diaminopyrimidine.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    diaminopyrimidine.keywords)

sulfonamide.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    sulfonamide.keywords)

quinolone.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    quinolone.keywords)

aminoglycoside.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    aminoglycoside.keywords)

macrolide.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    macrolide.keywords)


S3FigA <- make.confint.figure.panel(chloramphenicol.table,
                                    order.by.total.isolates,
                                    "chloramphenicol resistance")
S3FigB <- make.confint.figure.panel(tetracycline.table,
                                    order.by.total.isolates,
                                    "tetracycline resistance",
                                    no.category.label = TRUE)
S3FigC <- make.confint.figure.panel(MLS.table,
                                    order.by.total.isolates,
                                    "MLS resistance",
                                    no.category.label = TRUE)
S3FigD <- make.confint.figure.panel(multidrug.table,
                                    order.by.total.isolates,
                                    "multidrug resistance")
S3FigE <- make.confint.figure.panel(beta.lactam.table,
                                    order.by.total.isolates,
                                    "beta-lactam resistance",
                                    no.category.label = TRUE)
S3FigF <- make.confint.figure.panel(glycopeptide.table,
                                    order.by.total.isolates,
                                    "glycopeptide resistance",
                                    no.category.label = TRUE)
S3FigG <- make.confint.figure.panel(polypeptide.table,
                                    order.by.total.isolates,
                                    "polypeptide resistance")
S3FigH <- make.confint.figure.panel(diaminopyrimidine.table,
                                    order.by.total.isolates,
                                    "diaminopyrimidine resistance",
                                    no.category.label = TRUE)
S3FigI <- make.confint.figure.panel(sulfonamide.table,
                                    order.by.total.isolates,
                                    "sulfonamide resistance",
                                    no.category.label = TRUE)
S3FigJ <- make.confint.figure.panel(quinolone.table,
                                    order.by.total.isolates,
                                    "quinolone resistance")
S3FigK <- make.confint.figure.panel(aminoglycoside.table,
                                    order.by.total.isolates,
                                    "aminoglycoside resistance",
                                    no.category.label = TRUE)
S3FigL <- make.confint.figure.panel(macrolide.table,
                                    order.by.total.isolates,
                                    "macrolide resistance",
                                    no.category.label = TRUE)

S3Fig <- plot_grid(NULL, ## The nesting is to add a title.
                   plot_grid(
                       S3FigA, S3FigB, S3FigC,
                       S3FigD, S3FigE, S3FigF,
                       S3FigG, S3FigH, S3FigI,
                       S3FigJ, S3FigK, S3FigL,
                   rel_widths = c(1.5, 1, 1,
                                  1.5, 1, 1,
                                  1.5, 1, 1),
                   nrow = 4),
                   labels = c("Duplicated ARGs", ""),
                   ncol = 1,
                   rel_heights = c(0.05, 1))
ggsave("../results/S3Fig.pdf", S3Fig, height = 9.75, width=9.75)
########################################
## Figure S4: Single-copy ARGs.

chloramphenicol.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    chloramphenicol.keywords)

tetracycline.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    tetracycline.keywords)

MLS.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    MLS.keywords)

multidrug.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    multidrug.keywords)

beta.lactam.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    beta.lactam.keywords)

glycopeptide.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    glycopeptide.keywords)

polypeptide.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    polypeptide.keywords)

diaminopyrimidine.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    diaminopyrimidine.keywords)

sulfonamide.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    sulfonamide.keywords)

quinolone.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    quinolone.keywords)

aminoglycoside.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    aminoglycoside.keywords)

macrolide.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    macrolide.keywords)


S4FigA <- make.confint.figure.panel(chloramphenicol.control.table,
                                      order.by.total.isolates,
                                      "chloramphenicol resistance")
S4FigB <- make.confint.figure.panel(tetracycline.control.table,
                                      order.by.total.isolates,
                                      "tetracycline resistance",
                                      no.category.label = TRUE)
S4FigC <- make.confint.figure.panel(MLS.control.table,
                                      order.by.total.isolates,
                                      "MLS resistance",
                                      no.category.label = TRUE)
S4FigD <- make.confint.figure.panel(multidrug.control.table,
                                      order.by.total.isolates,
                                      "multidrug resistance")
S4FigE <- make.confint.figure.panel(beta.lactam.control.table,
                                      order.by.total.isolates,
                                      "beta-lactam resistance",
                                      no.category.label = TRUE)
S4FigF <- make.confint.figure.panel(glycopeptide.control.table,
                                      order.by.total.isolates,
                                      "glycopeptide resistance",
                                      no.category.label = TRUE)
S4FigG <- make.confint.figure.panel(polypeptide.control.table,
                                      order.by.total.isolates,
                                      "polypeptide resistance")
S4FigH <- make.confint.figure.panel(diaminopyrimidine.control.table,
                                      order.by.total.isolates,
                                      "diaminopyrimidine resistance",
                                      no.category.label = TRUE)
S4FigI <- make.confint.figure.panel(sulfonamide.control.table,
                                      order.by.total.isolates,
                                      "sulfonamide resistance",
                                      no.category.label = TRUE)
S4FigJ <- make.confint.figure.panel(quinolone.control.table,
                                    order.by.total.isolates,
                                    "quinolone resistance")
S4FigK <- make.confint.figure.panel(aminoglycoside.control.table,
                                      order.by.total.isolates,
                                      "aminoglycoside resistance",
                                      no.category.label = TRUE)
S4FigL <- make.confint.figure.panel(macrolide.control.table,
                                      order.by.total.isolates,
                                      "macrolide resistance",
                                      no.category.label = TRUE)

S4Fig <- plot_grid(NULL, ## The nesting is to add a title.
                   plot_grid(
                       S4FigA, S4FigB, S4FigC,
                       S4FigD, S4FigE, S4FigF,
                       S4FigG, S4FigH, S4FigI,
                       S4FigJ, S4FigK, S4FigL,
                   rel_widths = c(1.5, 1, 1,
                                  1.5, 1, 1,
                                  1.5, 1, 1),
                   nrow = 4),
                   labels = c("Single-copy ARGs", ""),
                   ncol = 1,
                   rel_heights = c(0.05, 1))
ggsave("../results/S4Fig.pdf", S4Fig, height = 9.75, width = 9.75)


##########################################################################
## Figure 3AB: Visualization of ARGs on plasmids and chromosomes.

categorize.as.MGE.ARG.or.other <- function(product) {
    if (is.na(product))
        return("Other function")
    else if (str_detect(product, antibiotic.keywords))
        return("ARG")
    else if (str_detect(product, IS.keywords))
        return("MGE")
    else
        return("Other function")
}


Fig3A.data <- duplicate.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation, Category) %>%
    summarize(Plasmid = sum(plasmid_count), Chromosome = sum(chromosome_count)) %>%
    pivot_longer(cols = c("Plasmid", "Chromosome"),
                 names_to = "Episome",
                 values_to = "Count") %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

Fig3B.data <- singleton.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation, Category) %>%
    summarize(Plasmid = sum(plasmid_count), Chromosome = sum(chromosome_count)) %>%
    pivot_longer(cols = c("Plasmid", "Chromosome"),
                 names_to = "Episome",
                 values_to = "Count") %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

Fig3A <- ggplot(Fig3A.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    facet_wrap(.~Episome) +
    theme_classic() +
    ggtitle("Distribution of duplicated genes") +
    xlab("Proportion of genes") +
    theme(legend.position="bottom") +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

Fig3legend <- get_legend(Fig3A)
Fig3A <- Fig3A + guides(fill = FALSE)

Fig3B <- ggplot(Fig3B.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    facet_wrap(.~Episome) +
    theme_classic() +
    ggtitle("Distribution of single-copy genes") +
    xlab("Proportion of genes") +
    guides(fill = FALSE) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

stackedbar.Fig3A <- ggplot(Fig3A.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity") +
    facet_wrap(.~Episome, scales = "free") +
    theme_classic() +
    ggtitle("Distribution of duplicated genes") +
    scale_x_continuous(labels=fancy_scientific) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

stackedbar.Fig3B <- ggplot(Fig3B.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity") +
    facet_wrap(.~Episome, scales = "free") +
    theme_classic() +
    ggtitle("Distribution of single-copy genes") +
    guides(fill = FALSE) +
    scale_x_continuous(labels=fancy_scientific) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.


Fig3AB <- plot_grid(Fig3A, Fig3B, Fig3legend, labels = c("A", "B"),
                  ncol = 1, rel_heights = c(1,1,0.25))
## Panel 3C will be added later. Needed for the full figure.

## This visualization is also useful.
stackedbar.Fig3AB <- plot_grid(stackedbar.Fig3A,
                               stackedbar.Fig3B, labels = c("A", "B"), ncol = 1)
ggsave("../results/stackedbar-Fig3AB.pdf", stackedbar.Fig3AB)

## remove the dataframes to save memory.
rm(Fig3A.data)
rm(Fig3B.data)
gc()

#########################
## S5 Figure.
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

S5FigA <- ggplot(both.chr.and.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("Both chromosome and plasmid") +
    theme(legend.position="bottom") +
    ylab("")

S5Fig.legend <- get_legend(S5FigA)
S5FigA <- S5FigA + guides(fill = FALSE)

S5FigB <- ggplot(just.chromosome.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("chromosome only") +
    guides(fill = FALSE) +
    ylab("")

S5FigC <- ggplot(just.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("plasmid only") +
    guides(fill = FALSE) +
    ylab("")

S5Fig <- plot_grid(NULL, S5FigA, S5FigB, S5FigC, S5Fig.legend, ncol = 1,
                   labels = c("Genomic distribution of duplicated genes", "A","B","C"),
                   rel_heights=c(0.2,1,1,1,0.25))
ggsave("../results/S5Fig.pdf", S5Fig)

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
## Use the data in Tables S4 and S5 to make Figure 2DEFGHI.
## The point of this figure is to show that the distribution of
## duplicated ARGs is not predicted by the distribution of single-copy ARGs
## in the ecological categories.

make.Fig2DEFGHI.df <- function(TableS1, TableS4, TableS5) {
    order.by.total_isolates <- TableS1$Annotation
    
    df <- full_join(TableS4, TableS5) %>%
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
        
    return(df)
}


Fig2DEFGHI.df <- make.Fig2DEFGHI.df(TableS1, TableS4, TableS5)
## Figure 2DEFGHI.
## Plot point estimates for the fraction of chromosomal genes that are
## the fraction of genes that are duplicated ARGs (panel D),
## the fraction of chromosomal genes that are duplicated ARGs (panel E),
## the fraction of plasmid genes that are duplicated ARGs (panel F),
## the fraction of genes that are single-copy ARGs (panel G).
## the fraction of chromosomal genes that are single-copy ARGs (panel H),
## the fraction of plasmid genes that are single-copy ARGs (panel I),

## Fig2D
Fig2D.df <- Fig2DEFGHI.df %>%
    mutate(p = total_duplicate_ARGs/(total_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_genes)) %>%
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = p - 1.96*se) %>%
    mutate(Right = p + 1.96*se) %>%
    ## truncate confidence limits to interval [0,1].
    rowwise() %>% mutate(Left = max(0, Left)) %>%
    rowwise() %>% mutate(Right = min(1, Right)) %>%
    select(Annotation, total_duplicate_ARGs, total_genes,
           p, Left, Right)


## Fig2E: the fraction of chromosomal genes that are duplicated ARGs.
Fig2E.df <- Fig2DEFGHI.df %>%
    mutate(p = chromosomal_duplicate_ARGs/(total_chromosomal_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_chromosomal_genes)) %>%
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = p - 1.96*se) %>%
    mutate(Right = p + 1.96*se) %>%
    ## truncate confidence limits to interval [0,1].
    rowwise() %>% mutate(Left = max(0, Left)) %>%
    rowwise() %>% mutate(Right = min(1, Right)) %>%
    select(Annotation, chromosomal_duplicate_ARGs, total_chromosomal_genes,
           p, Left, Right)

Fig2F.df <- Fig2DEFGHI.df %>%
    mutate(p = plasmid_duplicate_ARGs/(total_plasmid_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_plasmid_genes)) %>%
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = p - 1.96*se) %>%
    mutate(Right = p + 1.96*se) %>%
    ## truncate confidence limits to interval [0,1].
    rowwise() %>% mutate(Left = max(0, Left)) %>%
    rowwise() %>% mutate(Right = min(1, Right)) %>%
    select(Annotation, plasmid_duplicate_ARGs, total_plasmid_genes,
           p, Left, Right)

Fig2G.df <- Fig2DEFGHI.df %>%
    mutate(p = total_singleton_ARGs/(total_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_genes)) %>%
    ## and the Rule of Three to handle zeros.
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = p - 1.96*se) %>%
    mutate(Right = p + 1.96*se) %>%
    ## truncate confidence limits to interval [0,1].
    rowwise() %>% mutate(Left = max(0, Left)) %>%
    rowwise() %>% mutate(Right = min(1, Right)) %>%
    select(Annotation, total_singleton_ARGs, total_genes,
           p, Left, Right)

Fig2H.df <- Fig2DEFGHI.df %>%
    mutate(p = chromosomal_singleton_ARGs/(total_chromosomal_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_chromosomal_genes)) %>%
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = p - 1.96*se) %>%
    mutate(Right = p + 1.96*se) %>%
    ## truncate confidence limits to interval [0,1].
    rowwise() %>% mutate(Left = max(0, Left)) %>%
    rowwise() %>% mutate(Right = min(1, Right)) %>%
    select(Annotation, chromosomal_singleton_ARGs, total_chromosomal_genes,
           p, Left, Right)


Fig2I.df <- Fig2DEFGHI.df %>%
    mutate(p = plasmid_singleton_ARGs/(total_plasmid_genes)) %>%
    ## use the normal approximation for binomial proportion conf.ints
    mutate(se = sqrt(p*(1-p)/total_plasmid_genes)) %>%
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    mutate(Left = p - 1.96*se) %>%
    mutate(Right = p + 1.96*se) %>%
    ## truncate confidence limits to interval [0,1].
    rowwise() %>% mutate(Left = max(0, Left)) %>%
    rowwise() %>% mutate(Right = min(1, Right)) %>%
    select(Annotation, plasmid_singleton_ARGs, total_plasmid_genes,
           p, Left, Right)


make.Fig2DEFGHI.panel <- function(Table, order.by.total.isolates, title,
                            xlabel, no.category.label = FALSE) {
    panel <- Table %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total.isolates))) %>%
        ggplot(aes(y = Annotation, x = p)) +     
        geom_point(size=1) +
        ylab("") +
        xlab(xlabel) +
        scale_x_continuous(label=fancy_scientific) +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=0.2, size=0.2)
    
    if (no.category.label)
        panel <- panel +
            theme(axis.text.y=element_blank())
    
    return(panel)
}


Fig2D <- make.Fig2DEFGHI.panel(Fig2D.df, order.by.total.isolates,
                         "Duplicated ARGs",
                         "Proportion of all genes") +
    ## so that axis labels don't get cut off
    scale_x_continuous(label=fancy_scientific, limits = c(0,2.15e-4))
Fig2E <- make.Fig2DEFGHI.panel(Fig2E.df, order.by.total.isolates,
                         "Chromosome:\nDuplicated ARGs",
                         "Proportion of chromosomal genes",
                         no.category.label = TRUE)
Fig2F <- make.Fig2DEFGHI.panel(Fig2F.df, order.by.total.isolates,
                         "Plasmids:\nDuplicated ARGs",
                         "Proportion of plasmid genes",
                         no.category.label = TRUE)
Fig2G <- make.Fig2DEFGHI.panel(Fig2G.df, order.by.total.isolates,
                         "Single-copy ARGs",
                         "Proportion of all genes")
Fig2H <- make.Fig2DEFGHI.panel(Fig2H.df, order.by.total.isolates,
                         "Chromosome:\nSingle-copy ARGs",
                         "Proportion of chromosomal genes",
                         no.category.label = TRUE)
Fig2I <- make.Fig2DEFGHI.panel(Fig2I.df, order.by.total.isolates,
                         "Plasmids:\nSingle-copy ARGs",
                         "Proportion of plasmid genes",
                         no.category.label = TRUE)

Fig2DEFGHI.title <- title_theme <- ggdraw() +
    draw_label("Gene-level analysis",fontface="bold")

Fig2DEFGHI <- plot_grid(Fig2D, Fig2E, Fig2F, Fig2G, Fig2H, Fig2I,
                  labels = c('D','E','F','G','H','I'), nrow=2,
                  rel_widths = c(1.5, 1, 1, 1.5, 1, 1))

Fig2DEFGHI.with.title <- plot_grid(Fig2DEFGHI.title, Fig2DEFGHI,
                                   ncol = 1, rel_heights = c(0.2, 1))


## Now, make the complete Figure 2!
Fig2 <- plot_grid(Fig2ABC.with.title, Fig2DEFGHI.with.title,
                  ncol = 1, rel_heights = c(0.33,0.66))

ggsave("../results/Fig2.pdf", Fig2, height=9, width=10)

################################################################################
## Figure 3C: 
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

## A schematic figure in Illustrator to show the rationale goes into the
## Supplementary Figures.

make.selection.test.df <- function(duplicate.proteins, singleton.proteins,
                                   order.by.total.isolates,
                                   keywords, negate = FALSE) {

    if (negate) { ## then negate the str_detect for the keywords.
        duplicated.function.per.category <- duplicate.proteins %>%
            filter(!str_detect(.$product, keywords)) %>%
            group_by(Annotation) %>%
            summarize(function.duplicates = sum(count))
        
        singleton.function.per.category <- singleton.proteins %>%
            filter(!str_detect(.$product, keywords)) %>%
            group_by(Annotation) %>%
            summarize(function.singletons = sum(count))
    } else {
        duplicated.function.per.category <- duplicate.proteins %>%
            filter(str_detect(.$product, keywords)) %>%
            group_by(Annotation) %>%
            summarize(function.duplicates = sum(count))
        
        singleton.function.per.category <- singleton.proteins %>%
            filter(str_detect(.$product, keywords)) %>%
            group_by(Annotation) %>%
            summarize(function.singletons = sum(count))
    }

    duplicated.genes.per.category <- duplicate.proteins %>%
        group_by(Annotation) %>%
        summarize(gene.duplicates = sum(count))
    
    singleton.proteins.per.category <- singleton.proteins %>%
        group_by(Annotation) %>%
        summarize(singleton.proteins = sum(count))
    
    selection.test.df <- duplicated.function.per.category %>%
        full_join(duplicated.genes.per.category) %>%
        full_join(singleton.function.per.category) %>%
        full_join(singleton.proteins.per.category) %>%
        ## turn NAs to zeros.
        replace(is.na(.), 0) %>%
        mutate(p = function.duplicates / gene.duplicates) %>%
        mutate(q = function.singletons / singleton.proteins) %>%
        mutate(dup.singleton.ratio = p/q) %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total.isolates)))
    return(selection.test.df)
}


ARG.selection.test.df <- make.selection.test.df(
    duplicate.proteins,
    singleton.proteins,
    order.by.total.isolates,
    antibiotic.keywords) %>%
    select(Annotation, dup.singleton.ratio) %>%
    mutate(Category = "ARG")

MGE.selection.test.df <- make.selection.test.df(
    duplicate.proteins,
    singleton.proteins,
    order.by.total.isolates,
    IS.keywords) %>%
    select(Annotation, dup.singleton.ratio) %>%
    mutate(Category = "MGE")

other.selection.test.df <- make.selection.test.df(
    duplicate.proteins,
    singleton.proteins,
    order.by.total.isolates,
    antibiotic.or.IS.keywords, negate = TRUE) %>%
    select(Annotation, dup.singleton.ratio) %>%
    mutate(Category = "Other function")


big.selection.test.df <- ARG.selection.test.df %>%
    full_join(MGE.selection.test.df) %>%
    full_join(other.selection.test.df)

Fig3C <- big.selection.test.df %>%
    ggplot(aes(y = Annotation, x = dup.singleton.ratio, color = Category)) +
    geom_point() + theme_classic() +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    xlim(0,33) +
    xlab("proportion of genes in category among duplicate genes / proportion of genes in category among single-copy genes") +
    guides(color=FALSE) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

Fig3 <- plot_grid(Fig3AB, Fig3C, labels = c('','C'), nrow = 2, rel_heights=c(0.7,0.3))
ggsave("../results/Fig3.pdf", Fig3, width=10, height=7)

################################################################################
## Figure 6. A deterministic ODE model demonstrates that selection can
## drive the evolution of duplicated ARGs on plasmids.

## The panels of this figure are generated in my Pluto notebook:
## duplication-linear-ODE-model.jl.
################################################################################
## Figure 7: examples that indicate generality of our method.
## let's examine some other functions that we expect to be enriched in some, but
## not all ecological annotations.


## Tables for duplicate genes.
photosynthesis.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    "photosystem")
write.csv(x=photosynthesis.table, file="../results/TableS6.csv")

N2.fixation.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    "nitrogenase")
write.csv(x=N2.fixation.table, file="../results/TableS7.csv")

toxic.metal.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    "mercury|cadmium|arsen")
write.csv(x=toxic.metal.table, file="../results/TableS8.csv")

heme.table <- make.IsolateEnrichmentTable(
    gbk.annotation,
    duplicate.proteins,
    "heme")
write.csv(x=heme.table, file="../results/TableS9.csv")

## Tables for single-copy genes.
photosynthesis.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    "photosystem")

N2.fixation.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    "nitrogenase")

toxic.metal.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    "mercury|cadmium|arsen")

heme.control.table <- make.IsolateEnrichmentControlTable(
    gbk.annotation,
    singleton.proteins,
    "heme")

## tests for selection.
photosynthesis.selection.test.df <- make.selection.test.df(
    duplicate.proteins,
    singleton.proteins,
    order.by.total.isolates,
    "photosystem")

N2.fixation.selection.test.df <- make.selection.test.df(
    duplicate.proteins,
    singleton.proteins,
    order.by.total.isolates,
    "nitrogenase")

toxic.metal.selection.test.df <- make.selection.test.df(
    duplicate.proteins,
    singleton.proteins,
    order.by.total.isolates,
    "mercury|cadmium|arsen")

heme.selection.test.df <- make.selection.test.df(
    duplicate.proteins,
    singleton.proteins,
    order.by.total.isolates,
    "heme")

## plots of the tests for selection.
photosynthesis.selection.plot <- photosynthesis.selection.test.df %>%
    ggplot(aes(y = Annotation, x = dup.singleton.ratio)) +
    geom_point() + theme_classic() +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    xlab("proportion of function among duplicate genes / proportion of function among single-copy genes") +
    ggtitle("Photosynthesis")


N2.fixation.selection.plot <- N2.fixation.selection.test.df %>%
    ggplot(aes(y = Annotation, x = dup.singleton.ratio)) +
    geom_point() + theme_classic() +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    xlab("proportion of function among duplicate genes / proportion of function among single-copy genes") +
    ggtitle("Nitrogen fixation")

toxic.metal.selection.plot <- toxic.metal.selection.test.df %>%
    ggplot(aes(y = Annotation, x = dup.singleton.ratio)) +
    geom_point() + theme_classic() +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    xlab("proportion of function among duplicate genes / proportion of function among single-copy genes") +
    ggtitle("Toxic metal resistance")

heme.selection.plot <- heme.selection.test.df %>%
    ggplot(aes(y = Annotation, x = dup.singleton.ratio)) +
    geom_point() + theme_classic() +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    xlab("proportion of function among duplicate genes / proportion of function among single-copy genes") +
    ggtitle("Heme metabolism")

ggsave("../results/photosynthesis-selection-plot.pdf",
       height = 3, width = 7,
       photosynthesis.selection.plot)
ggsave("../results/nitrogenase-selection-plot.pdf",
       height = 3, width = 7,
       N2.fixation.selection.plot)
ggsave("../results/toxic-metal-selection-plot.pdf",
       height = 3, width = 7,
       toxic.metal.selection.plot)
ggsave("../results/heme-selection-plot.pdf",
       height = 3, width = 7,
       heme.selection.plot)

##########################################################################
## Make current version of Figure 7.
Fig7A <- make.confint.figure.panel(
    photosynthesis.table, order.by.total.isolates,
    "Duplicated photosystem genes")

Fig7B <- make.confint.figure.panel(
    photosynthesis.control.table, order.by.total.isolates,
    "Single-copy photosystem genes",
    no.category.label = TRUE)

Fig7C <- make.confint.figure.panel(
    N2.fixation.table, order.by.total.isolates,
    "Duplicated nitrogenases")

Fig7D <- make.confint.figure.panel(
    N2.fixation.control.table, order.by.total.isolates,
    "Single-copy nitrogenases",
    no.category.label = TRUE)

Fig7E <- make.confint.figure.panel(
    toxic.metal.table, order.by.total.isolates,
    "Duplicated toxic-metal resistance")

Fig7F <- make.confint.figure.panel(
    toxic.metal.control.table, order.by.total.isolates,
    "Single-copy toxic-metal resistance",
    no.category.label = TRUE)

Fig7G <- make.confint.figure.panel(
    heme.table, order.by.total.isolates,
    "Duplicated heme metabolism genes")

Fig7H <- make.confint.figure.panel(
    heme.control.table, order.by.total.isolates,
    "Single-copy heme metabolism genes",
    no.category.label = TRUE)

Fig7 <- plot_grid(Fig7A, Fig7B, Fig7C, Fig7D, Fig7E, Fig7F, Fig7G, Fig7H,
                  labels = c("A","B","C","D", "E", "F", "G", "H"),
                  ncol = 2,
                  rel_widths = c(1.4, 1, 1.4, 1, 1.4, 1, 1.4, 1))

ggsave(Fig7, file = "../results/Fig7.pdf", width = 8, height = 8)
##########################################################################
## run additional controls for taxonomy, as before.

## Control for Taxonomy (simpler control than for phylogeny)

## let's look at the taxonomic distribution of strains with duplicated function genes.
make.taxonomy.control.figure <- function(duplicate.proteins, gbk.annotation, keywords, title.text) {

    ## Let's mark the most commonly sampled genera, and recalculate a version
    ## of TableS1, and the associated figure.
    top.photosystem.genera <- c("Nostoc", "Synechococcus")
    top.N2.fixation.genera <- c("Rhizobium")
    top.toxic.metal.genera <- c("Klebsiella", "Citrobacter","Enterobacter",
                                "Escherichia", "Salmonella", "Pseudomonas",
                                "Staphylococcus")
    top.heme.genera <- c("Salmonella")

    gene.type <- keywords ## change for toxic metal.
    
    if (keywords == "photosystem") {
        top.func.genera <- top.photosystem.genera
    } else if (keywords == "nitrogenase") {
        top.func.genera <- top.N2.fixation.genera
    } else if (keywords == "mercury|cadmium|arsen") {
        top.func.genera <- top.toxic.metal.genera
        gene.type = "toxic metal"
    } else if (keywords == "heme") {
        top.func.genera <- top.heme.genera
    }
    
    duplicated.func.seq.genera.summary <- duplicate.proteins %>%
        tibble() %>% filter(str_detect(product,keywords)) %>%
        mutate(Genus = stringr::word(Organism, 1)) %>%
        group_by(Genus) %>%
        summarize(duplicated.func.count = n()) %>%
        arrange(desc(duplicated.func.count))
    
    duplicated.genera.seq.summary <- duplicate.proteins %>%
        tibble() %>%
        mutate(Genus = stringr::word(Organism, 1)) %>%
        group_by(Genus) %>%
        summarize(duplicated.seq.count = n()) %>%
        arrange(desc(duplicated.seq.count))
    
    all.genera.isolate.summary <- gbk.annotation %>%
        tibble() %>%
        mutate(Genus = stringr::word(Organism, 1)) %>%
        group_by(Genus) %>%
        summarize(genome.count = n()) %>%
        arrange(desc(genome.count))
    
    duplicated.genera.isolate.summary <- duplicate.proteins %>%
        filter(str_detect(.$product,keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Organism, Strain, Annotation) %>%
        distinct() %>%
        tibble() %>%
        mutate(Genus = stringr::word(Organism, 1)) %>%
        group_by(Genus) %>%
        summarize(duplicated.func.genome.count = n()) %>%
        arrange(desc(duplicated.func.genome.count))
    
    
    genera.isolate.comparison.df <- full_join(
        duplicated.genera.isolate.summary,
        all.genera.isolate.summary) %>%
        ## turn NAs to zeros.
        replace(is.na(.), 0) %>%
        mutate(percent.genomes.with.dup.funcs = duplicated.func.genome.count/genome.count) %>%
        mutate(in.top.func.genera = ifelse(Genus %in% top.func.genera, TRUE, FALSE))

    ## Hmmm... sampling biases could be problematic.
    FigA <- genera.isolate.comparison.df %>%
        ggplot(aes(x=sqrt(genome.count),
                   y = sqrt(duplicated.func.genome.count),
                   label = Genus,
                   color = in.top.func.genera)) +
        theme_classic() + geom_jitter() + geom_text_repel(fontface = "italic") +
        scale_color_manual(values=c("black", "red")) +
        guides(color=FALSE) +
        xlab("sqrt(Number of isolates)") +
        ylab("sqrt(Number of isolates with duplicated genes)")

    top.func.genera.isolates <- gbk.annotation %>%
        filter(Genus %in% top.func.genera)
    
    filtered.gbk.annotation <- gbk.annotation %>%
        filter(!(Genus %in% top.func.genera))
    
    filtered.duplicate.proteins <- duplicate.proteins %>%
        filter(!(Genus %in% top.func.genera))
    
    filtered.func.Table <- make.IsolateEnrichmentTable(filtered.gbk.annotation,
                                                       filtered.duplicate.proteins,
                                                       keywords)
    
    FigB <- make.confint.figure.panel(filtered.func.Table,
                                      order.by.total.isolates,
                                      "Duplicated genes after filtering top genera")

    ## Let's try an alternative strategy: downsample the data such that only one
    ## sample for each organism is allowed.
    single.organism.gbk.annotation <- gbk.annotation %>%
        group_by(Organism) %>%
        filter(row_number() == 1) %>% ## take the first one in the group.
        ungroup()
    
    single.organism.duplicate.proteins <- duplicate.proteins %>%
        filter(Annotation_Accession %in% single.organism.gbk.annotation$Annotation_Accession)
    
    single.organism.func.table <- make.IsolateEnrichmentTable(
    single.organism.gbk.annotation,
    single.organism.duplicate.proteins,
    keywords)

    FigC <- make.confint.figure.panel(
        single.organism.func.table, order.by.total.isolates,
    "Duplicated genes after downsampling species")

    title <- ggdraw() + draw_label(title.text, fontface='bold')
    
    FullFig <- plot_grid(title, FigA, FigB, FigC, labels = c("","A", "B", "C"),
                         nrow=4, rel_heights=c(0.1,1.5,1,1))
    return(FullFig)
}


S8Fig <- make.taxonomy.control.figure(
    duplicate.proteins, gbk.annotation, "photosystem",
    "Duplicated photosystem genes")
ggsave("../results/S8Fig.pdf", S8Fig, height=10)


S9Fig <- make.taxonomy.control.figure(
    duplicate.proteins, gbk.annotation, "nitrogenase",
    "Duplicated nitrogenase genes")
ggsave("../results/S9Fig.pdf",S9Fig, height=10)


S10Fig <- make.taxonomy.control.figure(
    duplicate.proteins, gbk.annotation, "mercury|cadmium|arsen",
    "Duplicated toxic metal resistance genes")
ggsave("../results/S10Fig.pdf",S10Fig, height=10)


S11Fig <- make.taxonomy.control.figure(
    duplicate.proteins, gbk.annotation, "heme",
    "Duplicated heme metabolism genes")
ggsave("../results/S11Fig.pdf", S11Fig, height=10)


## Make the single-organism filtered tables, in order to look at the
## numbers themselves.
single.organism.photosynthesis.table <- make.IsolateEnrichmentTable(
    single.organism.gbk.annotation,
    single.organism.duplicate.proteins,
    "photosystem")

single.organism.N2.fixation.table <- make.IsolateEnrichmentTable(
    single.organism.gbk.annotation,
    single.organism.duplicate.proteins,
    "nitrogenase")

single.organism.toxic.metal.table <- make.IsolateEnrichmentTable(
    single.organism.gbk.annotation,
    single.organism.duplicate.proteins,
    "mercury|cadmium|arsen")

single.organism.heme.table <- make.IsolateEnrichmentTable(
    single.organism.gbk.annotation,
    single.organism.duplicate.proteins,
    "heme")

##########################################################################
## Figure S6 and Figure S7:

## Lingchong suggested the following central point for the manuscript:
## The propensity for jumping to a plasmid is proportional to the propensity
## for gene duplication.

## Can we expand this analysis to 100,000s of genomes?
## Can we test this pattern using any given gene in a pangenome?
## For every type and function: expect to see a correlation between
## duplication and jumping to the plasmid.

## First let's make a figure for the correlation between isolates with duplicate
## ARGs and isolates with ARGs on plasmids.

make.plasmid.duplication.cor.plot <- function(duplicate.proteins,
                                              singleton.proteins,
                                              gbk.annotation,
                                              order.by.total.isolates,
                                              keywords,
                                              xlabel,
                                              ylabel) {

    ## count the number of isolates with duplicated functions in each category.
    function.category.counts <- duplicate.proteins %>%
        filter(str_detect(.$product, keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_duplicated_functions = n)

    isolates.with.plasmid.singleton.functions <- singleton.proteins %>%
        filter(plasmid_count > 0) %>%
        filter(str_detect(.$product, keywords)) %>%
        count(Annotation_Accession) %>%
        rename(singleton_plasmid_functions = n)

    isolates.with.plasmid.duplicate.functions <- duplicate.proteins %>%
        filter(plasmid_count > 0) %>%
        filter(str_detect(.$product, keywords)) %>%
        count(Annotation_Accession) %>%
        rename(duplicate_plasmid_functions = n)
    
    isolates.with.plasmid.functions.vec <- unique(c(
        isolates.with.plasmid.singleton.functions$Annotation_Accession,
        isolates.with.plasmid.duplicate.functions$Annotation_Accession))

    isolates.with.plasmid.functions.df <- gbk.annotation %>%
        filter(Annotation_Accession %in% isolates.with.plasmid.functions.vec) %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_plasmid_functions = n)
    
    ## join columns to make the new Table.
    Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(function.category.counts) %>%
        left_join(isolates.with.plasmid.functions.df) %>%
        mutate(isolates_with_duplicated_functions =
                   replace_na(isolates_with_duplicated_functions,0)) %>%
        mutate(p_dup = isolates_with_duplicated_functions/total_isolates) %>%
        ## use the normal approximation for binomial proportion conf.ints
        mutate(dup_se = sqrt(p_dup*(1-p_dup)/total_isolates)) %>%
        ## See Wikipedia reference:
        ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
        mutate(dup.Left = p_dup - 1.96*dup_se) %>%
        mutate(dup.Right = p_dup + 1.96*dup_se) %>%
        ## truncate confidence limits to interval [0,1].
        rowwise() %>% mutate(dup.Left = max(0, dup.Left)) %>%
        rowwise() %>% mutate(dup.Right = min(1, dup.Right)) %>%
        mutate(p_plasmid = isolates_with_plasmid_functions/total_isolates) %>%
        ## use the normal approximation for binomial proportion conf.ints
        mutate(plas_se = sqrt(p_plasmid*(1-p_plasmid)/total_isolates)) %>%
        ## See Wikipedia reference:
        ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
        mutate(plas.Left = p_plasmid - 1.96*plas_se) %>%
        mutate(plas.Right = p_plasmid + 1.96*plas_se) %>%
        ## truncate confidence limits to interval [0,1].
        rowwise() %>% mutate(plas.Left = max(0, plas.Left)) %>%
        rowwise() %>% mutate(plas.Right = min(1, plas.Right))
    
    
    plasmid_dup_correlation_plot <- Table %>%
        ggplot(aes(x = p_plasmid, y = p_dup, label = Annotation)) +
        geom_point() +
        ## plot both CIs.
        geom_errorbar(aes(ymin=dup.Left,ymax=dup.Right), width=0, size=0.2) +
        geom_errorbarh(aes(xmin=plas.Left,xmax=plas.Right), height=0, size=0.2) +
        geom_text_repel() + theme_classic() +
        xlab(xlabel) +
        ylab(ylabel)

    return(plasmid_dup_correlation_plot)
}


ARG_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    duplicate.proteins,
    singleton.proteins,
    gbk.annotation,
    order.by.total.isolates,
    antibiotic.keywords,
    "proportion of isolates with plasmid-borne ARGs",
    "proportion of isolates with duplicated ARGs") +
    ggtitle("ARGs")

photosynthesis_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    duplicate.proteins,
    singleton.proteins,
    gbk.annotation,
    order.by.total.isolates,
    "photosystem",
    "proportion of isolates with plasmid-borne photosystems",
    "proportion of isolates with duplicated photosystems") +
    ggtitle("Photosystems")

N2_fixation_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    duplicate.proteins,
    singleton.proteins,
    gbk.annotation,
    order.by.total.isolates,
    "nitrogenase",
    "proportion of isolates with plasmid-borne nitrogenases",
    "proportion of isolates with duplicated nitrogenases") +
    ggtitle("Nitrogenases")

toxic_metal_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    duplicate.proteins,
    singleton.proteins,
    gbk.annotation,
    order.by.total.isolates,
    "mercury|cadmium|arsen",
    "proportion of isolates with plasmid-borne toxic metal resistances",
    "proportion of isolates with duplicated toxic metal resistances") +
    ggtitle("Toxic metal resistance")

heme_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    duplicate.proteins,
    singleton.proteins,
    gbk.annotation,
    order.by.total.isolates,
    "heme",
    "proportion of isolates with plasmid-borne heme metabolism genes",
    "proportion of isolates with duplicated heme metabolism genes") +
    ggtitle("Heme metabolism")


S6Fig <- plot_grid(ARG_plasmid_dup_correlation_plot,
                   photosynthesis_plasmid_dup_correlation_plot,
                   N2_fixation_plasmid_dup_correlation_plot,
                   toxic_metal_plasmid_dup_correlation_plot,
                   heme_plasmid_dup_correlation_plot,
                   nrow = 3)
ggsave("../results/S6Fig.pdf", S6Fig, height=15, width=12)


## now, remake these plots, using the species-downsampled data.
single.organism.ARG_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    single.organism.duplicate.proteins,
    single.organism.singleton.proteins,
    single.organism.gbk.annotation,
    order.by.total.isolates,
    antibiotic.keywords,
    "proportion of isolates with plasmid-borne ARGs",
    "proportion of isolates with duplicated ARGs") +
    ggtitle("ARGs")

single.organism.photosynthesis_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    single.organism.duplicate.proteins,
    single.organism.singleton.proteins,
    single.organism.gbk.annotation,
    order.by.total.isolates,
    "photosystem",
    "proportion of isolates with plasmid-borne photosystems",
    "proportion of isolates with duplicated photosystems") +
    ggtitle("Photosystems")

single.organism.N2_fixation_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    single.organism.duplicate.proteins,
    single.organism.singleton.proteins,
    single.organism.gbk.annotation,
    order.by.total.isolates,
    "nitrogenase",
    "proportion of isolates with plasmid-borne nitrogenases",
    "proportion of isolates with duplicated nitrogenases") +
    ggtitle("Nitrogenases")

single.organism.toxic_metal_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    single.organism.duplicate.proteins,
    single.organism.singleton.proteins,
    single.organism.gbk.annotation,
    order.by.total.isolates,
    "mercury|cadmium|arsen",
    "proportion of isolates with plasmid-borne toxic metal resistances",
    "proportion of isolates with duplicated toxic metal resistances") +
    ggtitle("Toxic metal resistance")

single.organism.heme_plasmid_dup_correlation_plot <- make.plasmid.duplication.cor.plot(
    single.organism.duplicate.proteins,
    single.organism.singleton.proteins,
    single.organism.gbk.annotation,
    order.by.total.isolates,
    "heme",
    "proportion of isolates with plasmid-borne heme metabolism genes",
    "proportion of isolates with duplicated heme metabolism genes") +
    ggtitle("Heme metabolism")


S7Fig <- plot_grid(single.organism.ARG_plasmid_dup_correlation_plot,
                   single.organism.photosynthesis_plasmid_dup_correlation_plot,
                   single.organism.N2_fixation_plasmid_dup_correlation_plot,
                   single.organism.toxic_metal_plasmid_dup_correlation_plot,
                   single.organism.heme_plasmid_dup_correlation_plot,
                   nrow = 3)
ggsave("../results/S7Fig.pdf", S7Fig, height=15, width=12)

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
## function for doing TF-IDF for duplicate proteins on plasmids.
plasmid.duplicate.proteins <- duplicate.proteins %>% filter(plasmid_count >= 1)
make.plas.dup.annotation.freq.table <- partial(.f = .make.annotation.freq.table, duplicate.proteins)


## let's make a big table of product annotations per annotation category,
## for TF-IDF analysis.

## The analysis here closely follows the text mining example here:
## https://www.tidytextmining.com/tfidf.html

## run for duplicate proteins on plasmids.
big.plas.dup.prot.annotation.freq.table <- map_dfr(unique(plasmid.duplicate.proteins$Annotation),
                                              .f = make.plas.dup.annotation.freq.table) %>%
    ungroup() %>% ## have to ungroup before summing up all annotations in the table.
    mutate(total.annotation.count = sum(annotation.count))

plas.dup.prot.annotation.tf_idf <- big.plas.dup.prot.annotation.freq.table %>%
  bind_tf_idf(product, Annotation, annotation.count) %>%
  arrange(desc(tf_idf))

plas.top.dup.prot.annotation.tf_idf <- plas.dup.prot.annotation.tf_idf %>%
    group_by(Annotation) %>%
    slice_max(tf_idf, n = 5) %>%
    ungroup()

plas.dup.prot.annotation.tf_idf.plot <- plas.top.dup.prot.annotation.tf_idf %>%
    ggplot(aes(tf_idf, fct_reorder(product, tf_idf), fill = Annotation)) +
  geom_col(show.legend = TRUE) +
  facet_wrap(~Annotation, ncol = 2, scales = "free") +
  labs(x = "tf-idf", y = NULL)

ggsave("../results/plasmid-duplicate-protein-annotation-TF-IDF.pdf",
       plas.dup.prot.annotation.tf_idf.plot,
       height=21,width=21)


## run for all duplicate proteins.
big.dup.prot.annotation.freq.table <- map_dfr(
    unique(duplicate.proteins$Annotation),
    .f = make.dup.annotation.freq.table) %>%
    ## have to ungroup before summing up all annotations in the table.
    ungroup() %>% 
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
    ## have to ungroup before summing up all annotations in the table.
    ungroup() %>% 
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


#######################################################################################
## Figure 8.
## Let's analyze duplicated genes in the 12 GN0XXXX genomes that were sequenced with
## long-read technology (PacBio) by Vance Fowler's lab.
## Jon Bethke characterized the resistances of these strains, and additionally
## sequenced another 7 ESBL genomes with PacBio technology.

## 58,092 types of protein sequences in the 12 ESBL genomes.
Duke.ESBL.all.proteins <- data.table::fread("../results/Duke-ESBL-all-proteins.csv",
                                  drop="sequence")

## 57,263 single-copy protein sequences.
Duke.ESBL.singleton.proteins <- Duke.ESBL.all.proteins %>%
    filter(count == 1)

## 829 protein sequences have duplicates in the 12 ESBL genomes (not counting the duplicates)
Duke.ESBL.duplicate.proteins <- read.csv(
    "../results/Duke-ESBL-duplicate-proteins.csv") %>%
    select(-sequence)

## 6 of the 12 strains have duplicated ARGs.
## 10 ARG sequences have duplicates.
Duke.ESBL.duplicate.ARGs <- Duke.ESBL.duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords))
sum(Duke.ESBL.duplicate.ARGs$count) ## 23 duplicated ARGs in total.

## 215 MGE protein sequences are duplicated.
Duke.ESBL.duplicate.MGE.proteins <- Duke.ESBL.duplicate.proteins %>%
    filter(str_detect(.$product,IS.keywords))
sum(Duke.ESBL.duplicate.MGE.proteins$count) ## 547 duplicated MGE proteins in total.

## 208 unknown protein sequences are duplicated.
Duke.ESBL.duplicate.unknown.proteins <- Duke.ESBL.duplicate.proteins %>%
    ## some hypothetical proteins are "ISXX family insertion sequence hypothetical protein"
    ## so filter out those cases.
    filter(!str_detect(.$product,IS.keywords)) %>%
    filter(str_detect(.$product,unknown.protein.keywords))
sum(Duke.ESBL.duplicate.unknown.proteins$count) ## 450 duplicated unknown proteins in total.

## 396 remaining cases of duplicated proteins.
Duke.ESBL.remaining.duplicate.proteins <- Duke.ESBL.duplicate.proteins %>%
    filter(!str_detect(.$product,antibiotic.keywords)) %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    filter(!str_detect(.$product,unknown.protein.keywords))
sum(Duke.ESBL.remaining.duplicate.proteins$count) ## 821 other duplicate proteins in total.

################################################
###### Now, let's look at the singleton proteins.

## 558 singleton ARGs in the genomes.
Duke.ESBL.singleton.ARGs <- Duke.ESBL.singleton.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords))


## 3181 singleton MGE protein sequences in the genomes.
Duke.ESBL.singleton.MGE.proteins <- Duke.ESBL.singleton.proteins %>%
    filter(str_detect(.$product,IS.keywords))


## 6475 unknown singleton protein sequences in the genomes.
Duke.ESBL.singleton.unknown.proteins <- Duke.ESBL.singleton.proteins %>%
    ## some hypothetical proteins are "ISXX family insertion sequence hypothetical protein"
    ## so filter out those cases.
    filter(!str_detect(.$product,IS.keywords)) %>%
    filter(str_detect(.$product,unknown.protein.keywords))

## 47049 remaining cases of singleton proteins in the genomes.
Duke.ESBL.remaining.singleton.proteins <- Duke.ESBL.singleton.proteins %>%
    filter(!str_detect(.$product,antibiotic.keywords)) %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    filter(!str_detect(.$product,unknown.protein.keywords))

#####################################
## Now make Figure 8.

Fig8A.data <- Duke.ESBL.duplicate.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation_Accession, Category) %>%
    summarize(Count = sum(count))

Fig8B.data <- Duke.ESBL.singleton.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation_Accession, Category) %>%
    summarize(Count = sum(count))

Fig8A1 <- ggplot(Fig8A.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(legend.position="bottom") +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

Fig8legend <- get_legend(Fig8A1)
Fig8A1 <- Fig8A1 + guides(fill = FALSE)

Fig8A2 <- ggplot(Fig8A.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
    geom_bar(stat="identity", position = "fill") +
    theme_classic() +
    ## remove genome name labels.
    theme(axis.text.y=element_blank()) +
    guides(fill = FALSE) +
    xlab("Frequency") +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

Fig8A.title <- ggdraw() + draw_label("Distribution of duplicated genes in 12 ESBL-resistant isolates", fontface='bold')

Fig8A <- plot_grid(Fig8A.title,
                   plot_grid(Fig8A1, Fig8A2, labels = c("A",""),
                             nrow=1, rel_widths=c(2,1)),
                   nrow=2, rel_heights=c(0.2,2))


Fig8B1 <- ggplot(Fig8B.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
    geom_bar(stat="identity") +
    theme_classic() +
    guides(fill = FALSE) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

Fig8B2 <- ggplot(Fig8B.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
    geom_bar(stat="identity", position = "fill") +
    theme_classic() +
    ## remove genome name labels.
    theme(axis.text.y=element_blank()) +
    guides(fill = FALSE) +
    xlab("Frequency") +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.


Fig8B.title <- ggdraw() + draw_label("Distribution of single-copy genes in 12 ESBL-resistant isolates", fontface='bold')

Fig8B <- plot_grid(Fig8B.title,
                   plot_grid(Fig8B1, Fig8B2, labels = c("B",""),
                             nrow=1, rel_widths=c(2,1)),
                   nrow=2, rel_heights=c(0.2,2))

## now make the full Figure 8.
Fig8 <- plot_grid(Fig8A, Fig8B, Fig8legend, ncol = 1, rel_heights = c(2,2,0.25))
ggsave("../results/Fig8.pdf", Fig8, height = 7, width = 10)


################################################################################
## Analysis of chains of duplications, produced by join-duplications.py.
## Look at basic statistics
## for duplicated ARGs and associations with MGE genes,
## and re-calculate statistics for duplicated ARGs associated with MGEs,
## and duplicated ARGs that are not associated with MGEs.

joined.duplications <- read.csv("../results/joined-duplicate-proteins.csv") %>%
    tibble() %>%
    ## for numeric consistency, remove all duplications with NA product annotations.
    filter(!is.na(product))


make.ARG.MGE.region.contingency.table <- function(joined.duplications,
                                                  antibiotic.keywords,
                                                  IS.keywords) {

    ARG.joined.duplications <- joined.duplications %>%
        filter(str_detect(.$product,antibiotic.keywords))

    MGE.joined.duplications <- joined.duplications %>%
        filter(str_detect(.$product, IS.keywords))
    
    ## get the regions-- drop the sequence information.
    ## There are 756,165 regions in total.
    joined.regions <- joined.duplications %>%
        select(Annotation_Accession, Replicon_Accession, Replicon_type, region_index) %>%
        distinct()
    
    ARG.joined.regions <- ARG.joined.duplications %>%
        select(Annotation_Accession, Replicon_Accession, Replicon_type, region_index) %>%
        distinct()
    
    no.ARG.joined.regions <- anti_join(joined.regions, ARG.joined.regions)
    
    MGE.joined.regions <- MGE.joined.duplications %>%
        select(Annotation_Accession, Replicon_Accession, Replicon_type, region_index) %>%
        distinct()
    
    no.MGE.joined.regions <- anti_join(joined.regions, MGE.joined.regions)
    
    ## count regions (groups) that contain both ARGs and MGE genes.
    ## 3166 regions contain both ARGs and MGE genes.
    ARG.and.MGE.joined.regions <- inner_join(ARG.joined.regions, MGE.joined.regions)
    
    ## count regions (groups) that contain ARGs but no MGE genes.
    ## 2710 regions contain ARGs but no MGE genes.
    ARG.and.no.MGE.joined.regions <- inner_join(ARG.joined.regions, no.MGE.joined.regions)
    
    ## count regions (groups) that contain MGE genes but no ARGs.
    ## 572375 regions contain MGE genes but no ARGs.
    MGE.and.no.ARG.joined.regions <- inner_join(no.ARG.joined.regions, MGE.joined.regions)
    
    ## count regions (groups) that have neither MGE genes nor ARGs.
    ## 177363 regions contain neither MGE genes nor ARGs.
    no.MGE.and.no.ARG.joined.regions <- inner_join(no.ARG.joined.regions,no.MGE.joined.regions)
    
    ## now use a contingency table to test whether duplicated ARGs and duplicated MGEs
    ## are associated.
    joined.regions.contingency.table <- matrix(c(nrow(ARG.and.MGE.joined.regions),
                                                 nrow(ARG.and.no.MGE.joined.regions),
                                                 nrow(MGE.and.no.ARG.joined.regions),
                                                 nrow(no.MGE.and.no.ARG.joined.regions)),
                                               nrow = 2,
                                               dimnames = list(hasMGE = c("Yes","No"),
                                                               hasARG = c("Yes","No")))
    return(joined.regions.contingency.table)
}

joined.regions.contingency.table <- make.ARG.MGE.region.contingency.table(joined.duplications, antibiotic.keywords, IS.keywords)
fisher.test(joined.regions.contingency.table)

## the anti-correlation between duplicated MGEs and duplicated ARGs is stronger
## when just examining isolates from humans and livestock.

human.isolates <- gbk.annotation %>%
    filter(Annotation == "Human-host")

livestock.isolates <- gbk.annotation %>%
    filter(Annotation == "Livestock")

human.joined.duplications <- joined.duplications %>%
    filter(Annotation_Accession %in% human.isolates$Annotation_Accession)

livestock.joined.duplications <- joined.duplications %>%
    filter(Annotation_Accession %in% livestock.isolates$Annotation_Accession)

human.joined.regions.table <- make.ARG.MGE.region.contingency.table(human.joined.duplications, antibiotic.keywords, IS.keywords)

livestock.joined.regions.table <- make.ARG.MGE.region.contingency.table(livestock.joined.duplications, antibiotic.keywords, IS.keywords)


## Control analysis:
## Let's repeat this analysis, removing all duplicated regions that only contain MGE
## or unknown genes,
## to control for the possibility that this result is driven by the transposition of MGEs
## that don't have any passenger genes.
## This actually seems like the case.
## My feeling now is that it is probably better just to report the numbers for
## associations between ARGs and MGEs in the text, rather than making a figure
## about this point.

only.MGE.or.unknown.joined.regions <- joined.duplications %>%
    filter(str_detect(.$product, MGE.or.unknown.protein.keywords)) %>%
    group_by(Annotation_Accession, Replicon_Accession, Replicon_type, region_index) %>%
    mutate(num.MGE.genes = n()) %>%
    ## get the regions-- drop the sequence information.
    select(Annotation_Accession, Replicon_Accession, Replicon_type,
           region_index, num_proteins_in_region, num.MGE.genes) %>%
    distinct() %>%
    filter(num.MGE.genes == num_proteins_in_region) %>%
    mutate(MGE.region.filter.key = paste0(Annotation_Accession,Replicon_Accession,region_index))


not.just.MGE.joined.duplications <- joined.duplications %>%
    mutate(MGE.region.filter.key = paste0(Annotation_Accession,Replicon_Accession,region_index)) %>%
    filter(!(MGE.region.filter.key %in% only.MGE.or.unknown.joined.regions$MGE.region.filter.key))
write.csv("../results/no-only-MGE-or-unknown-joined-duplications.csv",
          x=not.just.MGE.joined.duplications)

not.just.MGE.joined.regions.contingency.table <- make.ARG.MGE.region.contingency.table(
    not.just.MGE.joined.duplications, antibiotic.keywords, IS.keywords)
fisher.test(not.just.MGE.joined.regions.contingency.table)


human.not.just.MGE.joined.duplications <- not.just.MGE.joined.duplications %>%
    filter(Annotation_Accession %in% human.isolates$Annotation_Accession)

livestock.not.just.MGE.joined.duplications <- not.just.MGE.joined.duplications %>%
    filter(Annotation_Accession %in% livestock.isolates$Annotation_Accession)

human.not.just.MGE.joined.regions.table <- make.ARG.MGE.region.contingency.table(
    human.not.just.MGE.joined.duplications, antibiotic.keywords, IS.keywords)

livestock.not.just.MGE.joined.regions.table <- make.ARG.MGE.region.contingency.table(
    livestock.not.just.MGE.joined.duplications, antibiotic.keywords, IS.keywords)


## let's make a figure in terms of percentages to make a simpler visualization of this pattern.
