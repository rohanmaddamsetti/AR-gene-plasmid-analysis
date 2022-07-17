## ARG-duplication-analysis.R by Rohan Maddamsetti.
## analyse the distribution of antibiotic resistance genes (ARGs)
## on chromosomes versus plasmids in  fully-sequenced genomes and plasmids
## in the NCBI Nucleotide database.

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)


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

## read in duplicate proteins with sequences, using a separate file.
## I want the sequence column for the duplicate genes,
## but not for the singletons, to save memory.
duplicate.proteins <- read.csv("../results/duplicate-proteins.csv") %>%
    ## now merge with gbk annotation.
    inner_join(gbk.annotation)

######## Lingchong asked for this control analysis.
## by default, don't count plasmid proteins as duplicates.
## The results are robust to this assumption; nothing changes.
COUNT.PLASMID.PROTEINS.AS.DUPLICATES <- FALSE

if (COUNT.PLASMID.PROTEINS.AS.DUPLICATES) {

    ## CRITICAL STEP: get plasmid proteins to count as duplicates.
    plasmid.proteins <- all.proteins %>%
        filter(plasmid_count >= 1) %>%
        inner_join(gbk.annotation)
    
    ## CRITICAL STEP: join plasmid proteins as duplicates.
    duplicate.proteins <- duplicate.proteins %>%
        left_join(plasmid.proteins)
    ## remove plasmid.proteins from memory once we are done with it.
    rm(plasmid.proteins)
       ## now get the singleton protein by filtering.
    singleton.proteins <- all.proteins %>%
        filter(count == 1) %>%
        ## proteins on plasmids do not count as singletons in this analysis.
        filter(plasmid_count == 0) %>%
        inner_join(gbk.annotation)

    } else { ## just get the singleton protein by filtering.
        singleton.proteins <- all.proteins %>%
            filter(count == 1) %>%
            inner_join(gbk.annotation)     
}
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

##########################################################################
## Figure 1 is a schematic of the analysis pipeline.
##########################################################################
## Code and data structures for Figure 2ABC.

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
        xlab("% of Isolates") +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=0.2, size=0.2)
    
    if (no.category.label)
        Fig.panel <- Fig.panel +
            theme(axis.text.y=element_blank())
    
    return(Fig.panel)
}

## Data structure for Figure 2A:
## normal-approximation confidence intervals for the percentage
## of isolates with duplicated ARGs.
TableS1 <- make.TableS1(gbk.annotation, duplicate.proteins)
## write Supplementary Table S1 to file.
write.csv(x=TableS1, file="../results/TableS1.csv")

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

gc() ## free memory after dealing with singleton data.

#########################################################################
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

## Data structure for Figure 2C.
TableS3 <- make.TableS3(gbk.annotation, duplicate.proteins)
## write TableS3 to file.
write.csv(x=TableS3, file="../results/TableS3.csv")

######################################################################
## Supplementary Figure S2: Control for Taxonomy (simpler control than for phylogeny)

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

## Let's mark the genera with the most duplicated ARGs, and recalculate a version
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
S2FigA <- genera.isolate.comparison.df %>%
    ggplot(aes(x=sqrt(genome.count),
               y = sqrt(duplicated.ARG.genome.count),
               label = Genus,
               color = in.top.ARG.genera)) +
    theme_classic() + geom_jitter() + geom_text_repel(fontface = "italic") +
    scale_color_manual(values=c("black", "red")) +
    guides(color="none") +
    xlab("sqrt(Number of isolates)") +
    ylab("sqrt(Number of isolates with D-ARGs)")

## 7,195 isolates are in the top ARG genera.
top.ARG.genera.isolates <- gbk.annotation %>%
    filter(Genus %in% top.ARG.genera)

## 11,509 isolates are in the remaining genera.
filtered.gbk.annotation <- gbk.annotation %>%
    filter(!(Genus %in% top.ARG.genera))

filtered.duplicate.proteins <- duplicate.proteins %>%
    filter(!(Genus %in% top.ARG.genera))

filtered.TableS1 <- make.TableS1(filtered.gbk.annotation, filtered.duplicate.proteins)

S2FigB <- make.confint.figure.panel(filtered.TableS1, order.by.total.isolates, "D-ARGs after\nfiltering top genera")

## let's check the results, just for the top ARG genera.
top.ARG.genera.duplicate.proteins <- duplicate.proteins %>%
    filter(Genus %in% top.ARG.genera)

top.ARG.genera.TableS1 <- make.TableS1(top.ARG.genera.isolates,
                                       top.ARG.genera.duplicate.proteins)

S2FigC <- make.confint.figure.panel(
    top.ARG.genera.TableS1, order.by.total.isolates,
    "D-ARGs in\ntop genera only",
    no.category.label = TRUE)

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

S2FigD <- make.confint.figure.panel(
    single.organism.TableS1, order.by.total.isolates,
    "D-ARGs after\ndownsampling species",
        no.category.label = TRUE)

S2FigBCD <- plot_grid(S2FigB, S2FigC, S2FigD, labels = c("B","C",'D'), nrow = 1, rel_widths = c(1.5, 1, 1, 1))


S2Fig <- plot_grid(S2FigA, S2FigBCD, labels = c("A", ""), nrow = 2, rel_heights = c(1.3,1))

ggsave("../results/S2Fig.pdf", S2Fig, width=7.5)

rm(single.organism.duplicate.proteins)
rm(single.organism.singleton.proteins)
gc()

##################################################################
## Figures S1 and S3. Proportion of isolates with duplicated or single-copy ARGs
## for 12 different antibiotic classes.
########################################
## Figure S1: Duplicated ARGs.

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


S1FigA <- make.confint.figure.panel(chloramphenicol.table,
                                    order.by.total.isolates,
                                    "chloramphenicol\nresistance")
S1FigB <- make.confint.figure.panel(tetracycline.table,
                                    order.by.total.isolates,
                                    "tetracycline\nresistance",
                                    no.category.label = TRUE)
S1FigC <- make.confint.figure.panel(MLS.table,
                                    order.by.total.isolates,
                                    "MLS\nresistance",
                                    no.category.label = TRUE)
S1FigD <- make.confint.figure.panel(multidrug.table,
                                    order.by.total.isolates,
                                    "multidrug\nresistance")
S1FigE <- make.confint.figure.panel(beta.lactam.table,
                                    order.by.total.isolates,
                                    "beta-lactam\nresistance",
                                    no.category.label = TRUE)
S1FigF <- make.confint.figure.panel(glycopeptide.table,
                                    order.by.total.isolates,
                                    "glycopeptide\nresistance",
                                    no.category.label = TRUE)
S1FigG <- make.confint.figure.panel(polypeptide.table,
                                    order.by.total.isolates,
                                    "polypeptide\nresistance")
S1FigH <- make.confint.figure.panel(diaminopyrimidine.table,
                                    order.by.total.isolates,
                                    "diaminopyrimidine\nresistance",
                                    no.category.label = TRUE)
S1FigI <- make.confint.figure.panel(sulfonamide.table,
                                    order.by.total.isolates,
                                    "sulfonamide\nresistance",
                                    no.category.label = TRUE)
S1FigJ <- make.confint.figure.panel(quinolone.table,
                                    order.by.total.isolates,
                                    "quinolone\nresistance")
S1FigK <- make.confint.figure.panel(aminoglycoside.table,
                                    order.by.total.isolates,
                                    "aminoglycoside\nresistance",
                                    no.category.label = TRUE)
S1FigL <- make.confint.figure.panel(macrolide.table,
                                    order.by.total.isolates,
                                    "macrolide\nresistance",
                                    no.category.label = TRUE)

S1Fig <- plot_grid(NULL, ## The nesting is to add a title.
                   plot_grid(
                       S1FigA, S1FigB, S1FigC,
                       S1FigD, S1FigE, S1FigF,
                       S1FigG, S1FigH, S1FigI,
                       S1FigJ, S1FigK, S1FigL,
                   rel_widths = c(1.5, 1, 1,
                                  1.5, 1, 1,
                                  1.5, 1, 1),
                   nrow = 4),
                   labels = c("D-ARGs", ""),
                   ncol = 1,
                   rel_heights = c(0.025, 1))

ggsave("../results/S1Fig.pdf", S1Fig, height = 10.5, width = 8)
########################################
## Figure S3: Single-copy ARGs.

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


S3FigA <- make.confint.figure.panel(chloramphenicol.control.table,
                                      order.by.total.isolates,
                                      "chloramphenicol\nresistance")
S3FigB <- make.confint.figure.panel(tetracycline.control.table,
                                      order.by.total.isolates,
                                      "tetracycline\nresistance",
                                      no.category.label = TRUE)
S3FigC <- make.confint.figure.panel(MLS.control.table,
                                      order.by.total.isolates,
                                      "MLS\nresistance",
                                      no.category.label = TRUE)
S3FigD <- make.confint.figure.panel(multidrug.control.table,
                                      order.by.total.isolates,
                                      "multidrug\nresistance")
S3FigE <- make.confint.figure.panel(beta.lactam.control.table,
                                      order.by.total.isolates,
                                      "beta-lactam\nresistance",
                                      no.category.label = TRUE)
S3FigF <- make.confint.figure.panel(glycopeptide.control.table,
                                      order.by.total.isolates,
                                      "glycopeptide\nresistance",
                                      no.category.label = TRUE)
S3FigG <- make.confint.figure.panel(polypeptide.control.table,
                                      order.by.total.isolates,
                                      "polypeptide\nresistance")
S3FigH <- make.confint.figure.panel(diaminopyrimidine.control.table,
                                      order.by.total.isolates,
                                      "diaminopyrimidine\nresistance",
                                      no.category.label = TRUE)
S3FigI <- make.confint.figure.panel(sulfonamide.control.table,
                                      order.by.total.isolates,
                                      "sulfonamide\nresistance",
                                      no.category.label = TRUE)
S3FigJ <- make.confint.figure.panel(quinolone.control.table,
                                    order.by.total.isolates,
                                    "quinolone\nresistance")
S3FigK <- make.confint.figure.panel(aminoglycoside.control.table,
                                      order.by.total.isolates,
                                      "aminoglycoside\nresistance",
                                      no.category.label = TRUE)
S3FigL <- make.confint.figure.panel(macrolide.control.table,
                                      order.by.total.isolates,
                                      "macrolide\nresistance",
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
                   labels = c("S-ARGs", ""),
                   ncol = 1,
                   rel_heights = c(0.025, 1))
ggsave("../results/S3Fig.pdf", S3Fig, height = 10.5, width = 8)

#########################
## S5 Figure.
## Analysis of duplicate pairs found just on chromosome, just on plasmid, or
## on both chromosomes and plasmids.

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
S5FigA <- S5FigA + guides(fill = "none")

S5FigB <- ggplot(just.chromosome.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("Chromosome only") +
    guides(fill = "none") +
    ylab("")

S5FigC <- ggplot(just.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("Plasmid only") +
    guides(fill = "none") +
    ylab("")

S5Fig <- plot_grid(NULL, S5FigA, S5FigB, S5FigC, S5Fig.legend, ncol = 1,
                   labels = c("Genomic distribution\n of D-genes", "A","B","C"),
                   rel_heights=c(0.35,1,1,1,0.25))
ggsave("../results/S5Fig.pdf", S5Fig, width=4, height=8)

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
## Use the data in Tables S4 and S5 to make the data structures for Figure 2DEFGHI.
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

## Finally -- make Figure 2ABC.
Fig2A <- make.confint.figure.panel(TableS1, order.by.total.isolates, "D-ARGs") +
        scale_x_continuous(breaks = c(0, 0.15))

Fig2B <- make.confint.figure.panel(TableS2, order.by.total.isolates,
                                   "S-ARGs", no.category.label=TRUE) +
    scale_x_continuous(breaks = c(0.85, 1.0))

Fig2C <- make.confint.figure.panel(TableS3, order.by.total.isolates,
                                   "All D-genes", no.category.label=TRUE) +
    scale_x_continuous(breaks = c(0.75, 0.95))

Fig2ABC.title <- title_theme <- ggdraw() +
    draw_label("Isolate-level analysis",fontface="bold")

Fig2ABC <- plot_grid(Fig2A, Fig2B, Fig2C, labels=c('A','B','C'),
                     rel_widths = c(1.5,1,1), nrow=1)

Fig2ABC.with.title <- plot_grid(Fig2ABC.title, Fig2ABC, ncol = 1, rel_heights = c(0.1, 1))

## the rest of Figure 2 -- Fig2DEFGHI -- requires Tables S4 and S5.
## I manually set axis labels so that they don't run into each other.
Fig2D <- make.Fig2DEFGHI.panel(Fig2D.df, order.by.total.isolates,
                         "\nD-ARGs",
                         "% of\nall genes") +
    scale_x_continuous(label=fancy_scientific, breaks = c(0, 2e-4), limits = c(0,2.2e-4))
Fig2E <- make.Fig2DEFGHI.panel(Fig2E.df, order.by.total.isolates,
                         "Chromosome:\nD-ARGs",
                         "% of\nchromosomal genes",
                         no.category.label = TRUE) +
    scale_x_continuous(label=fancy_scientific, breaks = c(0, 8e-5), limits = c(0,9e-5))
Fig2F <- make.Fig2DEFGHI.panel(Fig2F.df, order.by.total.isolates,
                         "Plasmids:\nD-ARGs",
                         "% of\nplasmid genes",
                         no.category.label = TRUE) +
    scale_x_continuous(label=fancy_scientific, breaks = c(0, 5e-3), limits = c(0,5.4e-3))
Fig2G <- make.Fig2DEFGHI.panel(Fig2G.df, order.by.total.isolates,
                         "\nS-ARGs",
                         "% of\nall genes") +
    scale_x_continuous(label=fancy_scientific, breaks = c(3.5e-3, 7e-3), limits = c(3.5e-3,7.4e-3))
Fig2H <- make.Fig2DEFGHI.panel(Fig2H.df, order.by.total.isolates,
                         "Chromosome:\nS-ARGs",
                         "% of\nchromosomal genes",
                         no.category.label = TRUE) +
    scale_x_continuous(label=fancy_scientific, breaks = c(4e-3, 6e-3))
Fig2I <- make.Fig2DEFGHI.panel(Fig2I.df, order.by.total.isolates,
                         "Plasmids:\nS-ARGs",
                         "% of\nplasmid genes",
                         no.category.label = TRUE) +
    scale_x_continuous(label=fancy_scientific, breaks = c(3e-3, 2e-2))

Fig2DEFGHI.title <- title_theme <- ggdraw() +
    draw_label("Gene-level analysis",fontface="bold")

Fig2DEFGHI <- plot_grid(Fig2D, Fig2E, Fig2F, Fig2G, Fig2H, Fig2I,
                  labels = c('D','E','F','G','H','I'), nrow=2,
                  rel_widths = c(1.5, 1, 1, 1.5, 1, 1))

Fig2DEFGHI.with.title <- plot_grid(Fig2DEFGHI.title, Fig2DEFGHI,
                                   ncol = 1, rel_heights = c(0.06, 1))

## Now, make the complete Figure 2!
Fig2 <- plot_grid(Fig2ABC.with.title, Fig2DEFGHI.with.title,
                  ncol = 1, rel_heights = c(0.3,0.7))

ggsave("../results/Fig2.pdf", Fig2, height=8, width=5.6)

##########################################################################
## Figure 3: Visualization of ARGs on plasmids and chromosomes, and evidence for selection.

## set up data structures for Figure 3AB.
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
    xlab("% of D-genes") +
    scale_x_continuous(breaks = c(0,1)) +
    theme(legend.position="bottom") +
    theme(strip.background = element_blank()) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

Fig3legend <- get_legend(Fig3A)
Fig3A <- Fig3A + guides(fill = "none")

Fig3B <- ggplot(Fig3B.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    facet_wrap(.~Episome) +
    theme_classic() +
    xlab("% of S-genes") +
    scale_x_continuous(breaks = c(0,1)) +
    guides(fill = "none") +
    theme(strip.background = element_blank()) +
    theme(axis.text.y=element_blank()) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

## remove the dataframes to save memory.
rm(Fig3A.data)
rm(Fig3B.data)
gc()

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
    xlab("(% of D-genes) / (% of S-genes)") +
    ##xlab(TeX("\\frac{% of D-genes}{% of S-genes}")) +
    guides(color = "none") +
    theme(axis.text.y=element_blank()) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

mainFig3 <- plot_grid(Fig3A, Fig3B, Fig3C, labels = c('A','B','C'), nrow = 1, rel_widths=c(1.4,1,1),
                  align = 'h', axis = 'tb')

Fig3 <- plot_grid(mainFig3, Fig3legend, ncol = 1, rel_heights = c(1,0.2))

ggsave("../results/Fig3.pdf", Fig3, width=8.5, height=3.5)

################################################################################
## Figure 4A. A deterministic ODE model demonstrates that selection can
## drive the evolution of duplicated ARGs on plasmids.

## The panels of this figure are generated in my Pluto notebook:
## duplication-linear-ODE-model.jl.
################################################################################
## Supplementary Figure S4.
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
## Now make Supplementary Figure S4.

S4FigA.data <- Duke.ESBL.duplicate.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation_Accession, Category) %>%
    summarize(Count = sum(count))

S4FigB.data <- Duke.ESBL.singleton.proteins %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
    group_by(Annotation_Accession, Category) %>%
    summarize(Count = sum(count))

S4FigA1 <- ggplot(S4FigA.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(legend.position="bottom") +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

S4Figlegend <- get_legend(S4FigA1)
S4FigA1 <- S4FigA1 + guides(fill = "none")

S4FigA2 <- ggplot(S4FigA.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
    geom_bar(stat="identity", position = "fill") +
    theme_classic() +
    ## remove genome name labels.
    theme(axis.text.y=element_blank()) +
    guides(fill = "none") +
    xlab("Frequency") +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

S4FigA.title <- ggdraw() + draw_label("Distribution of duplicated genes in 12 ESBL-resistant isolates", fontface='bold')

S4FigA <- plot_grid(S4FigA.title,
                   plot_grid(S4FigA1, S4FigA2, labels = c("A",""),
                             nrow=1, rel_widths=c(2,1)),
                   nrow=2, rel_heights=c(0.2,2))


S4FigB1 <- ggplot(S4FigB.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
    geom_bar(stat="identity") +
    theme_classic() +
    guides(fill = "none") +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

S4FigB2 <- ggplot(S4FigB.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
    geom_bar(stat="identity", position = "fill") +
    theme_classic() +
    ## remove genome name labels.
    theme(axis.text.y=element_blank()) +
    guides(fill = "none") +
    xlab("Frequency") +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.


S4FigB.title <- ggdraw() + draw_label("Distribution of single-copy genes in 12 ESBL-resistant isolates", fontface='bold')

S4FigB <- plot_grid(S4FigB.title,
                   plot_grid(S4FigB1, S4FigB2, labels = c("B",""),
                             nrow=1, rel_widths=c(2,1)),
                   nrow=2, rel_heights=c(0.2,2))

## now make the full Supplementary Figure S4.
S4Fig <- plot_grid(S4FigA, S4FigB, S4Figlegend, ncol = 1, rel_heights = c(2,2,0.25))
ggsave("../results/S4Fig.pdf", S4Fig, height = 7, width = 10)

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

