## Aim1-analysis.R by Rohan Maddamsetti.

## analyse the distribution of AR genes on chromosomes versus plasmids in
## fully-sequenced genomes and plasmids in the NCBI Nucleotide database.

## TODO: add a scale to Figure 1C to interpret the size of the circles
## (percentage of duplicated ARGs on plasmids)


## CRITICAL ANALYSIS TODO: look for evidence of recent diversification.
## In particular look more deeply at the result that AAA+ ATPases,
## and ATPases in general seem to be enriched in gene duplications.

## CRITICAL ANALYSIS TODO: Make sure numbers in genome.database,
## gbk.annotation, and all.proteins, duplicate.proteins, and
## singleton.proteins, in terms of number of isolates in each
## category, are COMPLETELY consistent with each other.

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)

## annotate source sequence as plasmid or chromosome.
genome.database <- read.csv("../results/AR-gene-duplication/chromosome-plasmid-table.csv")

## CRITICAL TODO: edit annotate-ecological-category.py to take care of the "blank" entries.
gbk.annotation <- as_tibble(read.csv("../results/AR-gene-duplication/computationally-annotated-gbk-annotation-table.csv")) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## get species name annotation from genome.database
    left_join(genome.database)


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

## import the 15GB file containing all proteins, including singletons.
## I can save a ton of memory if I don't import the sequence column,
## and by using the data.table package for import.
all.proteins <- data.table::fread("../results/AR-gene-duplication/all-proteins.csv",
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
all.duplicate.proteins <- read.csv("../results/AR-gene-duplication/duplicate-proteins.csv") %>% left_join(gbk.annotation)
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
rm(all.duplicate.proteins)
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
write.csv(x=TableS1, file="../results/AR-gene-duplication/TableS1.csv")

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
write.csv(x=ControlTable1, file="../results/AR-gene-duplication/ControlTable1.csv")

gc() ## free memory after dealing with singleton data.
###################################
###################################
## Tables to plot percentage of genes on plasmids, for Figure 1C.

## Table S2. Show number of duplicated genes on chromosomes, and number of
## duplicate genes on plasmids, for each category, for duplicated genes
## and duplicated AR genes. This table of raw data goes into the text. Then,
## sum over all categories for a 2x2 contingency table and report the result of a
## Fisher's exact test for asssociation between duplicated AR genes and plasmids.

make.TableS2 <- function(duplicate.proteins) {
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

TableS2 <- make.TableS2(duplicate.proteins)
## write Table S2 to file.
write.csv(x=TableS2,file="../results/AR-gene-duplication/TableS2.csv")

################
## Analysis of Table S2: Duplicate ARGs are associated with plasmids.

plasmid.chromosome.duplicate.ARG.contingency.test <- function(TableS2) {
    ## get values for Fisher's exact test.
    total.chr.AR.duplicates <- sum(TableS2$chromosomal_duplicate_ARGs)
    total.plasmid.AR.duplicates <- sum(TableS2$plasmid_duplicate_ARGs)

    total.chr.duplicates <- sum(TableS2$chromosomal_duplicate_genes)
    total.plasmid.duplicates <- sum(TableS2$plasmid_duplicate_genes)

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

plasmid.chromosome.duplicate.ARG.contingency.test(TableS2)

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

make.S1Fig <- function(ControlTable2,TableS2) {

    S1Fig.data <- full_join(ControlTable2, TableS2) %>%
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

S1Fig <- make.S1Fig(ControlTable2,TableS2)
ggsave("../results/AR-gene-duplication/S1Fig.pdf",S1Fig,width=9,height=11)

##################################################################################
## Figure 1 A & B: Diagram of the analysis workflow, made in Inkscape/Illustrator.
##################################################################################
## Figure 1C: main figure, showing enrichment of AR duplicates
## in human hosts and livestock.

make.Fig1C.df <- function(TableS1, ControlTable1, TableS2, ControlTable2, duplicate.gene.control.df) {
    ## join duplicate and singleton tables to make Fig 1C.

    ## have to remove p-values from the two tables, because
    ## the column names are the same, but the values are different
    ## (because these are two different tests).
    no.pval.TableS1 <- select(TableS1, -corrected.pval)
    no.pval.ControlTable1 <- select(ControlTable1, -corrected.pval)
    
    Fig1C.df <- no.pval.TableS1 %>%
        full_join(no.pval.ControlTable1) %>%
        full_join(TableS2) %>%
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
        mutate(plasmid_AR_duplicate_percent = plasmid_duplicate_ARGs/(plasmid_duplicate_ARGs + chromosomal_duplicate_ARGs + 0.1))

    return(Fig1C.df)
}

Fig1C.df <- make.Fig1C.df(TableS1, ControlTable1, TableS2, ControlTable2,duplicate.gene.control.df)

make.Fig1C <- function(Fig1C.df) {

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
        annotate("text", x = 65, y = 1.6, label = "Isolates with duplicated ARGs",
                 angle = 34, color = Fig1C.color.palette[1],size=3) +
        annotate("text", x = 65, y = 13.5, label = "Isolates with singleton ARGs",
                 angle = 34, color = Fig1C.color.palette[2],size=3) +
        annotate("text", x = 65, y = 22, label = "Isolates with duplicated genes",
                 angle = 34, color = "gray",size=3) +
        guides(size=FALSE) 
        
    return(Fig1C)
}

Fig1C <- make.Fig1C(Fig1C.df)
ggsave(Fig1C,file="../results/AR-gene-duplication/Fig1C.pdf",width=4,height=4)

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
ggsave(FigS2,file="../results/AR-gene-duplication/FigS2.pdf",width=6,height=6)

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

ggsave("../results/AR-gene-duplication/DI-index.pdf", DI.ARGs.to.all.fig)

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
write.csv(x=TableS2,file="../results/AR-gene-duplication/TableS2.csv")

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

write.csv(x=TableS3, file="../results/AR-gene-duplication/TableS3.csv")
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
    distinct()

HGT.candidates <- all.duplicate.proteins %>%
    group_by(sequence) %>%
    summarize(number.of.genomes = n(),
              total.copies = sum(count),
              chromosome.copies = sum(chromosome_count),
              plasmid.copies = sum(plasmid_count)) %>%
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

plasmid.only.HGT.candidates <- HGT.candidate.summary %>%
    filter(chromosome.copies == 0)

non.MGE.plasmid.only.HGT.candidates <- plasmid.only.HGT.candidates %>%
    filter(!str_detect(.$product,IS.keywords))

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
    filter(Annotation == manual.annot.string) %>%
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

## this function filters duplicate proteins by annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
make.seq.freq.table <- function(annot.string) {
    duplicate.proteins %>%
    filter(Annotation == annot.string) %>%
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
make.plas.annotation.freq.table <- function(annot.string) {
    duplicate.proteins %>%
        filter(plasmid_count >= 1) %>%
        filter(Annotation == annot.string) %>%
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

## this function filters duplicate proteins by annotation category,
## supplied as an argument, and summarizes by count of product annotation strings.
make.on.plas.seq.freq.table <- function(annot.string) {
    duplicate.proteins %>%
        filter(plasmid_count >= 1) %>%
        filter(Annotation == annot.string) %>%
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
