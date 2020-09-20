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
genome.database <- read.csv("../results/chromosome-plasmid-table.csv")
## This is 7046 strains.
length(unique(genome.database$Annotation_Accession))

##  while this is 7047 strains.
gbk.annotation <- as_tibble(read.csv("../data/manually-curated-gbk-annotation-table.csv")) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Manual_Annotation = replace_na(Manual_Annotation,"Unannotated"))

########################
## TEMPORARY HACK FOR SELF-CONSISTENCY:
gbk.annotation <- gbk.annotation %>% filter(Annotation_Accession %in% genome.database$Annotation_Accession)
########################
## IMPORTANT NOTE: This should be ~7300 genomes after updating the gbk annotation.
## not sure if chromosome-plasmid-table will also be larger.
length(unique(gbk.annotation$Annotation_Accession))
## there's one in gbk.annotation that is missnig from genome.database:
## GCA_900492195.1_T2.26MG-112.21_plasmid
## don't worry about it for now: it won't affect the analysis anyway.
cds.counts <- read.csv("../results/protein_db_CDS_counts.csv")

## number of host and isolation_source annotations in gbk_annotation.
## 647 unique host annotations.
length(unique(gbk.annotation$host))
## 1763 unique isolation_source annotations.
length(unique(gbk.annotation$isolation_source))

protein.db.metadata <- genome.database %>%
    left_join(gbk.annotation) %>%
    left_join(cds.counts)    

## We have 7046 isolates in the AR.results dataframe.
## This is an invariant that will help with debugging.
length(unique(protein.db.metadata$Annotation_Accession))

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


## import the 12GB file containing singletons.
## I can save a ton of memory if I don't import the sequence column,
## and by using the data.table package for import.
all.proteins <- data.table::fread("../results/duplicate-proteins.csv",
                                  drop="sequence") %>%
    ## now merge with gbk annotation.
    ## I am doing a left_join here, because I want the NA Manual_Accessions
    ## in order to predict where these unannotated strains come from.
    left_join(gbk.annotation) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Manual_Annotation = replace_na(Manual_Annotation,"Unannotated"))

duplicate.proteins <- all.proteins %>% filter(count > 1)
singleton.proteins <- all.proteins %>% filter(count == 1)

## free up memory by deallocating all.proteins,
rm(all.proteins)
## and running garbage collection.
gc()

#################
## TEMPORARY HACK FOR SELF-CONSISTENCY:
duplicate.proteins <- duplicate.proteins %>%
    filter(Annotation_Accession %in% genome.database$Annotation_Accession) %>%
    filter(Annotation_Accession %in% gbk.annotation$Annotation_Accession)

singleton.proteins <- singleton.proteins %>%
    filter(Annotation_Accession %in% genome.database$Annotation_Accession) %>%
    filter(Annotation_Accession %in% gbk.annotation$Annotation_Accession)

##################    
##########################################
## CRITICAL BUG TO FIX:
## there are 7,369 isolates which is more than the annotated genomes.
## in part, this discrepancy has to do with the manual annotation having been
## constructed on an older version of the gbk_annotation.csv, which was
## missing some strains. There are probably other bugs here to fix as well.
## DEBUG THESE ERRORS AND MAKE THESE NUMBERS CONSISTENT!!!
## for now, I have restricted these data by filtering on Annotation_Accession--
## see the temporary fixes in the code above.

problem.data <- duplicate.proteins %>%
    filter(!(Annotation_Accession %in% genome.database$Annotation_Accession)) %>%
    select(-count,-chromosome_count,-plasmid_count,-product) %>%
    distinct()
## there are 871 isolates with Annotation_Accession, but no metadata at all? Why?
length(unique(problem.data$Annotation_Accession))

length(unique(IS.removed.from.duplicate.proteins$Annotation_Accession))
## 5,772 isolate after filtered out IS. This includes Unannotated isolates.

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
    select(-count,-chromosome_count,-plasmid_count,-product) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(isolates_with_duplicate_genes = n()) %>%
    arrange(desc(isolates_with_duplicate_genes))

## Third col: count the number of isolates with duplicated AR genes in each category.
AR.category.counts <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(-count,-chromosome_count,-plasmid_count,-product) %>%
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

###########################################################################
## Positive control 1: Make a version of Table 1, examining the distribution
## of AR genes that have NOT duplicated.

## Second column:
## count the number of isolates with singleton AR genes in each category.

AR.singleton.category.counts <- singleton.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(-count,-chromosome_count,-plasmid_count,-product) %>%
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

## It is likely that I will put both of these tables into the Supplement
## for greatest transparency.

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
ControlTable2A <- raw.ControlTable2 %>%
    filter(Manual_Annotation != "Unannotated") %>%
    calc.expected.AR.singletons() %>% ## fourth column
    calc.AR.singleton.enrichment.pvals() ## fifth column

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

## p < 1e-150
fisher.test(Table4.contingency.table)
fisher.test(Table4.contingency.table)$p.value

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
## 4 Analyses based on 4 sets of query sequences:
## 1) genes that Yi gave me.

Yi.queries.AR.hits.E50 <- read.csv("../results/test-run-phmmer-hits-E50.csv") %>%
    arrange(desc(hits))
Yi.queries.AR.hits.E20 <- read.csv("../results/test-run-phmmer-hits-E20.csv")
## each row is a replicon (i.e. chromosome, plasmid).

## 2) Resfam database.
Resfam.queries.AR.hits.E50 <- read.csv("../results/Resfam-hmmsearch-hits.csv")

## 3) CARD database.
CARD.queries.AR.hits.E50 <- read.csv("../results/CARD-is_a-ontology-36009-phmmer-hits.csv")

## TODO: 4) Database curated by Carolyn Zhang.
################################################################################
## Analysis 1: genes that Yi gave me.

## Parameters:
## AR.hits: parsed HMMER results.
## outf.prefix: prefix for Figures and Tables written to disk.
## protein.db.metadata: dataframe containing metadata for all chromosomes
## and plasmids that were searched for hits.
do.aim1.analysis <- function(AR.hits, outf.prefix,  protein.db.metadata) {
    print("********************")
    print("analysis beginning.")
    AR.results <- protein.db.metadata %>%
        left_join(AR.hits) %>%
        mutate(hits = ifelse(is.na(hits),0,hits)) %>%
        ## this filter makes sure that all analyzed samples have
        ## manual annotation.
        filter(!is.na(Manual_Annotation)) %>% 
        filter(SequenceType %in% c("chromosome","plasmid")) %>% ## skip weird sequence types
        arrange(desc(hits))
    
    chromosome.df <- AR.results %>%
        filter(SequenceType == "chromosome") %>%
        group_by(Organism, Strain, Annotation_Accession, host, isolation_source, Manual_Annotation) %>%
        summarize(chromosome.hits = sum(hits), chromosome.CDS=sum(CDS_count)) %>%
        arrange(desc(chromosome.hits))
    
    ## keep strain identity information in the analysis,
    ## and compare difference in number of hits on chromosome or plasmid.
    plasmid.df <- AR.results %>%
        filter(SequenceType == "plasmid") %>%
        group_by(Organism, Strain, Annotation_Accession, host, isolation_source, Manual_Annotation) %>%
        summarize(plasmid.hits = sum(hits), plasmid.CDS=sum(CDS_count)) %>%
        arrange(desc(plasmid.hits))
    ## some plasmids in these data do not have matching chromosomes.
    
    ## for main analysis, want genomes with both chromosomes and plasmids.
    plasmids.that.match.chromosomes <- plasmid.df %>%
        filter(Annotation_Accession %in% chromosome.df$Annotation_Accession)
    
    combined.df <- full_join(plasmids.that.match.chromosomes, chromosome.df) %>%
        mutate(HumanHost = ifelse(Manual_Annotation == "Human-host",1,0)) %>%
        mutate(Soil = ifelse(Manual_Annotation == "Soil",1,0)) %>%
        filter(HumanHost | Soil) %>%
        mutate(Treatment = HumanHost - Soil) %>%
        ## next: can't be both HumanHost and Soil, or neither HumanHost or Soil.
        filter(Treatment %in% c(1,-1)) %>% 
        mutate(Treatment = ifelse(Treatment>0,"Human-host", "Soil")) %>%
        select(-HumanHost,-Soil) %>% ## drop temporary columns
        mutate(plasmid.hits = ifelse(is.na(plasmid.hits),0,plasmid.hits)) %>%
        mutate(chromosome.hits = ifelse(is.na(chromosome.hits),0,chromosome.hits)) %>%
        mutate(total.hits = plasmid.hits + chromosome.hits) %>%
        arrange(desc(total.hits))
    
    ## Question 1: among all isolates, what fraction have at least one AR locus?
    print("Fraction of isolates with at least one antibiotic resistance locus:")
    print(nrow(filter(combined.df,total.hits>0))/nrow(combined.df))
    
    resistant.combined.df <- combined.df %>%
        filter(total.hits > 0) %>%
        mutate(percent.on.plasmid = plasmid.hits/total.hits) %>%
        mutate(percent.on.chromosome = chromosome.hits/total.hits) %>%
        mutate(total.CDS = plasmid.CDS + chromosome.CDS) %>%
        mutate(percent.plasmid.CDS = plasmid.CDS/total.CDS) %>%
        mutate(percent.chromosome.CDS = chromosome.CDS/total.CDS) %>%
        ## observed fraction on plasmid - expected fraction on plasmid. 
        mutate(obs.minus.expect.frac.on.plasmid = percent.on.plasmid - percent.plasmid.CDS) %>%
        mutate(obs.minus.expect.frac.on.chromosome = percent.on.chromosome - percent.chromosome.CDS) 
    
    ## Question 2: Of those that have hits, what is the distribution of hits on
    ## chromosome or plasmid?
    AR.distribution.plot <- ggplot(
        resistant.combined.df,
        aes(x=plasmid.hits,
            y=chromosome.hits,
            color=Treatment)) +
        geom_jitter(width = 0.4) +
        theme_classic() +
        guides(color=FALSE)
    AR.distr.plot.name = paste0("../results/",outf.prefix,"-AR-distribution.pdf")
    ggsave(AR.distr.plot.name, AR.distribution.plot, width=6, height=6)
    
    
    human.host.combined.df <- resistant.combined.df %>%
        filter(Treatment == "Human-host")
    
    no.host.combined.df <- resistant.combined.df %>%
        filter(Treatment == "Soil")
    
    print("Mean fraction of AR genes on plasmids, in human host isolates:")
    print(mean(human.host.combined.df$percent.on.plasmid))
    ## mean plasmid percentage in human host isolates: 6.14%
    print("Mean fraction of AR genes on plasmids, in isolates without a host:")
    print(mean(no.host.combined.df$percent.on.plasmid))
    
    ## the difference in percentage on plasmid is NOT significantly different:
    ## Mann-Whitney U-test: p = 0.06082
    print("Mann-Whitney U-test 1: do isolates from human hosts have a larger fraction of AR genes on plasmids?")
    test.result1 <- wilcox.test(human.host.combined.df$percent.on.plasmid,
                                no.host.combined.df$percent.on.plasmid,
                                alternative="greater")
    print(test.result1)
    
    percentFig <- ggplot(resistant.combined.df, aes(x=Treatment,y=percent.on.plasmid)) +
        geom_boxplot() +
        geom_jitter() +
        theme_classic() +
        ylab("Percent of AR loci on plasmids")
    percentFig.name = paste0("../results/",outf.prefix,"-percentFig.pdf")
    ggsave(percentFig.name, percentFig, width=4, height=4)

    print("analysis finished.")
    return(resistant.combined.df)
}

## Analysis 1: genes that Yi gave me.

## at a -log10(E-value) threshold of 20, most hits are in soil bacteria,
## and most hits are chromosomally encoded.
## at a -log10(E-value) threshold of 50, 2 bacteria isolated perianally have
## plasmids that are highly enriched in antibiotic resistance genes, and >80%
## of resistance genes in the strain are on the plasmid (9/11 and 8/10).
resistance.df1 <- do.aim1.analysis(Yi.queries.AR.hits.E20, "Yi-testE20",  protein.db.metadata)
resistance.df2 <- do.aim1.analysis(Yi.queries.AR.hits.E50, "Yi-testE50",  protein.db.metadata)

Resfam.resistance.df <- do.aim1.analysis(Resfam.queries.AR.hits.E50,"ResfamE50", protein.db.metadata)
data.frame(Resfam.resistance.df) %>% arrange(desc(plasmid.hits))
## POSSIBLE BUG: NUMBERS ARE SUSPICIOUS IN THESE TOP HITS!
## CHECK TO SEE IF THERE IS SOMETHING WRONG HERE!
## see identical hit number and CDS numbers for Paraburkholderia caribensis.


CARD.resistance.df <- do.aim1.analysis(CARD.queries.AR.hits.E50,"CARD-E50", protein.db.metadata)
head(data.frame(CARD.resistance.df))
## POSSIBLE BUG: NUMBERS ARE SUSPICIOUS IN THESE TOP HITS!
## CHECK TO SEE IF THERE IS SOMETHING WRONG HERE!
## see identical hit number and CDS numbers for Paraburkholderia caribensis.


## My guess is that there will be far more chromosomally-encoded resistance in the
## complete set of Resfam and CARD datasets queries, because of the many varied
## mechanisms of resistance.

## What kinds of antibiotic resistance mechanisms do we expect to be enriched on plasmids?

## or rather, are any classes of antibiotic resistance genes enriched on plasmids in these
## data? Based on the CARD ontology?

## POTENTIAL TODO: Reduce Resfams dataset on antibiotic resistance genes that have
## known associations with mobile elements/transposons in the literature.


## Two choices here, to avoid circularity and faulty inference:
## 1) ask whether particular classes of resistance genes are enriched a priori (ask Yi).
## 2) describe the classes of AR resistance that best fit plasmid enrichment in the genome
##data.

## Perhaps I could then focus metagenomic analysis on those classes of genes,
## and ask whether they are encoded on plasmids in the metagenomic data.
