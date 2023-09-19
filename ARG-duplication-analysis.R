## ARG-duplication-analysis.R by Rohan Maddamsetti.
## analyse the distribution of antibiotic resistance genes (ARGs)
## on chromosomes versus plasmids in  fully-sequenced genomes and plasmids
## in the NCBI RefSeq database (dated March 26 2021).

## IMPORTANT NOTE: Julia needs to be accessible from the shell as called by the R function system()
## near the end of this script. This is for a relatively minor result; the big stuff doesn't depend on Julia
## being accessible from R.

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)


################################################################################
## Global variables that affect assumptions going into the analysis.

## For control analysis to check robustness of results to method used to score ARGs
## and MGE-associated genes.
## By default set to FALSE.
USE.CARD.AND.MOBILE.OG.DB <- FALSE

## by default, don't count plasmid proteins as duplicates.
## The results are robust to this assumption; nothing changes.
COUNT.PLASMID.PROTEINS.AS.DUPLICATES <- FALSE
################################################################################

fancy_scientific <- function(x) {
    ## function for plotting better axis labels.
    ## see solution here for nice scientific notation on axes.
    ## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

################################################################################
## Regular expressions used in this analysis.

## We build on the regular expressions used by Zeevi et al. (2019).
## Transposon: ‘transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon’
## plasmid: ‘relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation’
## phage: ‘capsid|phage|tail|head|tape measure|antiterminatio’
## other HGT mechanisms: ‘integrase|excision\S*|exo- nuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip’
## antibiotic resistance: ‘azole resistance|antibiotic resistance|TetR|tetracycline resistance|VanZ|betalactam\S*|beta-lactam|antimicrob\S*|lantibio\S*’.


## unknown protein keywords.
unknown.protein.keywords <- "unknown|Unknown|hypothetical|Hypothetical|Uncharacterized|Uncharacterised|uncharacterized|uncharacterised|DUF|unknow|putative protein in bacteria|Unassigned|unassigned"

## NOTE: some hypothetical proteins are "ISXX family insertion sequence hypothetical protein"
## so filter out those cases, when counting unknown proteins.

## match MGE genes using the following keywords in the "product" annotation
transposon.keywords <- "IS|transpos\\S*|insertion|Tra[A-Z]|Tra[0-9]|tra[A-Z]|conjugate transposon|Transpos\\S*|Tn[0-9]|tranposase|Tnp|Ins|ins"
plasmid.keywords <- "relax\\S*|conjug\\S*|mob\\S*|plasmid|type IV|chromosome partitioning|chromosome segregation|Mob\\S*|Plasmid|Rep|Conjug\\S*"
phage.keywords <- "capsid|phage|Tail|tail|head|tape measure|antiterminatio|Phage|virus|Baseplate|baseplate|coat|entry exclusion"
other.HGT.keywords <- "Integrase|integrase|excision\\S*|exonuclease|recomb|toxin|restrict\\S*|resolv\\S*|topoisomerase|reverse transcrip|intron|antitoxin|toxin|Toxin|Reverse transcriptase|hok|Hok|competence|addiction"


MGE.keywords <- paste(transposon.keywords, plasmid.keywords, phage.keywords, other.HGT.keywords, sep="|")
MGE.or.unknown.protein.keywords <- paste(MGE.keywords,unknown.protein.keywords,sep="|")

## Elongation Factor Tu (2 copies in most bacteria).
## \\b is a word boundary.
## see: https://stackoverflow.com/questions/62430498/detecting-whole-words-using-str-detect-in-r
EFTu.keywords <- "\\bTu | Tu\\b|-Tu\\b"

## antibiotic-specific keywords.
chloramphenicol.keywords <- "chloramphenicol|Chloramphenicol"
tetracycline.keywords <- "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
MLS.keywords <- "macrolide|lincosamide|streptogramin"
multidrug.keywords <- "Multidrug resistance|multidrug resistance|antibiotic resistance"
beta.lactam.keywords <- "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\\S*"
glycopeptide.keywords <- "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
polypeptide.keywords <- "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
diaminopyrimidine.keywords <- "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
sulfonamide.keywords <- "sulfonamide|Sul1|sul1|sulphonamide"
quinolone.keywords <- "quinolone|Quinolone|oxacin|qnr|Qnr"
aminoglycoside.keywords <- "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
macrolide.keywords <- "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
antimicrobial.keywords <- "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\\S*"


antibiotic.keywords <- paste(chloramphenicol.keywords, tetracycline.keywords, MLS.keywords, multidrug.keywords,
    beta.lactam.keywords, glycopeptide.keywords, polypeptide.keywords, diaminopyrimidine.keywords,
    sulfonamide.keywords, quinolone.keywords, aminoglycoside.keywords, macrolide.keywords, antimicrobial.keywords, sep="|")

antibiotic.or.MGE.keywords <- paste(MGE.keywords,antibiotic.keywords,sep="|")


categorize.as.MGE.ARG.or.other <- function(product) {
    if (is.na(product))
        return("Other function")
    else if (str_detect(product, antibiotic.keywords))
        return("ARG")
    else if (str_detect(product, MGE.keywords))
        return("MGE")
    else
        return("Other function")
}


################################################################################
## Set up the key data structures for the analysis:
## gbk.annotation, in particular.

## import the 37GB file containing all proteins, including singletons.
## I can save a ton of memory if I don't import the sequence column,
## and by using the data.table package for import.
all.proteins <- data.table::fread("../results/all-proteins.csv",
                                  drop="sequence")

## get the genomes with chromosomes smaller than some plasmid in the genome.
## we will remove these 12 genomes from the analysis.
bad.replicons <- read.csv("../results/bad-replicons.csv")

## get the genomes that passed assembly quality control (QC),
## and anti_join with bad.replicons to add an additional layer of QC.
## We will use this to filter episome.database and gbk.annotation.
QCed.genomes <- read.csv("../results/genome-assembly-metadata.csv") %>%
    anti_join(bad.replicons) %>%
    as_tibble()

## annotate source sequences as plasmid or chromosome.
episome.database <- read.csv("../results/chromosome-plasmid-table.csv") %>%
    as_tibble() %>%
    ## filter based on QCed.genomes.
    filter(Annotation_Accession %in% QCed.genomes$Annotation_Accession)

## Double-check that every plasmid is associated with a chromosome.
check.plasmid.accessions <- function(episome.database) {
    plasmid.accessions <- sort(filter(episome.database, SequenceType=='plasmid')$Annotation_Accession)
    chromosome.accessions <- sort(filter(episome.database, SequenceType=='chromosome')$Annotation_Accession)
    bad.plasmid.vec <- sapply(plasmid.accessions,function(x) ifelse(x %in% chromosome.accessions,0,1))
    if (sum(bad.plasmid.vec) > 0)
        return(1) ## at least one plasmid is not associated with a chromosome
    else
        return(0) ## all plasmids are associated with chromosomes.
}
## run the test.
check.plasmid.accessions(episome.database)

## I used the script cross-check-Hawkey2022-accessions.py
## to find 8 genomes in the Hawkey et al. 2022 dataset that are present in the main dataset
## in gbk.annotation.
## I then found another 5 based on an assertion statement  to make sure none of these
## genomes are analyzed in the main dataset.
##Let's omit these from the main analysis to preserve independence of the
## clinical validation data.

Hawkey.overlaps.in.gbk.annotation.vec <- c(
    ## These were found with cross-chrck-Hawkey2022-accessions.py.
    "GCF_016126855.1_ASM1612685v1",
    "GCF_015476295.1_ASM1547629v1",
    "GCF_015325925.1_ASM1532592v1",
    "GCF_015999425.1_ASM1599942v1",
    "GCF_015999405.1_ASM1599940v1",
    "GCF_017584065.1_ASM1758406v1",
    "GCF_904864465.1_INF298",
    "GCF_904864595.1_INF333",
    ## I caught these with an assertion.
    "GCF_904866485.1_KSB1_1H",
    "GCF_904863335.1_KSB1_5B",
    "GCF_904866485.1_KSB1_1H",
    "GCF_904863335.1_KSB1_5B",
    "GCF_904863225.1_KSB1_6J",
    "GCF_904866475.1_KSB2_9B",
    "GCF_904863435.1_KSB1_1B")

gbk.annotation <- read.csv(
    "../results/computationally-annotated-gbk-annotation-table.csv") %>%
    as_tibble() %>%
    ## filter based on QCed.genomes.
    filter(Annotation_Accession %in% QCed.genomes$Annotation_Accession) %>%
    ## remove the overlaps with the Hawkey et al. (2022) clinical validation data.
    filter(!(Annotation_Accession %in% Hawkey.overlaps.in.gbk.annotation.vec)) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## collapse Annotations into a smaller number of categories as follows:
    ## Marine, Freshwater --> Water
    ## Sediment, Soil, Terrestrial --> Earth
    ## Plants, Agriculture, Animals --> Plants & Animals
    mutate(Annotation = replace(Annotation, Annotation == "Marine", "Water")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Freshwater", "Water")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Sediment", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Soil", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Terrestrial", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Plants", "Plants & Animals")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Agriculture", "Plants & Animals")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Animals", "Plants & Animals")) %>%
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
    distinct() %>%
    ## And now remove all Unannotated genomes, since these are not analyzed
    ## at all in this first paper.
    filter(Annotation != "Unannotated") %>%
    ## and remove any strains (although none should fall in this category)
    ## that were not annotated by annotate-ecological-category.py.
    filter(Annotation != "blank")


## print gbk.annotation to file for cross-checking.
write.csv(gbk.annotation, file="../results/gbk-annotation-of-analyzed-complete-genomes.csv",
          row.names=FALSE,quote=F)

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

## score ARGs and MGE-associated genes
## using the CARD and mobileOG-db databases (>80% ID over >85% length of best hit).
all.proteins.in.CARD <- read.csv("../results/all-proteins-in-CARD.csv") %>%
    as_tibble() %>%
    ## now merge with gbk annotation.
    inner_join(gbk.annotation)

all.proteins.in.mobileOGdb <- read.csv("../results/all-proteins-in-mobileOG-db.csv") %>%
    as_tibble() %>%
    ## now merge with gbk annotation.
    inner_join(gbk.annotation)


## read in duplicate proteins with sequences, using a separate file.
## I want the sequence column for the duplicate genes,
## but not for the singletons, to save memory.
duplicate.proteins <- read.csv("../results/duplicate-proteins.csv") %>%
    ## now merge with gbk annotation.
    inner_join(gbk.annotation)

######## Control analysis for separate project:
## print out a table of ARGs found on plasmids.
plasmid.ARGs <- all.proteins %>%
    filter(plasmid_count >= 1) %>%
    inner_join(gbk.annotation) %>%
    filter(str_detect(.$product,antibiotic.keywords))
## and write out to file.
write.csv(plasmid.ARGs, "../results/all-plasmid-ARGs.csv",
          row.names=FALSE,quote=F)

######## Lingchong asked for this control analysis.
## by default, don't count plasmid proteins as duplicates.
## The results are robust to this assumption; nothing changes.
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
## Use CARD and MobileOG-db to check the precision and recall of keywords on
## duplicated ARGs and duplicated MGE-associated genes.

duplicate.proteins.in.CARD <- all.proteins.in.CARD %>%
    filter(count > 1)
duplicate.proteins.in.mobileOGdb <- all.proteins.in.mobileOGdb %>%
    filter(count > 1)

## Measure precision and recall of keywords, treating CARD hits as ground truth.
duplicate.ARGs.by.keyword <- duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords))

## True positives (CARD)
duplicate.proteins.in.CARD.matched.by.keywords <- duplicate.proteins.in.CARD %>%
    filter(str_detect(.$product,antibiotic.keywords))

## False negatives (CARD)
duplicate.proteins.in.CARD.not.matched.by.keywords <- duplicate.proteins.in.CARD %>%
    filter(!str_detect(.$product,antibiotic.keywords))

## false positives (CARD)
duplicate.ARGs.not.in.CARD <- duplicate.ARGs.by.keyword %>%
    filter(!(SeqID %in% duplicate.proteins.in.CARD$SeqID)) %>%
    arrange(product) %>% select(-sequence) %>% as_tibble()


## duplicated ARG precision = 0.931
length(duplicate.proteins.in.CARD.matched.by.keywords$SeqID)/(length(duplicate.ARGs.not.in.CARD$SeqID) + length(duplicate.proteins.in.CARD.matched.by.keywords$SeqID))

## duplicated ARG recall = 0.972
length(duplicate.proteins.in.CARD.matched.by.keywords$SeqID)/(length(duplicate.proteins.in.CARD.not.matched.by.keywords$SeqID) + length(duplicate.proteins.in.CARD.matched.by.keywords$SeqID))

##############
## Measure precision and recall of keywords, treating mobileOG-db hits as ground truth.
duplicate.MGE.genes.by.keyword <- duplicate.proteins %>%
    filter(str_detect(.$product,MGE.keywords))

## True positives (mobileOG-db)
duplicate.proteins.in.mobileOGdb.matched.by.keywords <- duplicate.proteins.in.mobileOGdb %>%
    filter(str_detect(.$product,MGE.keywords))

## False negatives (mobileOG-db)
duplicate.proteins.in.mobileOGdb.not.matched.by.keywords <- duplicate.proteins.in.mobileOGdb %>%
    filter(!str_detect(.$product,MGE.keywords))

## False positives (mobileOG-db)
duplicate.MGE.genes.not.in.mobileOGdb <- duplicate.MGE.genes.by.keyword %>%
    filter(!(SeqID %in% duplicate.proteins.in.mobileOGdb$SeqID)) %>%
    as_tibble()

## duplicated MGE genes precision = 0.786
length(duplicate.proteins.in.mobileOGdb.matched.by.keywords$SeqID)/(length(duplicate.MGE.genes.not.in.mobileOGdb$SeqID) + length(duplicate.proteins.in.mobileOGdb.matched.by.keywords$SeqID))

## duplicated MGE genes recall = 0.872
length(duplicate.proteins.in.mobileOGdb.matched.by.keywords$SeqID)/(length(duplicate.proteins.in.mobileOGdb.not.matched.by.keywords$SeqID) + length(duplicate.proteins.in.mobileOGdb.matched.by.keywords$SeqID))

##########################################################################
## IMPORTANT: global data structures used THROUGHOUT the entire data analysis.

if (USE.CARD.AND.MOBILE.OG.DB) {
    duplicate.proteins <- duplicate.proteins %>%
        ## add columns to duplicate proteins,
        ## based on whether they are in CARD or mobileOG-db or not.
        mutate(in.CARD = ifelse(SeqID %in%all.proteins.in.CARD$SeqID, 1, 0)) %>%
        mutate(in.mobileOGdb = ifelse(SeqID %in%all.proteins.in.mobileOGdb$SeqID, 1, 0)) %>%
        mutate(Category = ifelse(in.CARD, "ARG", ifelse(in.mobileOGdb, "MGE", "Other function")))

    singleton.proteins <- singleton.proteins %>%
        ## add columns to singleton proteins,
        ## based on whether they are in CARD or mobileOG-db or not.
        mutate(in.CARD = ifelse(SeqID %in%all.proteins.in.CARD$SeqID, 1, 0)) %>%
        mutate(in.mobileOGdb = ifelse(SeqID %in%all.proteins.in.mobileOGdb$SeqID, 1, 0)) %>%
        mutate(Category = ifelse(in.CARD, "ARG", ifelse(in.mobileOGdb, "MGE", "Other function")))
    
    duplicate.ARGs <- duplicate.proteins.in.CARD
    duplicate.MGE.genes <- duplicate.proteins.in.mobileOGdb 

    singleton.ARGs <- all.proteins.in.CARD %>%
        filter(count == 1)
        
    singleton.MGE.genes <- all.proteins.in.mobileOGdb %>%
        filter(count == 1)
} else {

    duplicate.proteins <- duplicate.proteins %>%
        mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other))
    
    singleton.proteins <- singleton.proteins %>%
        mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other))
    
    duplicate.ARGs <- duplicate.ARGs.by.keyword
    duplicate.MGE.genes <- duplicate.MGE.genes.by.keyword

    singleton.ARGs <- singleton.proteins %>%
        filter(str_detect(.$product, antibiotic.keywords))

    singleton.MGE.genes <- singleton.proteins %>%
        filter(str_detect(.$product, MGE.keywords))
}


##########################################################################
## Supplementary Data Files 3 and 4.

## Supplementary Data File 4:  Duplicated ARGs and their annotation.
## Need to run this code first, since S3DataFile construction depends on S4DataFile.
S4DataFile <- duplicate.ARGs %>%
    as_tibble() %>%
    arrange(Annotation_Accession)
write.csv(x=S4DataFile, file="../results/FileS4-Duplicated-ARGs.csv",quote=F,row.names=F)

## Supplementary Data File 3: All genomes, and whether or not they contain duplicated ARGs.
S3DataFile <- gbk.annotation %>%
    mutate(hasDuplicatedARGs = sapply(
               Annotation_Accession,
               function(x)
                   ifelse(x %in% S4DataFile$Annotation_Accession,TRUE,FALSE)))
write.csv(x=S3DataFile,
          file="../results/FileS3-Complete-Genomes-with-Duplicated-ARG-annotation.csv",
          quote=F, row.names=F)

##########################################################################
## Figure 1ABCD. A deterministic ODE model demonstrates that selection can
## drive the evolution of duplicated ARGs on plasmids.

## The panels of this figure are generated in my Pluto notebook:
## duplication-linear-ODE-model.jl.
##########################################################################
## Figure 1E is an analysis of the 9-day Tet50 evolution experiment.
##########################################################################
## Figure 2 is an analysis of the one-day Tet5 evolution experiment.
##########################################################################
## Figure 3 is a schematic of the analysis pipeline.
##########################################################################
## Code and data structures for Figure 4ABC.

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


make.TableS1 <- function(gbk.annotation, duplicate.ARGs) {

    ## count the number of isolates with duplicated ARGs in each category.
    ARG.category.counts <- duplicate.ARGs %>%
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

## Data structure for Figure 4A:
## normal-approximation confidence intervals for the percentage
## of isolates with duplicated ARGs.
TableS1 <- make.TableS1(gbk.annotation, duplicate.ARGs)
## write Supplementary Table S1 to file.
write.csv(x=TableS1, file="../results/TableS1.csv")

## Yi asked me to make these comparisons.
TableS1.chromosome.dups <- make.TableS1(gbk.annotation,
                                             filter(duplicate.ARGs, chromosome_count > 1))

TableS1.plasmid.dups <- make.TableS1(gbk.annotation,
                                             filter(duplicate.ARGs, plasmid_count > 1))

######################
## Table S2. Control: does the distribution of ARG singletons
## (i.e. genes that have NOT duplicated) follow the distribution
## of sampled isolates?

## No categories are enriched with ARG singletons,
## as most isolates have a gene that matches an antibiotic keyword.
## Animal-host isolates are depleted (perhaps due to aphid bacteria isolates?)

make.TableS2 <- function(gbk.annotation, singleton.ARGs) {

## count the number of isolates with singleton AR genes in each category.
    ARG.category.counts <- singleton.ARGs %>%
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

## This data frame will be used for Figure 4B.
TableS2 <- make.TableS2(gbk.annotation, singleton.ARGs)
## write TableS2 to file.
write.csv(x=TableS2, file="../results/TableS2.csv")
#########################################################################
## make Supplementary Figure S3.
TableS2.chromosome.only <- make.TableS2(gbk.annotation,
                                             filter(singleton.ARGs, plasmid_count == 0))

TableS2.plasmid.only <- make.TableS2(gbk.annotation,
                                             filter(singleton.ARGs, chromosome_count == 0))
gc() ## free memory after dealing with singleton data.

S3FigA <- make.confint.figure.panel(TableS1.chromosome.dups, order.by.total.isolates, "D-ARGs (only chromosome)")
S3FigB <- make.confint.figure.panel(TableS1.plasmid.dups, order.by.total.isolates, "D-ARGs (only plasmid)",
                                        no.category.label = TRUE)
S3FigC <- make.confint.figure.panel(TableS2.chromosome.only, order.by.total.isolates, "S-ARGs (only chromosome)")
S3FigD <- make.confint.figure.panel(TableS2.plasmid.only, order.by.total.isolates, "S-ARGs (only plasmid)",
                                    no.category.label = TRUE)

S3Fig <- plot_grid(S3FigA, S3FigB, S3FigC, S3FigD, labels = c("A", "B","C",'D'), nrow = 2, rel_widths = c(1.3, 1, 1.3, 1))

if (USE.CARD.AND.MOBILE.OG.DB) { ## then this is Supplementary Figure S15.
    ggsave("../results/S15Fig.pdf", S3Fig, width=6.5, height=3.75)
} else { ## This is supplementary Figure S3.
    ggsave("../results/S3Fig.pdf", S3Fig, width=6.5, height=3.75)
}

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

## Data structure for Figure 4C.
TableS3 <- make.TableS3(gbk.annotation, duplicate.proteins)
## write TableS3 to file.
write.csv(x=TableS3, file="../results/TableS3.csv")

######################################################################
## Supplementary Figure S4: Control for Taxonomy (simpler control than for phylogeny)

## let's look at the taxonomic distribution of strains with duplicated ARGs.
duplicated.ARG.seq.genera.summary <- duplicate.ARGs %>%
    as_tibble() %>%
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

duplicated.genera.isolate.summary <- duplicate.ARGs %>%
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
S4FigA <- genera.isolate.comparison.df %>%
    ggplot(aes(x=sqrt(genome.count),
               y = sqrt(duplicated.ARG.genome.count),
               label = Genus,
               color = in.top.ARG.genera)) +
    theme_classic() + geom_jitter() + geom_text_repel(fontface = "italic") +
    scale_color_manual(values=c("black", "red")) +
    guides(color="none") +
    xlab("sqrt(Number of isolates)") +
    ylab("sqrt(Number of isolates with D-ARGs)")

## 5,638 isolates are in the top ARG genera.
top.ARG.genera.isolates <- gbk.annotation %>%
    filter(Genus %in% top.ARG.genera)

## 13,300 isolates are in the remaining genera.
filtered.gbk.annotation <- gbk.annotation %>%
    filter(!(Genus %in% top.ARG.genera))

filtered.duplicate.ARGs <- duplicate.ARGs %>%
    filter(!(Genus %in% top.ARG.genera))

filtered.TableS1 <- make.TableS1(filtered.gbk.annotation, filtered.duplicate.ARGs)

S4FigB <- make.confint.figure.panel(filtered.TableS1, order.by.total.isolates, "D-ARGs after\nfiltering top genera")

## let's check the results, just for the top ARG genera.
top.ARG.genera.duplicate.ARGs <- duplicate.ARGs %>%
    filter(Genus %in% top.ARG.genera)

top.ARG.genera.TableS1 <- make.TableS1(top.ARG.genera.isolates,
                                       top.ARG.genera.duplicate.ARGs)

S4FigC <- make.confint.figure.panel(
    top.ARG.genera.TableS1, order.by.total.isolates,
    "D-ARGs in\ntop genera only",
    no.category.label = TRUE)

## let's downsample the data, using Assembly-dereplicator.
dereplicated.genomes <- read.table("../results/dereplicated-genomes.txt",header=FALSE) %>%
    rename(fasta_file = V1) %>%
    mutate(Annotation_Accession=str_replace(fasta_file, "_genomic.fna.gz","")) %>%
    as_tibble()

dereplicated.gbk.annotation <-  gbk.annotation %>%
    filter(Annotation_Accession %in% dereplicated.genomes$Annotation_Accession)

dereplicated.duplicate.ARGs <-  duplicate.ARGs %>%
    filter(Annotation_Accession %in% dereplicated.genomes$Annotation_Accession)

dereplicated.TableS1 <- make.TableS1(
    dereplicated.gbk.annotation, dereplicated.duplicate.ARGs)

S4FigD <- make.confint.figure.panel(
    dereplicated.TableS1, order.by.total.isolates,
    "D-ARGs after downsampling\nby Mash distance > 0.005",
        no.category.label = FALSE)

## Let's try an alternative strategy: downsample the data such that only one
## sample for each organism is allowed.

single.organism.gbk.annotation <- gbk.annotation %>%
    group_by(Organism) %>%
    filter(row_number() == 1) %>% ## take the first one in the group.
    ungroup()

single.organism.duplicate.ARGs <- duplicate.ARGs %>%
    filter(Annotation_Accession %in% single.organism.gbk.annotation$Annotation_Accession)

single.organism.TableS1 <- make.TableS1(
    single.organism.gbk.annotation, single.organism.duplicate.ARGs)

S4FigE <- make.confint.figure.panel(
    single.organism.TableS1, order.by.total.isolates,
    "D-ARGs after\ndownsampling species",
        no.category.label = TRUE)

S4FigBCDE <- plot_grid(S4FigB, S4FigC, S4FigD, S4FigE,
                       labels = c("B","C",'D','E'), nrow = 2, rel_widths = c(1.5, 1, 1, 1, 1))

S4Fig <- plot_grid(S4FigA, S4FigBCDE, labels = c("A", ""), nrow = 2, rel_heights = c(2,2))
ggsave("../results/S4Fig.pdf", S4Fig, height = 8, width=7.5)

##################################################################
## Figures S2 and S5. Proportion of isolates with duplicated or single-copy ARGs
## for 12 different antibiotic classes.
########################################
## Figure S2: Duplicated ARGs.

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


S2FigA <- make.confint.figure.panel(chloramphenicol.table,
                                    order.by.total.isolates,
                                    "chloramphenicol\nresistance")
S2FigB <- make.confint.figure.panel(tetracycline.table,
                                    order.by.total.isolates,
                                    "tetracycline\nresistance",
                                    no.category.label = TRUE)
S2FigC <- make.confint.figure.panel(MLS.table,
                                    order.by.total.isolates,
                                    "MLS\nresistance",
                                    no.category.label = TRUE)
S2FigD <- make.confint.figure.panel(multidrug.table,
                                    order.by.total.isolates,
                                    "multidrug\nresistance")
S2FigE <- make.confint.figure.panel(beta.lactam.table,
                                    order.by.total.isolates,
                                    "beta-lactam\nresistance",
                                    no.category.label = TRUE)
S2FigF <- make.confint.figure.panel(glycopeptide.table,
                                    order.by.total.isolates,
                                    "glycopeptide\nresistance",
                                    no.category.label = TRUE)
S2FigG <- make.confint.figure.panel(polypeptide.table,
                                    order.by.total.isolates,
                                    "polypeptide\nresistance")
S2FigH <- make.confint.figure.panel(diaminopyrimidine.table,
                                    order.by.total.isolates,
                                    "diaminopyrimidine\nresistance",
                                    no.category.label = TRUE)
S2FigI <- make.confint.figure.panel(sulfonamide.table,
                                    order.by.total.isolates,
                                    "sulfonamide\nresistance",
                                    no.category.label = TRUE)
S2FigJ <- make.confint.figure.panel(quinolone.table,
                                    order.by.total.isolates,
                                    "quinolone\nresistance")
S2FigK <- make.confint.figure.panel(aminoglycoside.table,
                                    order.by.total.isolates,
                                    "aminoglycoside\nresistance",
                                    no.category.label = TRUE)
S2FigL <- make.confint.figure.panel(macrolide.table,
                                    order.by.total.isolates,
                                    "macrolide\nresistance",
                                    no.category.label = TRUE)

S2Fig <- plot_grid(NULL, ## The nesting is to add a title.
                   plot_grid(
                       S2FigA, S2FigB, S2FigC,
                       S2FigD, S2FigE, S2FigF,
                       S2FigG, S2FigH, S2FigI,
                       S2FigJ, S2FigK, S2FigL,
                   rel_widths = c(1.5, 1, 1,
                                  1.5, 1, 1,
                                  1.5, 1, 1),
                   nrow = 4),
                   labels = c("D-ARGs", ""),
                   ncol = 1,
                   rel_heights = c(0.025, 1))

ggsave("../results/S2Fig.pdf", S2Fig, height = 9, width = 13)
########################################
## Figure S5: Single-copy ARGs.

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


S5FigA <- make.confint.figure.panel(chloramphenicol.control.table,
                                      order.by.total.isolates,
                                      "chloramphenicol\nresistance")
S5FigB <- make.confint.figure.panel(tetracycline.control.table,
                                      order.by.total.isolates,
                                      "tetracycline\nresistance",
                                      no.category.label = TRUE)
S5FigC <- make.confint.figure.panel(MLS.control.table,
                                      order.by.total.isolates,
                                      "MLS\nresistance",
                                      no.category.label = TRUE)
S5FigD <- make.confint.figure.panel(multidrug.control.table,
                                      order.by.total.isolates,
                                      "multidrug\nresistance")
S5FigE <- make.confint.figure.panel(beta.lactam.control.table,
                                      order.by.total.isolates,
                                      "beta-lactam\nresistance",
                                      no.category.label = TRUE)
S5FigF <- make.confint.figure.panel(glycopeptide.control.table,
                                      order.by.total.isolates,
                                      "glycopeptide\nresistance",
                                      no.category.label = TRUE)
S5FigG <- make.confint.figure.panel(polypeptide.control.table,
                                      order.by.total.isolates,
                                      "polypeptide\nresistance")
S5FigH <- make.confint.figure.panel(diaminopyrimidine.control.table,
                                      order.by.total.isolates,
                                      "diaminopyrimidine\nresistance",
                                      no.category.label = TRUE)
S5FigI <- make.confint.figure.panel(sulfonamide.control.table,
                                      order.by.total.isolates,
                                      "sulfonamide\nresistance",
                                      no.category.label = TRUE)
S5FigJ <- make.confint.figure.panel(quinolone.control.table,
                                    order.by.total.isolates,
                                    "quinolone\nresistance")
S5FigK <- make.confint.figure.panel(aminoglycoside.control.table,
                                      order.by.total.isolates,
                                      "aminoglycoside\nresistance",
                                      no.category.label = TRUE)
S5FigL <- make.confint.figure.panel(macrolide.control.table,
                                      order.by.total.isolates,
                                      "macrolide\nresistance",
                                      no.category.label = TRUE)

S5Fig <- plot_grid(NULL, ## The nesting is to add a title.
                   plot_grid(
                       S5FigA, S5FigB, S5FigC,
                       S5FigD, S5FigE, S5FigF,
                       S5FigG, S5FigH, S5FigI,
                       S5FigJ, S5FigK, S5FigL,
                   rel_widths = c(1.5, 1, 1,
                                  1.5, 1, 1,
                                  1.5, 1, 1),
                   nrow = 4),
                   labels = c("S-ARGs", ""),
                   ncol = 1,
                   rel_heights = c(0.025, 1))
ggsave("../results/S5Fig.pdf", S5Fig, height = 9, width = 8)

#########################
## S11 Figure.
## Analysis of duplicate pairs found just on chromosome, just on plasmid, or
## on both chromosomes and plasmids.

## let's look at cases of identical sequences on chromosomes and plasmids.
both.chr.and.plasmid.cases <- duplicate.proteins %>%
    filter(chromosome_count >= 1 & plasmid_count >= 1) %>%
    arrange(desc(count)) %>%
    tibble()

just.chromosome.cases <- duplicate.proteins %>%
    filter(chromosome_count >= 1 & plasmid_count == 0) %>%
    arrange(desc(count)) %>%
    tibble()

just.plasmid.cases <- duplicate.proteins %>%
    filter(chromosome_count == 0 & plasmid_count >= 1) %>%
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

S11FigA <- ggplot(both.chr.and.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("Both chromosome and plasmid") +
    theme(legend.position="bottom") +
    xlab("Frequency") +
    ylab("")

S11Fig.legend <- get_legend(S11FigA)
S11FigA <- S11FigA + guides(fill = "none")

S11FigB <- ggplot(just.chromosome.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("Chromosome only") +
    guides(fill = "none") +
    xlab("Frequency") +
    ylab("")

S11FigC <- ggplot(just.plasmid.summary,
                  aes(x = Count,
                      y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    theme_classic() +
    ggtitle("Plasmid only") +
    guides(fill = "none") +
    xlab("Frequency") +
    ylab("")

S11Fig <- plot_grid(NULL, S11FigA, S11FigB, S11FigC, S11Fig.legend, ncol = 1,
                   labels = c("Genomic distribution of D-genes", "A","B","C"),
                   rel_heights=c(0.35,1,1,1,0.25))
ggsave("../results/S11Fig.pdf", S11Fig, width=5, height=8)


## Run some statistics to explicitly test whether
## duplicated genes encoded solely on plasmids are more likely to encode antibiotic resistance
## and functions other than those associated with mobile genetic elements,
## in comparison to both duplicated genes encoded solely on the chromosome,
## and duplicated genes encoded on plasmids and the chromosome, as is seen in S11Fig.

category.summed.both.chr.and.plasmid <- both.chr.and.plasmid.summary %>%
    group_by(Category) %>%
    summarize(summed_count = sum(Count))

category.summed.just.chromosome <- just.chromosome.summary %>%
    group_by(Category) %>%
    summarize(summed_count = sum(Count))

category.summed.just.plasmid <- just.plasmid.summary %>%
    group_by(Category) %>%
    summarize(summed_count = sum(Count))

both.chr.and.plasmid.test.vec <- c(
    filter(category.summed.both.chr.and.plasmid,Category=='ARG')$summed_count
    + filter(category.summed.both.chr.and.plasmid,Category=='Other function')$summed_count,
    filter(category.summed.both.chr.and.plasmid,Category=='MGE')$summed_count)

just.plasmid.test.vec <- c(
    filter(category.summed.just.plasmid,Category=='ARG')$summed_count
    + filter(category.summed.just.plasmid,Category=='Other function')$summed_count,
    filter(category.summed.just.plasmid,Category=='MGE')$summed_count)

just.chromosome.test.vec <- c(
    filter(category.summed.just.chromosome,Category=='ARG')$summed_count
    + filter(category.summed.just.chromosome,Category=='Other function')$summed_count,
    filter(category.summed.just.chromosome,Category=='MGE')$summed_count)

## test 1: compare proportions between just.plasmid and just.chromosome.
binom.test(just.plasmid.test.vec,p=just.chromosome.test.vec[1]/sum(just.chromosome.test.vec))
binom.test(just.plasmid.test.vec,p=just.chromosome.test.vec[1]/sum(just.chromosome.test.vec))$p.value

## test 2: compare proportions between just.plasmid and both plasmid and chromosome.
binom.test(just.plasmid.test.vec,p=both.chr.and.plasmid.test.vec[1]/sum(both.chr.and.plasmid.test.vec))
binom.test(just.plasmid.test.vec,p=both.chr.and.plasmid.test.vec[1]/sum(both.chr.and.plasmid.test.vec))$p.value


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
    duplicate.chromosome.ARGs.count <- duplicate.ARGs %>%
        group_by(Annotation) %>%
        summarize(chromosomal_duplicate_ARGs = sum(chromosome_count))
    
    ## Column 4
    duplicate.plasmid.ARGs.count <- duplicate.ARGs %>%
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
    singleton.chromosome.ARGs.count <- singleton.ARGs %>%
        group_by(Annotation) %>%
        summarize(chromosomal_singleton_ARGs = sum(chromosome_count))
    gc() ## free memory when dealing with singleton.proteins.
    
    ## Column 4
    singleton.plasmid.ARGs.count <- singleton.ARGs %>%
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
## Use the data in Tables S4 and S5 to make the data structures for Figure 4DEFGHI.
## The point of this figure is to show that the distribution of
## duplicated ARGs is not predicted by the distribution of single-copy ARGs
## in the ecological categories.

make.Fig4DEFGHI.df <- function(TableS1, TableS4, TableS5) {
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


Fig4DEFGHI.df <- make.Fig4DEFGHI.df(TableS1, TableS4, TableS5)
## Figure 4DEFGHI.
## Plot point estimates for the fraction of chromosomal genes that are
## the fraction of genes that are duplicated ARGs (panel D),
## the fraction of chromosomal genes that are duplicated ARGs (panel E),
## the fraction of plasmid genes that are duplicated ARGs (panel F),
## the fraction of genes that are single-copy ARGs (panel G).
## the fraction of chromosomal genes that are single-copy ARGs (panel H),
## the fraction of plasmid genes that are single-copy ARGs (panel I),

## Fig4D
Fig4D.df <- Fig4DEFGHI.df %>%
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


## Fig4E: the fraction of chromosomal genes that are duplicated ARGs.
Fig4E.df <- Fig4DEFGHI.df %>%
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

Fig4F.df <- Fig4DEFGHI.df %>%
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

Fig4G.df <- Fig4DEFGHI.df %>%
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

Fig4H.df <- Fig4DEFGHI.df %>%
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


Fig4I.df <- Fig4DEFGHI.df %>%
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


make.Fig4DEFGHI.panel <- function(Table, order.by.total.isolates, title,
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


## Finally -- make Figure 4ABC.
## Throughout, add special scales for Figure 4 but not for Supplementary Figure S14.
Fig4A <- make.confint.figure.panel(TableS1, order.by.total.isolates, "D-ARGs")
if (! USE.CARD.AND.MOBILE.OG.DB) {
    Fig4A <- Fig4A +
        scale_x_continuous(breaks = c(0, 0.15), limits = c(0,0.16))
}

Fig4B <- make.confint.figure.panel(TableS2, order.by.total.isolates,
                                   "S-ARGs", no.category.label=TRUE)
if (! USE.CARD.AND.MOBILE.OG.DB) {
    Fig4B <- Fig4B +
        scale_x_continuous(breaks = c(0.85, 1.0), limits = c(0.85,1))
}

Fig4C <- make.confint.figure.panel(TableS3, order.by.total.isolates,
                                   "All D-genes", no.category.label=TRUE)
if (! USE.CARD.AND.MOBILE.OG.DB) {
    Fig4C <- Fig4C +
        scale_x_continuous(breaks = c(0.75, 0.95), limits = c(0.75, 1.0))
}

Fig4ABC.title <- title_theme <- ggdraw() +
    draw_label("Isolate-level analysis",fontface="bold")

Fig4ABC <- plot_grid(Fig4A, Fig4B, Fig4C, labels=c('A','B','C'),
                     rel_widths = c(1.5,1,1), nrow=1)

Fig4ABC.with.title <- plot_grid(Fig4ABC.title, Fig4ABC, ncol = 1, rel_heights = c(0.1, 1))

## the rest of Figure 4 -- Fig4DEFGHI -- requires Tables S4 and S5.
## I manually set axis labels so that they don't run into each other.
Fig4D <- make.Fig4DEFGHI.panel(Fig4D.df, order.by.total.isolates,
                         "\nD-ARGs",
                         "Proportion of\nall genes") +
    ## This scale is for both Figure 4D and S14DFig.
    scale_x_continuous(label=fancy_scientific, breaks = c(0, 2e-4), limits = c(0,2.5e-4))


Fig4E <- make.Fig4DEFGHI.panel(Fig4E.df, order.by.total.isolates,
                         "Chromosome:\nD-ARGs",
                         "Proportion of\nchromosomal genes",
                         no.category.label = TRUE)
if (! USE.CARD.AND.MOBILE.OG.DB) {
    Fig4E <- Fig4E +
        scale_x_continuous(label=fancy_scientific, breaks = c(0, 8e-5), limits = c(0,1.2e-4))
}

Fig4F <- make.Fig4DEFGHI.panel(Fig4F.df, order.by.total.isolates,
                         "Plasmids:\nD-ARGs",
                         "Proportion of\nplasmid genes",
                         no.category.label = TRUE)
if (! USE.CARD.AND.MOBILE.OG.DB) {
    Fig4F <- Fig4F +
         scale_x_continuous(label=fancy_scientific, breaks = c(0, 5e-3), limits = c(0,6e-3))
}

Fig4G <- make.Fig4DEFGHI.panel(Fig4G.df, order.by.total.isolates,
                         "\nS-ARGs",
                         "Proportion of\nall genes")         
if (! USE.CARD.AND.MOBILE.OG.DB) {
    Fig4G <- Fig4G +
        scale_x_continuous(label=fancy_scientific, breaks = c(2e-3, 4e-3), limits = c(1.5e-3, 4.5e-3))
}

Fig4H <- make.Fig4DEFGHI.panel(Fig4H.df, order.by.total.isolates,
                         "Chromosome:\nS-ARGs",
                         "Proportion of\nchromosomal genes",
                         no.category.label = TRUE) +
    ## This scale is for S14DFig.
    scale_x_continuous(label=fancy_scientific, breaks = c(2e-3, 4e-3, 6e-3), limits = c(0, 6.25e-3))
if (! USE.CARD.AND.MOBILE.OG.DB) {
    Fig4H <- Fig4H +
        scale_x_continuous(label=fancy_scientific, breaks = c(2e-3, 3e-3), limits = c(1.5e-3, 3.5e-3))
}

Fig4I <- make.Fig4DEFGHI.panel(Fig4I.df, order.by.total.isolates,
                         "Plasmids:\nS-ARGs",
                         "Proportion of\nplasmid genes",
                         no.category.label = TRUE) +
    ## This scale is for S14IFig.
    scale_x_continuous(label=fancy_scientific, breaks = c(3e-3, 2e-2), limits = c(0,2.5e-2))
if (! USE.CARD.AND.MOBILE.OG.DB) {
    Fig4I <- Fig4I +
        scale_x_continuous(label=fancy_scientific, breaks = c(3e-3, 2e-2))
}

Fig4DEFGHI.title <- title_theme <- ggdraw() +
    draw_label("Gene-level analysis",fontface="bold")

Fig4DEFGHI <- plot_grid(Fig4D, Fig4E, Fig4F, Fig4G, Fig4H, Fig4I,
                  labels = c('D','E','F','G','H','I'), nrow=2,
                  rel_widths = c(1.5, 1, 1, 1.5, 1, 1))

Fig4DEFGHI.with.title <- plot_grid(Fig4DEFGHI.title, Fig4DEFGHI,
                                   ncol = 1, rel_heights = c(0.06, 1))

## Now, make the complete Figure 4!
Fig4 <- plot_grid(Fig4ABC.with.title, Fig4DEFGHI.with.title,
                  ncol = 1, rel_heights = c(0.3,0.7))

if (! USE.CARD.AND.MOBILE.OG.DB) { ## Make main figure 4.
    ggsave("../results/Fig4.pdf", Fig4, height=6.25, width=6.25)
} else { ## Make supplementary Figure S14.
    ggsave("../results/S14Fig.pdf", Fig4, height=6.25, width=8)
}

##########################################################################
## Figure 5: Visualization of ARGs on plasmids and chromosomes, and evidence for selection.
## set up data structures for Figure 5AB.
Fig5A.data <- duplicate.proteins %>%
    group_by(Annotation, Category) %>%
    summarize(Plasmid = sum(plasmid_count), Chromosome = sum(chromosome_count)) %>%
    pivot_longer(cols = c("Plasmid", "Chromosome"),
                 names_to = "Episome",
                 values_to = "Count") %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

Fig5B.data <- singleton.proteins %>%
    group_by(Annotation, Category) %>%
    summarize(Plasmid = sum(plasmid_count), Chromosome = sum(chromosome_count)) %>%
    pivot_longer(cols = c("Plasmid", "Chromosome"),
                 names_to = "Episome",
                 values_to = "Count") %>%
    mutate(Annotation = factor(
               Annotation,
               levels = rev(order.by.total.isolates)))

Fig5A <- ggplot(Fig5A.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    facet_wrap(.~Episome) +
    theme_classic() +
    xlab("Proportion of D-genes") +
    scale_x_continuous(breaks = c(0,1)) +
    theme(legend.position="bottom") +
    theme(strip.background = element_blank()) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

Fig5legend <- get_legend(Fig5A)
Fig5A <- Fig5A + guides(fill = "none")

Fig5B <- ggplot(Fig5B.data, aes(x = Count, y = Annotation, fill = Category)) +
    geom_bar(stat="identity", position = "fill", width = 0.95) +
    facet_wrap(.~Episome) +
    theme_classic() +
    xlab("Proportion of S-genes") +
    scale_x_continuous(breaks = c(0,1)) +
    guides(fill = "none") +
    theme(strip.background = element_blank()) +
    theme(axis.text.y=element_blank()) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

## let's get the actual numbers of duplicated MGE genes in the plasmid and chromosome,
## to report in the main text.
category.summed.D.genes <- Fig5A.data %>%
    group_by(Category, Episome) %>%
    summarize(summed_count = sum(Count))

category.summed.S.genes <- Fig5B.data %>%
    group_by(Category, Episome) %>%
    summarize(summed_count = sum(Count))


## Figure 5C: 
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
                                   order.by.total.isolates, category.string) {

    duplicated.function.per.category <- duplicate.proteins %>%
        filter(Category == category.string) %>%
        group_by(Annotation) %>%
        summarize(function.duplicates = sum(count))
    
    singleton.function.per.category <- singleton.proteins %>%
        filter(Category == category.string) %>%
        group_by(Annotation) %>%
        summarize(function.singletons = sum(count))

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
                   levels = rev(order.by.total.isolates))) %>%
        mutate(Category = category.string)
    
    return(selection.test.df)
}


ARG.selection.test.df <- make.selection.test.df(
    duplicate.proteins, singleton.proteins,
    order.by.total.isolates, "ARG") %>%
    select(Annotation, Category, dup.singleton.ratio)

MGE.selection.test.df <- make.selection.test.df(
    duplicate.proteins, singleton.proteins,
    order.by.total.isolates, "MGE") %>%
    select(Annotation, Category, dup.singleton.ratio)

other.selection.test.df <- make.selection.test.df(
    duplicate.proteins, singleton.proteins,
    order.by.total.isolates, "Other function") %>%
    select(Annotation, Category, dup.singleton.ratio)

big.selection.test.df <- ARG.selection.test.df %>%
    full_join(MGE.selection.test.df) %>%
    full_join(other.selection.test.df)


Fig5C <- big.selection.test.df %>%
    ggplot(aes(y = Annotation, x = log(dup.singleton.ratio), color = Category)) +
    geom_point() + theme_classic() +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    xlab("log(% of D-genes / % of S-genes)") +
    ##xlab(TeX("\\frac{% of D-genes}{% of S-genes}")) +
    guides(color = "none") +
    theme(axis.text.y=element_blank()) +
    ylab("") ## remove the redundant "Annotation" label on the y-axis.

################################################################################
## Write Fig5 ABC to file. Fig 5D is a diagram of a workflow, and Figure 5E is produced at the end
## of this script.
Fig5.top.panels <- plot_grid(Fig5A, Fig5B, Fig5C, labels = c('A','B','C'), nrow = 1, rel_widths=c(1.4,1,1),
              align = 'h', axis = 'tb')
Fig5ABC <- plot_grid(Fig5.top.panels, Fig5legend, ncol = 1, rel_heights = c(1,0.2))

if (USE.CARD.AND.MOBILE.OG.DB) { ## Then save as a Supplementary Figure S16.
    ggsave("../results/S16Fig.pdf", Fig5ABC, width=8.75, height=4)
} else { ## save as a main Figure.
    ggsave("../results/Fig5ABC.pdf", Fig5ABC, width=8.75, height=4)
}


################################################################################
## compare linkage between D-ARGs and MGE-genes and S-ARGs and MGE-genes.
## These calculations go into the manuscript, but no specific figure.

ARG.MGE.adjacency.data <- read.csv("../results/ARG-MGE-adjacency-counts.csv") %>%
    pivot_longer(cols=everything(), names_to="original_column", values_to = "Count") %>%
    mutate(ARG.class = ifelse(str_detect(original_column,"dARG"),"D-ARGs","S-ARGs")) %>%
    mutate(next.to.MGE = ifelse(str_detect(original_column,"not"),FALSE,TRUE))

ARG.MGE.adjacency.total <- ARG.MGE.adjacency.data %>% group_by(ARG.class) %>%
    summarize(ARG.total=sum(Count))

ARGs.next.to.MGEs <- filter(ARG.MGE.adjacency.data,
                             next.to.MGE==TRUE) %>%
    rename(next.to.MGE.count=Count) %>%
    select(ARG.class, next.to.MGE.count)

ARG.MGE.adjacency.df <- full_join(ARG.MGE.adjacency.total, ARGs.next.to.MGEs) %>%
    mutate(percent.next.to.MGE = next.to.MGE.count/ARG.total)

## formally calculate statistical significance with a binomial test.
ARG.MGE.adjacency.statistic <- binom.test(
    x=filter(ARG.MGE.adjacency.df,ARG.class=="D-ARGs")$next.to.MGE.count, ## x = 4651
    n=filter(ARG.MGE.adjacency.df,ARG.class=="D-ARGs")$ARG.total, ## n = 9836
    p = filter(ARG.MGE.adjacency.df,ARG.class=="S-ARGs")$percent.next.to.MGE)


################################################################################
## Analysis of clinical antibiotic resistant isolates.

make.clinical.genomes.D.ARG.Figure <- function(clinical.duplicate.proteins, clinical.singleton.proteins, panelA.title) {
    
    FigA.data <- clinical.duplicate.proteins %>%
        mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
        group_by(Annotation_Accession, Category) %>%
        summarize(Count = sum(count))

    FigB.data <- clinical.singleton.proteins %>%
        mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other)) %>%
        group_by(Annotation_Accession, Category) %>%
        summarize(Count = sum(count))

    FigA1 <- ggplot(FigA.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
        geom_bar(stat="identity") +
        theme_classic() +
        theme(legend.position="top") +
        ylab("") ## remove the redundant "Annotation" label on the y-axis.

    Figlegend <- get_legend(FigA1)
    FigA1 <- FigA1 + guides(fill = "none")

    FigA2 <- ggplot(FigA.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
        geom_bar(stat="identity", position = "fill") +
        theme_classic() +
        ## remove genome name labels.
        theme(axis.text.y=element_blank()) +
        guides(fill = "none") +
        xlab("Frequency") +
        ylab("") ## remove the redundant "Annotation" label on the y-axis.

    FigA.title <- ggdraw() + draw_label(panelA.title, fontface='bold')

    FigA <- plot_grid(FigA.title, Figlegend,
                   plot_grid(FigA1, FigA2, labels = c("A",""),
                             nrow=1, rel_widths=c(2,1)),
                   nrow=3, rel_heights=c(0.1,0.1,2))

    FigB1 <- ggplot(FigB.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
        geom_bar(stat="identity") +
        theme_classic() +
        guides(fill = "none") +
        ## remove genome name labels.
        theme(axis.text.y=element_blank()) +
        ylab("") ## remove the redundant "Annotation" label on the y-axis.

    FigB2 <- ggplot(FigB.data, aes(x = Count, y = Annotation_Accession, fill = Category)) +
        geom_bar(stat="identity", position = "fill") +
        theme_classic() +
        ## remove genome name labels.
        theme(axis.text.y=element_blank()) +
        guides(fill = "none") +
        xlab("Frequency") +
        ylab("") ## remove the redundant "Annotation" label on the y-axis.

    FigB.title <- ggdraw() + draw_label("Corresponding distribution\nof single-copy genes", fontface='bold')

    FigB <- plot_grid(FigB.title,
                      plot_grid(FigB1, FigB2, labels = c("B",""),
                                nrow=1, rel_widths=c(1,1)),
                      nrow=2, rel_heights=c(0.2,2))
    
    FullFig <- plot_grid(FigA, FigB, ncol = 2, rel_widths = c(1.5,1))
    return(FullFig)
}

################################################################################
## Supplementary Figure S6.
## Let's analyze duplicated genes in the 12 GN0XXXX genomes that were sequenced with
## long-read technology (PacBio) by Vance Fowler's lab.
## Jon Bethke characterized the resistances of these strains, and additionally
## sequenced another 7 ESBL genomes with PacBio technology.

## 58,092 types of protein sequences in the 12 ESBL genomes.
Duke.ESBL.all.proteins <- data.table::fread("../results/Duke-ESBL-all-proteins.csv",
                                  drop="sequence")

## assert that this is an independent dataset.
stopifnot(nrow(filter(Duke.ESBL.all.proteins, Annotation_Accession %in% gbk.annotation$Annotation_Accession)) == 0)

## 57,263 single-copy protein sequences.
Duke.ESBL.singleton.proteins <- Duke.ESBL.all.proteins %>%
    filter(count == 1)

## 829 protein sequences have duplicates in the 12 ESBL genomes (not counting the duplicates)
Duke.ESBL.duplicate.proteins <- read.csv(
    "../results/Duke-ESBL-duplicate-proteins.csv") %>%
    select(-sequence)

Duke.ESBL.duplicate.ARGs <- Duke.ESBL.duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords))
sum(Duke.ESBL.duplicate.ARGs$count)

## 6/12 strains have duplicate ARGs.
length(unique(Duke.ESBL.duplicate.ARGs$Annotation_Accession))
length(unique(Duke.ESBL.singleton.proteins$Annotation_Accession))

## Now make Supplementary Figure S6.
S6Fig <- make.clinical.genomes.D.ARG.Figure(Duke.ESBL.duplicate.proteins, Duke.ESBL.singleton.proteins,
                                            "Distribution of duplicated genes in 12 ESBL-resistant isolates\nfrom Duke Hospital")
ggsave("../results/S6Fig.pdf", S6Fig, height = 9, width = 10)

################################################################################
## Supplementary Figure S7.
## Let's analyze duplicated genes in 46 genomes that were sequenced with
## long-read technology by BARNARDS group in Nature Microbiology (2022).

BARNARDS.all.proteins <- data.table::fread("../results/BARNARDS-all-proteins.csv",
                                  drop="sequence")

## assert that this is an independent dataset.
stopifnot(nrow(filter(BARNARDS.all.proteins, Annotation_Accession %in% gbk.annotation$Annotation_Accession)) == 0)

BARNARDS.singleton.proteins <- BARNARDS.all.proteins %>%
    filter(count == 1)

BARNARDS.duplicate.proteins <- read.csv(
    "../results/BARNARDS-duplicate-proteins.csv") %>%
    select(-sequence)

BARNARDS.duplicate.ARGs <- BARNARDS.duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords))

## 23/46 strains have duplicate ARGs.
length(unique(BARNARDS.duplicate.ARGs$Annotation_Accession))
length(unique(BARNARDS.singleton.proteins$Annotation_Accession))

## Now make Supplementary Figure S7.
S7Fig <- make.clinical.genomes.D.ARG.Figure(BARNARDS.duplicate.proteins, BARNARDS.singleton.proteins,
                                            "Distribution of duplicated genes in 46 genomes\nfrom the BARNARDS study")
ggsave("../results/S7Fig.pdf", S7Fig, height = 14, width = 8)

################################################################################
## Supplementary Figure S8.
## Let's analyze duplicated genes in 149 genomes that were sequenced with
## long-read technology by  Dantas group in mSystems (2022).

Dantas.all.proteins <- data.table::fread("../results/Mahmud2022-all-proteins.csv",
                                  drop="sequence")
## assert that this is an independent dataset.
stopifnot(nrow(filter(Dantas.all.proteins, Annotation_Accession %in% gbk.annotation$Annotation_Accession)) == 0)

Dantas.singleton.proteins <- Dantas.all.proteins %>%
    filter(count == 1)

Dantas.duplicate.proteins <- read.csv(
    "../results/Mahmud2022-duplicate-proteins.csv") %>%
    select(-sequence)

Dantas.duplicate.ARGs <- Dantas.duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords))

## 36/149 strains have duplicate ARGs.
length(unique(Dantas.duplicate.ARGs$Annotation_Accession))
length(unique(Dantas.singleton.proteins$Annotation_Accession))

## Now make Supplementary Figure S8.
S8Fig <- make.clinical.genomes.D.ARG.Figure(Dantas.duplicate.proteins, Dantas.singleton.proteins,
                                            "Distribution of duplicated genes in 149 genomes\nfrom Barnes-Jewish Hospital (Mahmud et al. 2022)")

ggsave("../results/S8Fig.pdf", S8Fig, height = 18, width = 11)

#######################################################
## Supplementary Figure S9.
## Let's analyze duplicated genes in 114 complete genomes that were sequenced with
## long-read technology by Hawkey et al. in Genome Medicine (2022).
## I downloaded the RefSeq accessions, since most of the Genbank sequences didn't
## have any gene annotations.
Hawkey.all.proteins <- data.table::fread("../results/Hawkey2022-all-proteins.csv",
                                         drop="sequence")

## assert that this is an independent dataset.
stopifnot(nrow(filter(Hawkey.all.proteins, Annotation_Accession %in% gbk.annotation$Annotation_Accession)) == 0)

Hawkey.singleton.proteins <- Hawkey.all.proteins %>%
    filter(count == 1)

Hawkey.duplicate.proteins <- read.csv(
    "../results/Hawkey2022-duplicate-proteins.csv") %>%
    select(-sequence)

Hawkey.duplicate.ARGs <- Hawkey.duplicate.proteins %>%
    filter(str_detect(.$product,antibiotic.keywords))
sum(Hawkey.duplicate.ARGs$count)

## 20/114 strains have duplicate ARGs.
length(unique(Hawkey.duplicate.ARGs$Annotation_Accession))
length(unique(Hawkey.singleton.proteins$Annotation_Accession))

## Now make Supplementary Figure S9.
S9Fig <- make.clinical.genomes.D.ARG.Figure(Hawkey.duplicate.proteins, Hawkey.singleton.proteins,
                                             "Distribution of duplicated genes in 114 genomes\nfrom an Australian ICU (Hawkey et al. 2022)")
ggsave("../results/S9Fig.pdf", S9Fig, height = 14, width = 11)

################################################################################
## Now do the analysis of enrichment across all these clinical datasets.

## 6/12 strains have duplicate ARGs.
## 23/46 strains have duplicate ARGs.
## 36/149 strains have duplicate ARGs.
## 20/114 strains have duplicate ARGs.

## clinical resistance genomes are even more enriched with ARGs than the baseline dataset.
## the null comes from Table S1 row for humans.
binom.test(x=(6+23+36+20),n=(12+46+149+114),p=0.142)

################################################################################
## Supplementary Figure 11. Analysis of copy number in the genomes from Hawkey et al. (2022).
##    The Hawkey et al. 2022 paper specifically focuses on ESBL resistance,
##    so let's focus on beta-lactamases, and can compare to other kinds of resistances in these data.

ARG.copy.number.data <- read.csv("../results/Hawkey2022_ARG_copy_numbers.csv") %>%
    mutate(beta.lactam.resistance = ifelse(str_detect(product,beta.lactam.keywords), TRUE, FALSE))

beta.lactam.ARGs <- filter(ARG.copy.number.data, beta.lactam.resistance==TRUE)
non.beta.lactam.ARGs <- filter(ARG.copy.number.data, beta.lactam.resistance==FALSE)

chromosome.plasmid.copy.number.data <- read.csv("../results/Hawkey2022_chromosome_plasmid_copy_numbers.csv") %>%
    mutate(has.ARG = ifelse(SeqID %in% ARG.copy.number.data$SeqID, TRUE, FALSE)) %>%
    mutate(has.beta.lactamase = ifelse(SeqID %in% beta.lactam.ARGs$SeqID, TRUE, FALSE)) %>%
    ## 0 == no ARG, 1 == has ARG, 2 == has beta-lactamase.
    mutate(ARG.classification = has.ARG + has.beta.lactamase) %>%
    mutate(ARG.classification = as.factor(ARG.classification)) %>%
    mutate(`Plasmid class` = recode(ARG.classification, `0` = "No ARGs",
                                    `1` = "Non-beta-lactamase ARGs",
                                    `2` = "Beta-lactamases")) %>%
    ## remove outlier points with very low coverage.
    filter(CopyNumber > 0.5)

## beta-lactamases have higher copy number compared to other ARGs in these strains.
wilcox.test(beta.lactam.ARGs$CopyNumber, non.beta.lactam.ARGs$CopyNumber)$p.value

## Plasmids with ARGs actually have lower copy numbers than
## plasmids without ARGs.
plasmid.copy.number.data <- chromosome.plasmid.copy.number.data %>%
    filter(SeqType == "plasmid") %>%
    arrange(CopyNumber)

ARG.plasmid.data <- plasmid.copy.number.data %>%
    filter(has.ARG==TRUE)

no.ARG.plasmid.data <- plasmid.copy.number.data %>%
    filter(has.ARG == FALSE)

wilcox.test(ARG.plasmid.data$CopyNumber, no.ARG.plasmid.data$CopyNumber)$p.value

## Supplementary Figure S10.
plasmid.copy.number.plot <- ggplot(plasmid.copy.number.data,
                                   aes(x=log10(CopyNumber),
                                       fill=`Plasmid class`)) +
    geom_histogram(bins=50) +
    theme_classic() +
    xlab("log10(Plasmid copy number)")  +
    theme(legend.position="top")

ggsave("../results/S10Fig-Hawkey2022-plasmid-copy-number.pdf",
       plasmid.copy.number.plot,height=5.75,width=5.75)

################################################################################
## Analysis of chains of duplications, produced by join-duplications.py.
## Look at basic statistics
## for duplicated ARGs and associations with MGE genes,
## and re-calculate statistics for duplicated ARGs associated with MGEs,
## and duplicated ARGs that are not associated with MGEs.

all.joined.duplications <- read.csv("../results/joined-duplicate-proteins.csv") %>%
    tibble() %>%
    ## for numeric consistency, remove all duplications with NA product annotations.
    filter(!is.na(product))


make.ARG.MGE.region.contingency.table <- function(joined.duplications,
                                                  antibiotic.keywords,
                                                  MGE.keywords) {

    ARG.joined.duplications <- joined.duplications %>%
        filter(str_detect(.$product,antibiotic.keywords))

    MGE.joined.duplications <- joined.duplications %>%
        filter(str_detect(.$product, MGE.keywords))
        
    ## get the regions-- drop the sequence information.
    ## There are 773,051 regions in total.
    joined.regions <- joined.duplications %>%
        select(Annotation_Accession, Replicon_Accession, Replicon_type,
               region_index, region_length, num_proteins_in_region, region_start, region_end) %>%
        distinct()
    
    ARG.joined.regions <- ARG.joined.duplications %>%
        select(Annotation_Accession, Replicon_Accession, Replicon_type,
               region_index, region_length, num_proteins_in_region, region_start, region_end) %>%
        distinct()
    
    no.ARG.joined.regions <- anti_join(joined.regions, ARG.joined.regions)
    
    MGE.joined.regions <- MGE.joined.duplications %>%
        select(Annotation_Accession, Replicon_Accession, Replicon_type,
               region_index, region_length, num_proteins_in_region, region_start, region_end) %>%
        distinct()
    
    no.MGE.joined.regions <- anti_join(joined.regions, MGE.joined.regions)
    
    ## count regions (groups) that contain both ARGs and MGE genes.
    ## 3,536 regions contain both ARGs and MGE genes.
    ARG.and.MGE.joined.regions <- inner_join(ARG.joined.regions, MGE.joined.regions)
    
    ## count regions (groups) that contain ARGs but no MGE genes.
    ## 2,551 regions contain ARGs but no MGE genes.
    ARG.and.no.MGE.joined.regions <- inner_join(ARG.joined.regions, no.MGE.joined.regions)
    
    ## count regions (groups) that contain MGE genes but no ARGs.
    ## 624,335 regions contain MGE genes but no ARGs.
    MGE.and.no.ARG.joined.regions <- inner_join(no.ARG.joined.regions, MGE.joined.regions)
    
    ## count regions (groups) that have neither MGE genes nor ARGs.
    ## 142,629 regions contain neither MGE genes nor ARGs.
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

joined.regions.contingency.table <- make.ARG.MGE.region.contingency.table(all.joined.duplications, antibiotic.keywords, MGE.keywords)

fisher.test(joined.regions.contingency.table)

#######################################################################################
## What is the relative contribution of segmental duplications to transpositions for duplicated ARGs in this dataset?
## This analysis shows that segmental duplications play a relatively small role compared to MGEs.

## get all regions containing duplicated ARGs (6,087 of these).
joined.regions.containing.ARGs <- all.joined.duplications %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    select(Annotation_Accession, Replicon_Accession, Replicon_type,
           region_index, region_length, num_proteins_in_region, region_start, region_end) %>%
    distinct() %>%
    ## index these distinct regions containing duplicated ARGs
    mutate(dupARG_region_index = row_number())

## get all associated information, including sequences.
joined.duplications.containing.ARGs <- joined.regions.containing.ARGs %>%
    ## join in this way to keep the dupARG_region_index column.
    left_join(all.joined.duplications)
## write to file.
write.csv(x=joined.duplications.containing.ARGs,
          file="../results/joined-duplications-containing-ARGs.csv")

## from looking at this spreadsheet, several segmental duplications contain several ARGs as well as transposases etc.
## so counting statistics per ARG may not be the best unit for counting independent units.

## let's count the number of these joined duplications containing ARGs that represent segmental duplications,
## and let's get the repeat number of these segmental duplications too.

count.duplicates.within.region <- function(region.df) {
    region.df %>%
        group_by(
            ## general region metadata
            Annotation_Accession, Replicon_Accession, Replicon_type,
            region_index, region_length, num_proteins_in_region, region_start, region_end,
            ## keep the dupARG_region_index column to keep track of what region each ARG is in.
            dupARG_region_index,
            ## actual sequence data
            protein_length, product, sequence) %>%
        summarize(dup.count = n()) %>%
        ungroup()
}


dup.ARG.regions.with.segmental.dup.data <- joined.duplications.containing.ARGs %>%
    split(.$dupARG_region_index) %>%
    map_dfr(count.duplicates.within.region)

## just look at ARGs within these regions with segmental dup data.
dup.ARGs.within.regions.with.segmental.dup.data <- dup.ARG.regions.with.segmental.dup.data %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    select(Annotation_Accession, Replicon_Accession, Replicon_type,
           region_length, num_proteins_in_region, protein_length, product, sequence, dup.count, dupARG_region_index) %>%
    distinct()
## write to file.
write.csv(x=dup.ARGs.within.regions.with.segmental.dup.data,
          file="../results/joined-regions-just-ARGs-and-dup-counts.csv")

## 406 ARGs have multiple copies in the same region of duplicated genes.
segmental.ARG.duplications <- dup.ARGs.within.regions.with.segmental.dup.data %>%
    filter(dup.count > 1)
## These 406 ARGs are found in 237 unique regions of duplicated genes.
length(unique(segmental.ARG.duplications$dupARG_region_index))
## so 237 out of the 6,087 joined regions containing duplicated ARGs have multiple copies of some ARG.


#################################################
## What MGE functions are associated with duplicated ARGs?
## To answer this question, let's look at MGE functions that are jointly duplicated in the same region as
## duplicated ARGs.

## 8,449 genes with MGE-associated functions here.
MGEs.in.joined.duplications.containing.ARGs <- joined.duplications.containing.ARGs %>%
    filter(str_detect(.$product,MGE.keywords))
## write to file.
write.csv(x=MGEs.in.joined.duplications.containing.ARGs,
          file="../results/joined-regions-MGE-functions-associated-with-dupARGs.csv")

## 5,541 transposase sequences here (out of the 8,449).
transposase.in.joined.duplications.containing.ARGs <- joined.duplications.containing.ARGs %>%
    filter(str_detect(.$product,"transposase"))

## 1,046 integrases here (out of the 8,449).
integrase.in.joined.duplications.containing.ARGs <- joined.duplications.containing.ARGs %>%
    filter(str_detect(.$product,"integrase"))

################################################################################
## Examine the frequency of transposase sequences that are found with ARGs, and
## host range of these transposases.

annotated.ARG.associated.transposases <- transposase.in.joined.duplications.containing.ARGs %>%
    left_join(gbk.annotation) %>%
    ## Annotate the genera.
    mutate(Genus = stringr::word(Organism, 1))

## write relevant columns  to file, so that I can cluster these transposons,
## allowing k-mismatches from most common sequence within a cluster, using Julia.
annotated.ARG.associated.transposases %>%
    select(product, Genus, sequence) %>%
    rename(product_annotation = product) %>%
write.csv(file="../results/transposases-in-dup-regions-with-ARGs.csv", quote=FALSE,row.names=FALSE)

## also make a file to cluster E. coli transposases.
Ecoli.ARG.associated.transposases <- annotated.ARG.associated.transposases %>%
    filter(str_detect(.$Organism, "Escherichia coli")) %>%
    group_by(product, sequence) %>%
    summarize(count = n()) %>%
    arrange(desc(count))
## and write it to file.
write.csv(x=Ecoli.ARG.associated.transposases,
          file="../results/Ecoli-transposases-in-dup-regions-with-ARGs.csv")

####################################################################
## now cluster the transposases using the Julia script cluster-transposases.jl.
print("Running: julia cluster-transposases.jl")
system("julia cluster-transposases.jl")
print("Completed: julia cluster-transposases.jl.")
####################################################################

## now read in the clustered transposons made by cluster-transposases.jl.
clustered.ARG.associated.transposases <- read.csv(
    "../results/merged_transposases-in-dup-regions-with-ARGs.csv") %>%
    ## remove "NA" Genera
    filter(!is.na(Genus))

## now let's examine the clustered sequences in order to make a rank order list
## and see the distribution across genus for each sequence.

##make a rank column by total count of the sequence across Genus.
clustered.ARG.associated.transposase.ranks <- clustered.ARG.associated.transposases %>%
    group_by(sequence) %>% summarize(total.count = sum(count)) %>%
    arrange(desc(total.count)) %>%
    mutate(rank = row_number()) %>%
    ungroup()

clustered.ARG.associated.transposases.with.ranks <- clustered.ARG.associated.transposases %>%
    full_join(clustered.ARG.associated.transposase.ranks) %>%
    arrange(rank, desc(count))

full.plot.of.clustered.ARG.associated.transposases <- clustered.ARG.associated.transposases.with.ranks %>%
    ggplot(aes(x=rank, fill=Genus, y=count)) +
    geom_bar(position="stack", stat="identity") +
    theme_classic() +
    theme(legend.position="bottom")

ggsave("../results/S13Fig.pdf",
       full.plot.of.clustered.ARG.associated.transposases,
       width=8, height=7)

## filter on just the top ranks.
top.clustered.ARG.associated.transposases <- clustered.ARG.associated.transposases.with.ranks %>%
    filter(rank <= 10)
## write this out, so that I can make a table of these transposase classes by hand.
write.csv(x=top.clustered.ARG.associated.transposases, file="../results/top-clustered-ARG-associated-transposases.csv")

## by hand, I made a table of rank to IS Class for these top 10 transposase classes,
## by using blastp searches of the top sequence in the cluster against the ISFinder database:
## https://www-is.biotoul.fr/blast.php.
top.transposase.annotations <- read.csv("../data/Top10-ARG-associated-transposase-annotations.csv")
## now add these data as a column.
annotated.top.clustered.ARG.associated.transposases <- full_join(
    top.clustered.ARG.associated.transposases,
    top.transposase.annotations) %>%
    ## order the labels according to their rank.
    mutate(Transposase = factor(Transposase,levels=unique(Transposase), ordered=TRUE))
    

plot.of.top.clustered.ARG.associated.transposases <- annotated.top.clustered.ARG.associated.transposases %>%
    ggplot(aes(x=Transposase, fill=Genus, y=count)) +
    geom_bar(position="stack", stat="identity") +
    theme_classic()  +
    theme(legend.position="bottom", legend.title=element_blank(),
          legend.text = element_text(size = 8)) +
    theme(axis.text.x  = element_text(angle=45,vjust=0.5))

## This will be Figure 5E.
ggsave("../results/top-clustered-ARG-transposons.pdf",
       plot.of.top.clustered.ARG.associated.transposases, height = 5, width=6)

# let's examine just the unique sequences in order to make a rank order list
## and see the distribution across genus for each sequence.

ARG.associated.transposase.sequences <-  annotated.ARG.associated.transposases %>%
    group_by(Genus, sequence) %>% summarize(count.per.genus = n()) %>% arrange(desc(count.per.genus)) %>%
    ## remove rows with missing Genus.
    filter(!is.na(Genus))

##make a rank column by total count of the sequence across Genus.
ARG.associated.transposase.ranks <- annotated.ARG.associated.transposases %>%
    group_by(sequence) %>% summarize(total.count = n()) %>%
    arrange(desc(total.count)) %>%
    mutate(rank = row_number()) %>%
    ungroup()
## and join back to get the Genus information.
ARG.associated.transposase.sequences.with.ranks <- ARG.associated.transposase.sequences %>%
    full_join(ARG.associated.transposase.ranks)

ARG.associated.transposase.rank.plot1 <- ARG.associated.transposase.sequences.with.ranks %>%
    ggplot(aes(x=rank, fill=Genus, y=count.per.genus)) +
    geom_bar(position="stack", stat="identity") +
    theme_classic()

ggsave("../results/ARG-transposon-rank-plot1.pdf",ARG.associated.transposase.rank.plot1)

## filter on just the top ranks.
top.ARG.associated.transposase.sequences.with.ranks <- ARG.associated.transposase.sequences.with.ranks %>%
    filter(rank < 10)

ARG.associated.transposase.rank.plot2 <- top.ARG.associated.transposase.sequences.with.ranks %>%
    ggplot(aes(x=rank, fill=Genus, y=count.per.genus)) +
    geom_bar(position="stack", stat="identity") +
    theme_classic()

ggsave("../results/ARG-transposon-rank-plot2.pdf",ARG.associated.transposase.rank.plot2)

