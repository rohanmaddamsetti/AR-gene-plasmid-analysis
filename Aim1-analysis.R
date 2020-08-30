## Aim1-analysis.R by Rohan Maddamsetti.

## analyse the distribution of AR genes on chromosomes versus plasmids in
## fully-sequenced genomes and plasmids in the NCBI Nucleotide database.

## TODO: weigh hits by the number of genes on the chromosome or the plasmids.
## TODO: expand analysis to additional query sets of antibiotic resistance genes:
## 1) Resfams 2) CARD 3) Carolyn's annotated resistance genes.

## CRITICAL TODO: re-annotate manually-curated-gbk-annotation-table.csv
## using the updated gbk-annotation table. This will add ~300 extra genomes
## to the analysis.

## CRITICAL ANALYSIS TODO: look for evidence of recent diversification.
## In particular look more deeply at the result that AAA+ ATPases,
## and ATPases in general seem to be enriched in gene duplications.

library(tidyverse)
library(cowplot)

## annotate source sequence as plasmid or chromosome.
genome.database <- read.csv("../results/chromosome-plasmid-table.csv")
## This is 7046 strains.
length(unique(genome.database$Annotation_Accession))

gbk.annotation <- as_tibble(read.csv("../data/manually-curated-gbk-annotation-table.csv"))
##  while this is 7047 strains.

## IMPORTANT NOTE: This should be ~7300 genomes after updating the gbk annotation.
## not sure if chromosome-plasmid-table will also be larger.

length(unique(gbk.annotation$Annotation_Accession))
## there's one in gbk.annotation that is missnig from genome.database:
## GCA_900492195.1_T2.26MG-112.21_plasmid
## don't worry about it: it won't affect the analysis anyway.
cds.counts <- read.csv("../results/protein_db_CDS_counts.csv")

## number of host and isolation_source annotations in gbk_annotation.
## 647 unique host annotations.
length(unique(gbk.annotation$host))
## 1763 unique isolation_source annotations.
length(unique(gbk.annotation$isolation_source))

## Filter out genomes where the Manual_Annotation field is NA.
good.gbk.annotation <- gbk.annotation %>%
    filter(!is.na(Manual_Annotation))

length(unique(good.gbk.annotation$Annotation_Accession))
## 5123 strains have good annotation.

protein.db.metadata <- genome.database %>%
    left_join(good.gbk.annotation) %>%
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
antibiotic.keywords <- "lactamase|chloramphenicol|quinolone"

## IMPORTANT TODO: Use the same regular expressions used by Zeevi et al. (2019).
## Transposon: ‘transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon’
## plasmid: ‘relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation’
## phage: ‘capsid|phage|tail|head|tape measure|antiterminatio’
## other HGT mechanisms: ‘integrase|excision\S*|exo- nuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip’
## antibiotic resistance: ‘azole resistance|antibiotic resistance|TetR|tetracycline resist- ance|VanZ|betalactam\S*|beta-lactam|antimicrob\S*|lantibio\S*’.

naive.HGT.data <- read.csv("../results/naive-HGT.csv") %>%
    ## now merge with good gbk annotation.
    ## I am doing a left_join here, because I want the NA Manual_Accessions
    ## in order to predict where these unannotated strains come from.
    left_join(good.gbk.annotation)


## for current analysis, remove any strains with no manual annotation.
## in future work, consider using the NA set to see whether I can predict
## what environment they come from, based on which genes are duplicated.
good.naive.HGT.data <- naive.HGT.data %>%
    filter(!is.na(Manual_Annotation))

## TODO: examine data with and without IS.
IS.removed.from.naive.HGT.data <- naive.HGT.data %>%
    filter(!str_detect(.$product,IS.keywords))
## TODO: examine correlations between IS count and other genes
## within each genome to find genes embedded in transposons.

good.IS.removed.from.naive.HGT.data <- IS.removed.from.naive.HGT.data %>%
    filter(!is.na(Manual_Annotation))

AR.naive.HGT.data <- naive.HGT.data %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    arrange(desc(count))

good.AR.naive.HGT.data <- AR.naive.HGT.data %>%
    filter(!is.na(Manual_Annotation))

## let's look at cases of identical sequences on chromosomes and plasmids.
both.chr.and.plasmid.cases <- good.naive.HGT.data %>%
    filter(chromosome_count >= 1 & plasmid_count >= 1) %>%
    arrange(desc(count))

## make a nested tibble to explore correlations between IS
## and other genes in those genomes.
transposable.groups <- both.chr.and.plasmid.cases %>%
    filter(!is.na(Manual_Annotation)) %>%
    group_by(Annotation_Accession) %>% nest()

just.chromosome.cases <- good.naive.HGT.data %>%
    filter(chromosome_count >= 1 & plasmid_count == 0) %>%
    arrange(desc(count))

just.plasmid.cases <- good.naive.HGT.data %>%
    filter(chromosome_count == 0 & plasmid_count >= 1) %>%
    arrange(desc(count))

##########################################
## Important basic statistics on the data.

## 198,641 duplicated genes, including IS elements and hypothetical proteins.
length(naive.HGT.data$Annotation_Accession)
## 76,974 duplicated genes, excluding IS elements and hypothetical proteins.
length(IS.removed.from.naive.HGT.data$Annotation_Accession)

## calculate the number of isolates in the naive.HGT.data
length(unique(naive.HGT.data$Annotation_Accession))
## 7,369 isolates in the native.HGT.data. This includes strains with
## no good gbk annotation (i.e Manual_Annotation == NA)


## CRITICAL BUG TO FIX:
## there are 7,369 isolates which is more than the annotated genomes.
## in part, this discrepancy has to do with the manual annotation having been
## constructed on an older version of the gbk_annotation.csv, which was
## missing some strains. There are probably other bugs here to fix as well.
## DEBUG THESE ERRORS AND MAKE THESE NUMBERS CONSISTENT!!!
problem.data <- naive.HGT.data %>%
    filter(!(Annotation_Accession %in% genome.database$Annotation_Accession)) %>%
    select(-sequence,-count,-chromosome_count,-plasmid_count,-product) %>%
    distinct()
## there are 871 isolates with Annotation_Accession, but no metadata at all? Why?
length(unique(problem.data$Annotation_Accession))

length(unique(good.naive.HGT.data$Annotation_Accession))

length(unique(IS.removed.from.naive.HGT.data$Annotation_Accession))
## 6,578 isolate after filtered out IS. This includes strains with
## no good gbk annotation (i.e Manual_Annotation == NA)


## 123,539 duplicated genes, including IS elements and hypothetical proteins.
length(good.naive.HGT.data$Annotation_Accession)
## 44,434 duplicated genes, excluding IS elements and hypothetical proteins.
length(good.IS.removed.from.naive.HGT.data$Annotation_Accession)


## calculate the number of isolates in the AR.naive.HGT.data
length(unique(AR.naive.HGT.data$Annotation_Accession))
## 585 isolates with duplicate AR genes.
length(unique(good.AR.naive.HGT.data$Annotation_Accession))
## 311 isolates with duplicate AR genes in the well-annotated set.

###############
## calculate number of isolates in each category.
category.counts <- good.naive.HGT.data %>%
    ## next two lines is to count isolates rather than genes
    select(-count,-chromosome_count,-plasmid_count,-product,-sequence) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(category_count = n()) %>%
    arrange(desc(category_count))

## calculate number of isolates in each category, with AR genes
AR.category.counts <- good.AR.naive.HGT.data %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(-count,-chromosome_count,-plasmid_count,-product,-sequence) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(category_count = n()) %>%
    arrange(desc(category_count))

## Of these, what percentages of ARGs are in the chromosome, plasmid, or both? 
AR.chr.category.counts <- just.chromosome.cases %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(-count,-chromosome_count,-plasmid_count,-product,-sequence) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(category_count = n()) %>%
    arrange(desc(category_count))

AR.plasmid.category.counts <- just.plasmid.cases %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(-count,-chromosome_count,-plasmid_count,-product,-sequence) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(category_count = n()) %>%
    arrange(desc(category_count))

AR.chr.and.plasmid.category.counts <- both.chr.and.plasmid.cases %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    ## next two lines is to count isolates rather than genes
    select(-count,-chromosome_count,-plasmid_count,-product,-sequence) %>%
    distinct() %>%
    group_by(Manual_Annotation) %>%
    summarize(category_count = n()) %>%
    arrange(desc(category_count))

## What's the copy number of each ARG (in the chromosome or in the plasmid)
## in these isolates?

AR.both.chr.and.plasmid.cases <- both.chr.and.plasmid.cases %>%
    filter(str_detect(.$product,antibiotic.keywords))

AR.chr.cases <- just.chromosome.cases %>%
    filter(str_detect(.$product,antibiotic.keywords))

AR.plasmid.cases <- just.plasmid.cases %>%
    filter(str_detect(.$product,antibiotic.keywords))

#########################################################################################
## Figures that Lingchong asked me to make.

## 1A. Stacked bar(or pie chart) of all isolates regardless of duplication
## I think a table is a better way of showing these data.
## make a table in addition.

Fig1A.data <- good.gbk.annotation %>%
    ## adding a dummy column to make a single stacked bar.
    mutate(Isolate="isolate")

Fig1A.with.legend <- ggplot(data=Fig1A.data,
                aes(x=Isolate, fill=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    ggtitle("All isolates")

Fig1A <- Fig1A.with.legend +
    guides(fill=FALSE)

## 1B. Stacked bar for isolates with some sort of duplication

Fig1B.data <- Fig1A.data %>%
    filter(Annotation_Accession %in% good.naive.HGT.data$Annotation_Accession)

Fig1B <- ggplot(data=Fig1B.data,
                aes(x=Isolate, fill=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    ggtitle("Isolates with duplications") +
    guides(fill=FALSE)

## 1C. Stacked bar for isolates with duplications, excluding mobile elements
Fig1C.data <- Fig1A.data %>%
    filter(Annotation_Accession %in% good.IS.removed.from.naive.HGT.data$Annotation_Accession)

Fig1C <- ggplot(data=Fig1C.data,
                aes(x=Isolate, fill=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    ggtitle("Isolates with duplications, excluding MGEs") +
    guides(fill=FALSE)

## 1D. Stacked bar for isolates with duplication of ARGs
Fig1D.data <- Fig1A.data %>%
    filter(Annotation_Accession %in% good.AR.naive.HGT.data$Annotation_Accession)

Fig1D <- ggplot(data=Fig1D.data,
                aes(x=Isolate, fill=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    ggtitle("Isolates with duplicated ARGs") +
    guides(fill=FALSE)

Fig1.legend <- cowplot::get_legend(Fig1A.with.legend)
Fig1.panels <- plot_grid(Fig1A,Fig1B,Fig1C,Fig1D,
                  labels=c('A','B','C','D'),ncol=1)

## By changing the y-axis scale for each plot while maintaining its height,
## Fig. 1 also shows the following information: 
## Bar charts for percentage of duplication of any gene
## Bar charts for percentage of duplication of any gene, excluding MGE.
## Bar charts for percentage of duplication for ARGs
Fig1 <- plot_grid(Fig1.panels,Fig1.legend,ncol=2)
ggsave("../results/Fig1.pdf",Fig1)

## Fig2A. Average copy number of duplicated genes

Fig2A.data <- good.naive.HGT.data %>%
    group_by(Manual_Annotation) %>%
    summarize(average.copy.num =  mean(count),
              average.chromosome.copies = mean(chromosome_count),
              average.plasmid.copies = mean(plasmid_count))
## Fig2B. Average copy number of duplicated ARG genes

Fig2B.data <- good.AR.naive.HGT.data %>%
    group_by(Manual_Annotation) %>%
    summarize(average.AR.copy.num =  mean(count),
              average.AR.chromosome.copies = mean(chromosome_count),
              average.AR.plasmid.copies = mean(plasmid_count))

## Fig3. Show average number of genes on just chromosome, just plasmid, or on both,
## per category (normalize by number of genomes in each category).

Fig3A.data <- just.chromosome.cases %>%
    group_by(Manual_Annotation) %>%
    summarize(average.just.chromosome.copy.num =  mean(count))

Fig3B.data <- just.plasmid.cases %>%
    group_by(Manual_Annotation) %>%
    summarize(average.just.plasmid.copy.num =  mean(count))

Fig3C.data <- both.chr.and.plasmid.cases %>%
    group_by(Manual_Annotation) %>%
    summarize(average.on.both.chromosome.and.plasmid.copy.num =  mean(count),
              average.chromosome.copies = mean(chromosome_count),
              average.plasmid.copies = mean(plasmid_count))

## Fig4. Show average number of genes on just chromosome, just plasmid, or on both,
## per category (normalize by number of genomes in each category).
## Exclude MGEs.

Fig4A.data <- just.chromosome.cases %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(average.just.chromosome.copy.num =  mean(count))

Fig4B.data <- just.plasmid.cases %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(average.just.plasmid.copy.num =  mean(count))

Fig4C.data <- both.chr.and.plasmid.cases %>%
    filter(!str_detect(.$product,IS.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(average.on.both.chromosome.and.plasmid.copy.num =  mean(count),
              average.chromosome.copies = mean(chromosome_count),
              average.plasmid.copies = mean(plasmid_count))


## Fig5. Show average number of AR genes on just chromosome, just plasmid, or on both,
## per category (normalize by number of genomes in each category).

Fig5A.data <- just.chromosome.cases %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(average.just.chromosome.copy.num =  mean(count))

Fig5B.data <- just.plasmid.cases %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(average.just.plasmid.copy.num =  mean(count))

Fig5C.data <- both.chr.and.plasmid.cases %>%
    filter(str_detect(.$product,antibiotic.keywords)) %>%
    group_by(Manual_Annotation) %>%
    summarize(average.on.both.chromosome.and.plasmid.copy.num =  mean(count),
              average.chromosome.copies = mean(chromosome_count),
              average.plasmid.copies = mean(plasmid_count))




#########################################################################################
## older plots.

## TODO: make count a column, and add a row for chromosome and plasmid.
## Then facet on this.
## THIS IS COOL! Some genomes have very high-copy number duplications.
both.chr.plasmid.distribution.plot <- ggplot(both.chr.and.plasmid.cases,
                                             aes(x=chromosome_count,y=plasmid_count,
                                                 color=Manual_Annotation)) +
    geom_jitter() +
    geom_rug() +
    facet_wrap(.~Manual_Annotation) +
    theme_classic() +
    ggtitle("both chromosome and plasmid")

AR.both.chr.plasmid.distribution.plot <- ggplot(AR.both.chr.and.plasmid.cases,
                                                aes(x=chromosome_count,y=plasmid_count,
                                                    color=Manual_Annotation)) +
    geom_jitter() +
    geom_rug() +
    facet_wrap(.~Manual_Annotation) +
    theme_classic() +
    ggtitle("both chromosome and plasmid")

AR.chr.distribution.plot <- ggplot(AR.chr.cases,
                                   aes(x=chromosome_count,
                                       fill=Manual_Annotation)) +
    geom_histogram() +
    theme_classic() +
    ggtitle("just chromosome")

AR.plasmid.distribution.plot <- ggplot(AR.plasmid.cases,
                                       aes(x=plasmid_count,
                                           fill=Manual_Annotation)) +
    geom_histogram() +
    theme_classic() +
    ggtitle("just plasmid")

############################################

both.chr.plasmid.metadata.plot <- ggplot(both.chr.and.plasmid.cases,
                                            aes(x=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill=FALSE) +
    ggtitle("both chromosome and plasmid")

just.chromosome.metadata.plot <- ggplot(just.chromosome.cases,
                                            aes(x=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill=FALSE) +
    ggtitle("just on chromosome")

just.plasmid.metadata.plot <- ggplot(just.plasmid.cases,
                                            aes(x=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill=FALSE) +
    ggtitle("just on plasmid")

full.metadata.panels <- plot_grid(just.chromosome.metadata.plot,
                                   just.plasmid.metadata.plot,
                                   both.chr.plasmid.metadata.plot,
                                   labels=c('A','B','C'),
                                   ncol=1)

full.metadata.title <- ggdraw() + 
  draw_label(
    "location of recently duplicated genes",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

full.metadata.plot <- plot_grid(full.metadata.title,
                                full.metadata.panels,
                                ncol=1,
                                ## rel_heights values control vertical title margins
                                rel_heights = c(0.1, 1))

ggsave("../results/all-recent-duplicates-metadata.pdf",full.metadata.plot,height=12,width=7)
##############
## now just examine genes with the antibiotic keywords that I chose above.
both.chr.and.plasmid.antibiotic.cases <- both.chr.and.plasmid.cases %>%
    filter(str_detect(.$product,antibiotic.keywords))

just.chromosome.antibiotic.cases <- just.chromosome.cases %>%
    filter(str_detect(.$product,antibiotic.keywords))

just.plasmid.antibiotic.cases <- just.plasmid.cases %>%
    filter(str_detect(.$product,antibiotic.keywords))


AR.both.chr.plasmid.metadata.plot <- ggplot(both.chr.and.plasmid.antibiotic.cases,
                                            aes(x=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill=FALSE) +
    ylim(0,120) +
    ggtitle("both chromosome and plasmid")

AR.just.chromosome.metadata.plot <- ggplot(just.chromosome.antibiotic.cases,
                                            aes(x=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill=FALSE) +
    ylim(0,120) +
    ggtitle("just on chromosome")

AR.just.plasmid.metadata.plot <- ggplot(just.plasmid.antibiotic.cases,
                                            aes(x=Manual_Annotation)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill=FALSE) +
    ggtitle("just on plasmid")

full.AR.metadata.panels <- plot_grid(AR.just.chromosome.metadata.plot,
                                   AR.just.plasmid.metadata.plot,
                                   AR.both.chr.plasmid.metadata.plot,
                                   labels=c('A','B','C'),
                                   nrow=1)

full.AR.metadata.title <- ggdraw() + 
  draw_label(
    "location of recently duplicated antibiotic resistance genes",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

full.AR.metadata.plot <- plot_grid(full.AR.metadata.title,
                                   full.AR.metadata.panels,
                                   ncol=1,
                                   ## rel_heights values control vertical title margins
                                   rel_heights = c(0.1, 1))

ggsave("../results/AR-recent-duplicates-metadata.pdf",full.AR.metadata.plot,height=7,width=12)

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


    ## TODO: FIGURE OUT AN APPROPRIATE NORMALIZATION BY NUMBER
    ## OF GENES ON PLASMIDS VERSUS THE CHROMOSOME.
    
    ## now weight by fraction of genes on plasmids.
    ##print("Weighted mean fraction of AR genes on plasmids, in human host isolates:")
    ##print(mean(human.host.combined.df$weighted.percent.on.plasmid))
    ##print("Mean fraction of AR genes on plasmids, in isolates without a host:")
    ##print(mean(no.host.combined.df$weighted.percent.on.plasmid))

    ## the difference in percentage on plasmid is NOT significantly different:
    ## Mann-Whitney U-test: p = 0.06082
    ##print("Mann-Whitney U-test 2: do isolates from human hosts have a larger fraction of AR genes on plasmids, weighted by number of CDS of plasmids?")
    ##test.result2 <- wilcox.test(human.host.combined.df$weighted.percent.on.plasmid,
    ##                            no.host.combined.df$weighted.percent.on.plasmid,
   ##                             alternative="greater")
    ##print(test.result2)

   ##     print("Mann-Whitney U-test 3: do isolates from human hosts have a larger fraction of AR genes on chromosomes")
   ## test.result2 <- wilcox.test(human.host.combined.df$weighted.percent.on.plasmid,
  ##                              no.host.combined.df$weighted.percent.on.plasmid,
 ##                               alternative="less")
    ##   print(test.result3)
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
