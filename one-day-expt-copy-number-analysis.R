##one-day-expt-copy-number-analysis.R by Rohan Maddamsetti.

## 1) use xml2 to get negative binomial fit from
## breseq output summary.html. This is H0 null distribution of 1x coverage.

## 2) Find intervals longer than max.read.len that reject H0 coverage in genome.
##    at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.

## 3) Do a more rigorous test for each region. take positions in the region separated by more than max.read.len,
## and determine the probability that all are independently significant under the null, compared to
## a corrected bonferroni. The max.read.len ensures positions cannot be spanned by a single Illumina read.

## 4) Estimate copy number by dividing mean coverage in each region by the mean
##   of the H0 1x coverage distribution.

## 5) return copy number and boundaries for each significant amplification.

library(tidyverse)
library(xml2)
library(cowplot)
library(data.table)
library(dtplyr)
library(assertthat)
library(viridis)

## Bioconductor dependencies
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

#' parse the summary.html breseq output file, and return the mean and relative variance
#' of the negative binomial fit to the read coverage distribution, returned as a
#' data.frame with columns {mean, relative.variance}.
#' NOTE: this code has only been tested on the summary file
#' output by breseq 0.35.0. It will fail on breseq 0.37 and later, which uses the term "relative variance".

coverage.nbinom.from.html <- function(breseq.output.dir, sample.has.plasmid=TRUE) {
    summary.html.f <- file.path(breseq.output.dir, "output", "summary.html")
    tree <- read_html(summary.html.f)
    ## print text in the table 'Reference Sequence Information.
    query <- '//table[./tr/th[contains(text(),"fit dispersion")]]'
    table <- xml_find_first(tree,query)
    table.data <- xml_find_all(table,'./tr/td')
    chromosome.avg <- as.numeric(xml_text(table.data[5]))
    chromosome.relative.variance <- as.numeric(xml_text(table.data[6]))
    ## all samples should have these data.
    coverage.df <- data.frame('Sample' = basename(breseq.output.dir),
                              'mean'=c(chromosome.avg),
                              'relative.variance'=c(chromosome.relative.variance),
                              'variance'=c(chromosome.avg * chromosome.relative.variance),
                              'replicon'=c("chromosome"))
    if (sample.has.plasmid) {
            plasmid.avg <- as.numeric(xml_text(table.data[21]))
            plasmid.relative.variance <- as.numeric(xml_text(table.data[22]))
            plasmid.coverage.df <- data.frame('Sample' = basename(breseq.output.dir),
                                              'mean' = plasmid.avg,
                                              'relative.variance' = plasmid.relative.variance,
                                              'variance' = plasmid.avg * plasmid.relative.variance,
                                              'replicon' = "plasmid")
            ## now join the plasmid coverage data.
            coverage.df <- rbind(coverage.df, plasmid.coverage.df)
    }
    return(coverage.df)
}

#' get the maximum length of a sequencing read from the summary.html breseq
#' output file.
max.readlen.from.html <- function(breseq.output.dir) {
    summary.html.f <- file.path(breseq.output.dir, "output", "summary.html")
    tree <- read_html(summary.html.f)
    ## print text in the table 'Read File Information.
    query <- '//table[./tr/th[contains(text(),"longest")]]'
    table <- xml_find_first(tree,query)
    table.data <- xml_find_all(table,'./tr/td')
    readlen.index <- length(table.data) - 1
    max.readlen <- xml_integer(xml_find_all(table.data[readlen.index],".//b//text()"))
    return(max.readlen)
}


#' Find intervals longer than max.read.len that reject H0 coverage in genome.
#' at an uncorrected alpha = 0.05. This is to have generous predicted boundaries for amplifications.
#' Then do a more rigorous test for each region. take positions in the region separated by more than max.read.len,
#' and determine the probability that all are independently significant under the null, compared to
#' a corrected bonferroni. max.read.len ensures positions cannot be spanned by a single Illumina read.
#' Estimate copy number by dividing mean coverage in each region by the mean of the H0 1x coverage distribution.
#' return mean copy number, and boundaries for each region that passes the amplification test.

find.K12.candidate.amplifications <- function(breseq.output.dir, gnome) { #gnome is not a misspelling.

    gnome <- as.character(gnome)
    print(gnome)
    ## Use xml2 to get negative binomial fit and relative variance from
    ## breseq output summary.html. This is H0 null distribution of 1x coverage.
    nbinom.fit <- coverage.nbinom.from.html(breseq.output.dir) %>%
        filter(replicon=="chromosome")
    
    ## Use xml2 to get max read length from summary.html.
    max.read.len <- max.readlen.from.html(breseq.output.dir)
    genome.length <- 4641652 ## length of K-12 MG1655 reference.
    my.size.parameter <- nbinom.fit$mean^2/(nbinom.fit$variance - nbinom.fit$mean)
    
    alpha <- 0.05
   
    uncorrected.threshold <- qnbinom(p=alpha, mu=nbinom.fit$mean, size=my.size.parameter, lower.tail=FALSE)
    
    genome.coverage.file <- file.path(breseq.output.dir,"08_mutation_identification", "NC_000913.coverage.tab")
    
    ## use dtplyr for speed!
    genome.coverage <- lazy_dt(fread(genome.coverage.file)) %>%
        select(position,unique_top_cov,unique_bot_cov) %>% mutate(coverage=unique_top_cov+unique_bot_cov)
    
    ## find candidate amplifications that pass the uncorrected threshold.
    candidate.amplifications <- genome.coverage %>%
        filter(coverage > uncorrected.threshold) %>%
        ## now finally turn into a tibble.
        as_tibble()
    
    ## calculate intervals of candidate amplifications.
    boundaries <- candidate.amplifications %>%
        mutate(left.diff=position - lag(position)) %>%
        mutate(right.diff=lead(position) - position) %>%
        ## corner case: check for the NA values at the endpoints and set them as boundaries.
        mutate(is.right.boundary=is.na(right.diff)|ifelse(right.diff>1,TRUE,FALSE)) %>%
        mutate(is.left.boundary=is.na(left.diff)|ifelse(left.diff>1,TRUE,FALSE)) %>%
        filter(is.left.boundary==TRUE | is.right.boundary==TRUE)
 
    left.boundaries <- filter(boundaries,is.left.boundary==TRUE) %>%
        arrange(position)
        
    right.boundaries <- filter(boundaries,is.right.boundary==TRUE) %>%
        arrange(position)
    
    assert_that(nrow(left.boundaries) == nrow(right.boundaries))
    
    ## helper higher-order function to get min, max, mean coverage of each segment.
    get.segment.coverage <- function(left.bound,right.bound,coverage.table,funcx) {
        seg <- coverage.table %>% filter(position>left.bound) %>% filter(position<right.bound)
        return(funcx(seg$coverage))
    }
    
    amplified.segments <- data.frame(left.boundary=left.boundaries$position,right.boundary=right.boundaries$position) %>%
        ## filter out intervals less than 2 * max.read.len.
        mutate(len=right.boundary-left.boundary) %>%
    filter(len>(2*max.read.len)) %>% mutate(amplication.index=row_number())

    ## return empty dataframe  if there are no significant amplified segments.
    if (nrow(amplified.segments) == 0) return(data.frame())

    amplified.segments <- amplified.segments %>%
        ## find min, max, and mean coverage of each amplified segment.
        group_by(left.boundary,right.boundary) %>%
        summarise(coverage.min=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,min),
                  coverage.max=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,max),
                  coverage.mean=get.segment.coverage(left.boundary,right.boundary,candidate.amplifications,mean)) %>%
        mutate(len=right.boundary-left.boundary) %>%
        mutate(copy.number.min=coverage.min/nbinom.fit$mean,copy.number.max=coverage.max/nbinom.fit$mean,
               copy.number.mean=coverage.mean/nbinom.fit$mean) %>%
        ## annotate with the sample name.
        mutate(Sample=as.character(gnome))

    return(amplified.segments)
}


find.K12.chromosomal.amplifications <- function(breseq.output.dir, gnome) { #gnome is not a misspelling.
  
    amplified.segments <- find.K12.candidate.amplifications(breseq.output.dir, gnome)
    ## handle the case that there are no amplified.segments (empty dataframe).
    if (nrow(amplified.segments) == 0) return(amplified.segments)
    
    ## Use xml2 to get negative binomial fit and relative variance from
    ## breseq output summary.html. This is H0 null distribution of 1x coverage.
    nbinom.fit <- coverage.nbinom.from.html(breseq.output.dir) %>%
        filter(replicon=="chromosome")  
    ## Use xml2 to get max read length from summary.html.
    max.read.len <- max.readlen.from.html(breseq.output.dir)
    genome.length <- 4641652 ## length of K-12 MG1655 reference.
    my.size.parameter <- nbinom.fit$mean^2/(nbinom.fit$variance - nbinom.fit$mean)
    alpha <- 0.05

    ## divide alpha by the number of tests for the bonferroni correction.
    bonferroni.alpha <- alpha/(genome.length + sum(amplified.segments$len))
    corrected.threshold <- qnbinom(p = bonferroni.alpha, mu = nbinom.fit$mean, size = my.size.parameter, lower.tail=FALSE)
    
    ## This is my test: take the probability of the minimum coverage under H0 to the power of the number of
    ## uncorrelated sites in the amplification (sites more than max.read.len apart). Then see if this is smaller than the
    ## bonferroni corrected p-value for significance..
    significant.amplifications <- amplified.segments %>%
        mutate(pval=(pnbinom(q = coverage.min,
                             mu = nbinom.fit$mean,
                             size = my.size.parameter,
                             lower.tail=FALSE))^(len%/%max.read.len)) %>%
        mutate(is.significant=ifelse(pval < bonferroni.alpha,TRUE,FALSE)) %>%
        filter(is.significant==TRUE) %>%
        mutate(bonferroni.corrected.pval=pval*alpha/bonferroni.alpha)
    
    return(significant.amplifications)
}


annotate.sample.amplifications <- function(sample.amplifications) {

    ancestor.gff <- unique(sample.amplifications$gff_path)
    
    ## create the IRanges object.
    amp.ranges <- IRanges(sample.amplifications$left.boundary,
                          sample.amplifications$right.boundary)
    ## Turn into a GRanges object in order to find overlaps with K12 genes.
    g.amp.ranges <- GRanges("NC_000913", ranges=amp.ranges)
    ## and add the data.frame of sample.amplifications as metadata.
    mcols(g.amp.ranges) <- sample.amplifications
    
    ## find the genes within the amplifications.
    ancestor.gff.data <- import.gff(ancestor.gff)
    ancestor.Granges <- as(ancestor.gff.data, "GRanges")
    
    ancestor.genes <- ancestor.Granges[ancestor.Granges$type == 'gene']
    ## find overlaps between annotated genes and amplifications.
    hits <- findOverlaps(ancestor.genes,g.amp.ranges,ignore.strand=FALSE)
    
    ## take the hits, the ancestor annotation, and the amplifications,
    ## and produce a table of genes found in each amplication.
    
    hits.df <- data.frame(query.index=queryHits(hits),subject.index=subjectHits(hits))
    
    query.df <- data.frame(query.index=seq_len(length(ancestor.genes)),
                           gene=ancestor.genes$Name,locus_tag=ancestor.genes$ID,
                           start=start(ranges(ancestor.genes)),end=end(ranges(ancestor.genes)))
    
    subject.df <- bind_cols(data.frame(subject.index=seq_len(length(g.amp.ranges))),data.frame(mcols(g.amp.ranges)))
    
    amplified.genes.df <- left_join(hits.df,query.df) %>% left_join(subject.df) %>%
        ## if gene is NA, replace with locus_tag. have to change factors to strings!
        mutate(gene = ifelse(is.na(gene),as.character(locus_tag),as.character(gene)))
    
    return(amplified.genes.df)
}


annotate.amplifications <- function(amps.with.ancestors) {
    amps.with.ancestors %>% split(.$Sample) %>%
        map_dfr(.f = annotate.sample.amplifications)    
}


plot.amp.segments <- function(annotated.amps) {
    
    labeled.annotated.amps <- annotated.amps %>%
        mutate(log2.copy.number.mean=log2(copy.number.mean)) %>%
        mutate(left.boundary.MB = left.boundary/1000000) %>%
        mutate(right.boundary.MB = right.boundary/1000000)
    
    ## order the genes by start to get axes correct on heatmap.
    labeled.annotated.amps$gene <- with(labeled.annotated.amps, reorder(gene, start))
    ## reverse the order of genomes to make axes consistent with stacked barplot.
    labeled.annotated.amps$Sample <- factor(labeled.annotated.amps$Sample)
    labeled.annotated.amps$Sample <- factor(labeled.annotated.amps$Sample,
                                            levels=rev(levels(labeled.annotated.amps$Sample)))
    
    segmentplot <- ggplot(
        labeled.annotated.amps,
        aes(x=left.boundary.MB,
            xend=right.boundary.MB,
            y=Sample,
            yend=Sample,
            color=copy.number.mean,
            size=20,
            frame=Plasmid)) +
        geom_segment() +
        ## draw vertical lines at acrBAR starts.
        geom_vline(size=0.2,
                   color = 'red',
                   linetype = 'dashed',
                   xintercept = c(480478/1000000, 484843/1000000, 484985/1000000)) +
        ## draw lines at boundaries of the MD1 deletion in "Engineering a Reduced Escherichia coli Genome"
        geom_vline(size=0.2,
                   color = 'purple',
                   linetype = 'dashed',
                   xintercept = c(262914/1000000, 324588/1000000)) + 
        ## draw vertical lines at boundaries of MD12 deletion
        ## in Table 1 of the the paper "Engineering a Reduced Escherichia coli Genome"
        geom_vline(size=0.2,
                   color = 'green',
                   linetype = 'dashed',
                   ##  coordinates of cryptic prophage DL12 (MD12) deletion.
                   xintercept = c(563979/1000000, 585280/1000000)) +
        xlab("Genomic position (Mb)") +
        ylab("") +
        scale_color_viridis(name="copy number",option="plasma") +
        facet_wrap(~Plasmid,nrow=2, scales = "free_y") +
        theme_classic(base_family='Helvetica') +
        guides(size= "none") +
        theme(legend.position="bottom") +
        theme(axis.ticks=element_line(size=0.1))
    return(segmentplot)
}

#######################################
## Analysis time!

## assert that we are in the src directory, such that
## proj.dir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("ARG-duplications","src")))
projdir <- file.path("..")

## get metadata for all the evolved population metagenomes.
metagenome.metadata <- read.csv("../data/one-day-expt-evolved-sample-metadata.csv")

mixedpop.output.dir <- file.path(projdir, "results", "one-day-expt-genome-analysis", "one-day-expt-evolved-pops")
all.mixedpops <- list.files(mixedpop.output.dir,pattern='^RM')
all.mixedpop.paths <- sapply(all.mixedpops, function(x) file.path(mixedpop.output.dir,x))
mixedpop.input.df <- data.frame(Sample=all.mixedpops, path=all.mixedpop.paths) %>%
    ## skip the two clone samples for now.
    inner_join(metagenome.metadata) %>%
    ## skip the pUC samples.
    filter(Plasmid != "pUC")

## get metadata for the ancestral clones
ancestralclone.metadata <- read.csv("../data/one-day-expt-ancestral-sample-metadata.csv") %>%
    ##skip the pUC samples.
    filter(Plasmid != "pUC")

## get corresponding inputs for the ancestral clones.
ancestralclone.output.dir <- file.path(projdir, "results", "one-day-expt-genome-analysis")
all.ancestralclones <- list.files(ancestralclone.output.dir, pattern='^remapped-RM')
## we have to cut the "remapped-" prefix from the sample names.
remapped.prefix.len <- 10
all.ancestralclone.samples <- substring(all.ancestralclones, remapped.prefix.len)

all.ancestralclone.paths <- sapply(all.ancestralclones, function(x) file.path(ancestralclone.output.dir,x))

ancestralclone.input.df <- data.frame(
    Sample=all.ancestralclone.samples,
    path=all.ancestralclone.paths) %>%
    inner_join(ancestralclone.metadata)

ancestral.clones.df <- ancestralclone.metadata %>%
    ## these clones are the ancestors, so add a Ancestor column to sample ID.
    ## IMPORTANT: ancestral.clones.df is needed downstream is some places,
    ## in other places ancestralclone.metadata is what is needed.
   dplyr::rename(Ancestor = Sample) %>%
    select(-SampleType) %>%
    mutate(gff_name = paste0(Ancestor, ".gff3")) %>%
    mutate(gff_path = file.path(projdir, "results", "one-day-expt-genome-analysis", gff_name))

######################################################################
## Plot the plasmid/chromosome and transposon/chromosome ratio in each sample.

## Get the actual coverage for the B30 and B59 transposons.
## This is calculated by get-one-day-expt-transposon-coverage.py.
transposon.coverage.file <- file.path(projdir, "results", "one-day-expt-genome-analysis", "transposon-coverage.csv")
transposon.coverage.df <- read.csv(transposon.coverage.file) %>%
    ## let's add metadata and rename columns for compatibility.
    dplyr::rename(mean = TransposonCoverage) %>%
    dplyr::mutate(replicon = "transposon")

## data for the ancestral samples.
ancestral.transposon.coverage.df <- ancestralclone.metadata %>%
    left_join(transposon.coverage.df)

ancestral.replicon.coverage.df <- map_dfr(.x = ancestralclone.input.df$path, .f = coverage.nbinom.from.html) %>%
    ## fix the names of the samples.
    mutate(Sample = substring(Sample, remapped.prefix.len)) %>%
    ## add metadata.
    full_join(ancestralclone.metadata) %>%
    ## select the columns that matter in this analysis.
    select(Sample, SampleType, Transposon, Plasmid, replicon, mean) %>%
    ## add rows for the transposon coverage data.
    full_join(ancestral.transposon.coverage.df) %>%
    ## set sensible default values for these columns
    ## for comparison to the evolved populations.
    mutate(Population = 1) %>%
    mutate(Tet = 0)

ancestral.replicon.coverage.ratio.df <- ancestral.replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Transposon, Plasmid, Population, Tet) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome),
              plasmids.per.chromosome = (mean_plasmid/mean_chromosome),
              transposons.per.plasmid = (mean_transposon/mean_plasmid)) %>%
    pivot_longer(cols = c(transposons.per.chromosome,plasmids.per.chromosome,transposons.per.plasmid),
                 names_to = "ratio_type", values_to = "ratio") %>%
    ## add a Day column to compare to evolved samples.
    mutate(Day = 0)

## data for the evolved samples.
evolved.transposon.coverage.df <- metagenome.metadata %>%
    left_join(transposon.coverage.df) %>%
    ## IMPORTANT: remove the pUC samples.
    filter(Plasmid != "pUC") %>%
    ## we don't need this when joining to evolved.replicon.coverage.df.
    select(-SampleType)

evolved.replicon.coverage.df <- map_dfr(.x = mixedpop.input.df$path, .f = coverage.nbinom.from.html) %>%
    inner_join(metagenome.metadata) %>%
    ## I am not examining dispersion or variance at this point.
    select(Sample, mean, replicon, Transposon, Plasmid, Population, Tet) %>%
    ## add rows for the transposon coverage data.
    full_join(evolved.transposon.coverage.df)    

evolved.replicon.coverage.ratio.df <- evolved.replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Transposon, Plasmid, Population, Tet) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome),
              plasmids.per.chromosome = (mean_plasmid/mean_chromosome),
              transposons.per.plasmid = (mean_transposon/mean_plasmid)) %>%
    ungroup() %>%
    pivot_longer(cols = c(transposons.per.chromosome,plasmids.per.chromosome,transposons.per.plasmid),
                 names_to = "ratio_type", values_to = "ratio") %>%
    ## add a Day column to compare to ancestral samples.
    mutate(Day = 1)

Tet5.ratio.plot <- evolved.replicon.coverage.ratio.df %>%
    filter(Tet == 5) %>%
    ggplot(aes(y = ratio, x = Plasmid, color = ratio_type, shape = Transposon)) +
    geom_point() +
    theme_classic() +
    facet_wrap(ratio_type~Transposon, scales = "free") +
    theme(strip.background = element_blank()) +
    ggtitle("5 ug/mL tetracycline, Day 1") +
    guides(color= "none", shape = "none")

Tet0.ratio.plot <- evolved.replicon.coverage.ratio.df %>%
    filter(Tet == 0) %>%
    ggplot(aes(y = ratio, x = Plasmid, color = ratio_type, shape = Transposon)) +
    geom_point() +
    theme_classic() +
    facet_wrap(ratio_type~Transposon, scales = "free") +
    theme(strip.background = element_blank()) +
    ggtitle("0 ug/mL tetracycline, Day 1") +
    guides(color = "none", shape = "none")

ratio.figure <- plot_grid(Tet5.ratio.plot, Tet0.ratio.plot, labels=c('A','B'),nrow=2)
## The warning messages in saving the plot correspond to the missing ratios
## for the 'No plasmid' strains for which transposon.per.plasmid
## and plasmid.per.chromosome are undefined (NA).
ggsave("../results/one-day-expt-coverage-ratios.pdf", ratio.figure)

## let's write out the table too.
write.csv(evolved.replicon.coverage.ratio.df, "../results/one-day-expt-plasmid-transposon-coverage-ratios.csv",
          quote=F, row.names=FALSE)

################################################################################
## Figure 2A.

## This is the big data frame for making Figure 2A.
ancestral.and.evolved.replicon.coverage.ratio.df <- full_join(
    ancestral.replicon.coverage.ratio.df, evolved.replicon.coverage.ratio.df) %>%
    ungroup()

Fig2A.df <- ancestral.and.evolved.replicon.coverage.ratio.df %>%
    filter(ratio_type == "transposons.per.chromosome") %>%
    ## update the names of the Plasmid factor for a prettier plot.
    mutate(Plasmid = fct_recode(as.factor(Plasmid),
                                `No plasmid` = "None",
                                p15A = "p15A")) %>%
    ## update the names of the Transposon factor for a prettier plot.
    mutate(Transposon = fct_recode(as.factor(Transposon),
                                   Inactive = "B59",
                                   Active = "B30")) %>%
    ## turn Tet, Population, Day into discrete factors for plotting.
    mutate(Tet = as.factor(Tet)) %>%
    mutate(Population = as.factor(Population)) %>%
    mutate(Day = as.factor(Day))

## Write Figure 2A Source Data.
write.csv(Fig2A.df, "../results/Source-Data/Fig2A-Source-Data.csv", row.names=FALSE, quote=FALSE)

## Make Figure 2A.
Fig2A <- ggplot(Fig2A.df, aes(x = Day, y = ratio, color = Transposon, shape = Tet)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    theme(strip.background = element_blank()) +
    ylab("tetA-transposons per chromosome")
ggsave("../results/Fig2A.pdf", Fig2A, width=4.5, height=3.25)

########################################################
## Heatmap figure for examining amplifications.

## Find chromosomal amplifications in B59 Tet 5 samples, and annotate with their ancestor gff file.
B59.Tet5.mixedpops <- mixedpop.input.df %>%
    filter(Transposon == "B59") %>%
    filter(Tet == 5)

amps.with.ancestors <- map2_df(B59.Tet5.mixedpops$path,
                               B59.Tet5.mixedpops$Sample,                              
                               ##find.K12.candidate.amplifications) %>% ## for uncorrected p-values
                               find.K12.chromosomal.amplifications) %>% ## for corrected p-values
    ungroup() %>%
    left_join(metagenome.metadata) %>%
    select(-SampleType) %>%
    left_join(ancestral.clones.df)

annotated.amps <- annotate.amplifications(amps.with.ancestors)

## make a figure of the significant amplifications found with the bonferroni method.
amp.segment.plot <- plot.amp.segments(annotated.amps)
## show the plot.
ggsave("../results/B59-amplification-heatmap.pdf", amp.segment.plot)

parallel.amplified.genes <- annotated.amps %>%
    group_by(gene, locus_tag, start, end) %>%
    summarize(parallel.amplifications = n()) %>%
    arrange(desc(parallel.amplifications)) %>%
    filter(parallel.amplifications > 2)

## yaiT is a pseudogene in K12 that is found in fragments.
yaiT.amps <- annotated.amps %>%
    filter(gene == "yaiT")

acrABR.amps <- annotated.amps %>%
    filter(str_detect(gene, "acr"))

################################################################################
## Figure 2D.
## Make a figure of normalized chromosomal coverage per gene in E. coli K12+B59 samples.
## show K12+B30 samples too, as an additional control.
##These samples do not have an active transposase, and there is evidence of
##copy number variation in a very large region surround the native efflux pump operon
##acrAB.

## for each sample, calculate mean coverage per protein-coding gene.
get.protein.coverage.per.sample <- function(row.df) { 

    breseq.output.dir <- row.df$path
    gnome <- row.df$Sample #gnome is not a misspelling.
    
    gnome <- as.character(gnome)
    print(gnome)
   
    genome.coverage.file <- file.path(breseq.output.dir, "08_mutation_identification", "NC_000913.coverage.tab")
    ## use dtplyr for fast input.
    genome.coverage.df <- fread(genome.coverage.file) %>%
        lazy_dt() %>%
        select(position,unique_top_cov,unique_bot_cov) %>%
        mutate(coverage=unique_top_cov+unique_bot_cov) %>%
        ## immediately turn into a data.frame when we are done.
        ## just reading the data into memory is a major bottleneck.
        as.data.frame()

    get.coverage.per.protein <- function(my.protein) {        
        relevant.coverage <- genome.coverage.df %>%
            dplyr::filter(position >= my.protein$start) %>%
            dplyr::filter(position <= my.protein$end)
        
        mean.coverage <- base::mean(relevant.coverage$coverage) 
        my.protein.coverage.df <- mutate(my.protein, mean_coverage=mean.coverage)
        return(my.protein.coverage.df)
    }

    ## get the proteins in the reference genome.
    reference.gff <- row.df$gff_path
    ref.gff.data <- import.gff(reference.gff)
    ref.proteins <- ref.gff.data[ref.gff.data$type == 'CDS']
    ## AFAIK Granges and GrangeLists are not designed for split-apply-combine.
    ## trying is incredibly slow and truly painful... convert into a vanilla R data frame.
    ref.protein.df <- data.frame(ref = ref.proteins) %>%
        select(ref.start, ref.end, ref.width, ref.Note, ref.Name, ref.Alias, ref.ID) %>%
        dplyr::rename(Gene = ref.Name) %>%
        dplyr::rename(Alias = ref.Alias) %>%
        dplyr::rename(ID = ref.ID) %>%
        dplyr::rename(start = ref.start) %>%
        dplyr::rename(end = ref.end) %>%
        dplyr::rename(gene_length = ref.width) %>%
        dplyr::rename(product = ref.Note) %>%
        ## this next line filters out the tetA and CmR genes on the transposon and plasmid.
        dplyr::filter(!is.na(gene_length)) %>%
        ## IMPORTANT: some genes have multiple protein products.
        ## However, each protein corresponds to unique row in ref.protein.df.
        ## To address this issue, add a protein.index column, and
        ## split by protein.index rather than Gene (which is non-unique).
        mutate(protein.index = row_number())

    ## now get the coverage for each annotated protein in the reference genome.
    protein.coverage.df <- ref.protein.df %>%
        split(.$protein.index) %>%
        map_dfr(.f = get.coverage.per.protein) %>%
        ## this is a critical step: we need to label with the sample,
        mutate(Sample = gnome) %>%
        ## and add the sample metadata.
        left_join(row.df)

    return(protein.coverage.df)
}


evolved.protein.coverage.input.df <- mixedpop.input.df %>%
    left_join(metagenome.metadata) %>%
    select(-SampleType) %>%
    left_join(ancestral.clones.df)

## get coverage per gene for all evolved samples.
evolved.sample.protein.coverage.df <- evolved.protein.coverage.input.df %>%
    split(.$Sample) %>%
    map_dfr(.f = get.protein.coverage.per.sample)

## Make data structure for Figure 2D.
Fig2D.df <- evolved.sample.protein.coverage.df %>%
    mutate(is.B59.Tet5 = ifelse(Transposon == "B59", ifelse(Tet == 5, TRUE, FALSE), FALSE)) %>%
    relocate(Sample, Population, Transposon, Tet, Plasmid) %>%
    ## update the names of the Transposon factor for a prettier plot.
    mutate(Transposon_factor = fct_recode(as.factor(Transposon),
                                          `Tn5+` = "B30",
                                          `Tn5-` = "B59")) %>%
    mutate(Tet_factor = fct_recode(as.factor(Tet),
                                   `Tet 0` = "0",
                                   `Tet 5` = "5")) %>%
    ## update the names of the Plasmid factor for a prettier plot.
    mutate(Plasmid_factor = fct_recode(as.factor(Plasmid),
                                       `No plasmid` = "None",
                                       p15A = "p15A")) %>%
    ## unite the Transposon_factor, Tet_factor columns together.
    unite("Treatment", Transposon_factor:Tet_factor, sep="\n", remove = FALSE) %>%
    ## the product and Alias columns are lists!!
    ## turn these columns into character vectors,
    ## and handle the cases for tetA, KanR, CmR, and Tn5 genes,
    ## so that we can save this dataframe as Source Data.
    mutate(product = as.character(product)) %>%
    mutate(product = ifelse(Gene %in% c("tetA", "CmR", "KanR", "Tn5"), Gene, product)) %>%
    mutate(Alias = as.character(Alias)) %>%
    mutate(Alias = ifelse(Gene %in% c("tetA", "CmR", "KanR", "Tn5"), Gene, Alias)) %>%
    mutate(ID = ifelse(Gene %in% c("tetA", "CmR", "KanR", "Tn5"), Gene, ID))

## Write Figure 2D Source Data to file.
## IMPORTANT: one column includes a newline character, so we need quote=TRUE.
write.csv(Fig2D.df, "../results/Source-Data/Fig2D-Source-Data.csv", row.names = FALSE, quote=TRUE)
 
## make Figure 2D.
Fig2D <- ggplot(Fig2D.df, aes(x=start, y=mean_coverage, color = is.B59.Tet5)) +
    facet_grid(Treatment~Plasmid_factor) +
    geom_point(size=0.01,alpha=0.05) +
    theme_classic() +
    ## draw vertical lines at acrBAR starts.
    geom_vline(size=0.2,
               color = 'red',
               linetype = 'dashed',
               xintercept = c(480478, 484843, 484985)) +
    scale_x_continuous(name = "Genomic position (Mb)", breaks = c(0,1000000,2000000,3000000,4000000),
                       labels = c(0,1,2,3,4)) +
    scale_y_continuous(name = "Illumina read coverage", breaks = c(0,250,500), limits = c(0,575)) +
    theme(axis.text.y = element_text(size = 5)) +
    guides(color = "none") +
    scale_color_manual(values=c("gray", "blue"))
## 2 outlier points with very high coverage are removed from this plot
## by setting max coverage to 575x in this figure.
ggsave("../results/Fig2D.pdf", Fig2D, height=3,width=7)

## ALSO: let's look at genes that have no unique coverage-- these are repeats!
repeated.genes.by.coverage <- Fig2D.df %>%
    filter(mean_coverage < 50) %>%
    filter(Tet == 5) %>%
    filter(Transposon == "B59") %>%
    filter(start > 200000) %>%
    filter(end < 1000000) %>%
    select(Sample, Population, Treatment, start, end, gene_length, product, Gene, Alias, ID, mean_coverage) %>%
    arrange(start)
 
