##copy-number-analysis.R by Rohan Maddamsetti.

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
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)
library(DT)          # prettier data.frame output
library(data.table)  # faster fread()
library(dtplyr)      # dplyr works with data.table now.
library(cowplot)     # layout figures nicely.
library(xml2)
library(assertthat)

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
    transposon.avg <- as.numeric(xml_text(table.data[13]))
    transposon.relative.variance <- as.numeric(xml_text(table.data[14]))
    ## all samples should have these data.
    coverage.df <- data.frame('Sample' = basename(breseq.output.dir),
                              'mean'=c(chromosome.avg, transposon.avg),
                              'relative.variance'=c(chromosome.relative.variance, transposon.relative.variance),
                              'variance'=c(chromosome.avg * chromosome.relative.variance,
                                           transposon.avg * transposon.relative.variance),
                              'replicon'=c("chromosome", "transposon"))
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
find.K12.chromosomal.amplifications <- function(breseq.output.dir, gnome) { #gnome is not a misspelling.
    
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
        mutate(len=right.boundary-left.boundary) %>% filter(len>(2*max.read.len)) %>% mutate(amplication.index=row_number())

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
               copy.number.mean=coverage.mean/nbinom.fit$mean)
    
    ##print(data.frame(amplified.segments))
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
        mutate(Sample=as.character(gnome)) %>%
        mutate(bonferroni.corrected.pval=pval*alpha/bonferroni.alpha)
    
    return(significant.amplifications)
}


annotate.sample.amplifications <- function(sample.amplifications, ancestor.gff) {

    ancestor.gff <- unique(sample.amplifications$gff_path)
    
    ## create the IRanges object.
    amp.ranges <- IRanges(sample.amplifications$left.boundary,
                          sample.amplifications$right.boundary)
    ## Turn into a GRanges object in order to find overlaps with NEB5-alpha genes.
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


plot.amp.segments <- function(annotated.amps,clone.labels) {
    
    ## for annotated.amps and clone.labels to play nicely with each other.
    clone.labels$Name <- as.character(clone.labels$Name)
    
    labeled.annotated.amps <- left_join(annotated.amps,clone.labels,by=c("Genome" = 'Name')) %>%
        select(-query.index,-subject.index,-is.significant,-SampleType, -Population) %>%
        mutate(log.pval=log(bonferroni.corrected.pval)) %>%
        mutate(log2.copy.number.mean=log2(copy.number.mean)) %>%
        mutate(left.boundary.MB = left.boundary/1000000) %>%
        mutate(right.boundary.MB = right.boundary/1000000) ##%>%
##        mutate(Genome.Class=recode(Environment,
##                                   DM0 = "DM0-evolved genomes",
##                                   DM25 = "DM25-evolved genomes"))
    
    ## order the genes by start to get axes correct on heatmap.
    labeled.annotated.amps$gene <- with(labeled.annotated.amps, reorder(gene, start))
    ## reverse the order of genomes to make axes consistent with stacked barplot.
    labeled.annotated.amps$Genome <- factor(labeled.annotated.amps$Genome)
    labeled.annotated.amps$Genome <- factor(labeled.annotated.amps$Genome,
                                            levels=rev(levels(labeled.annotated.amps$Genome)))
    
    segmentplot <- ggplot(
        labeled.annotated.amps,
        aes(x=left.boundary.MB,
            xend=right.boundary.MB,
            y=Genome,
            yend=Genome,
            color=log2.copy.number.mean,
            size=20,
            frame=Genome.Class)) +
        geom_segment() +
        xlab("Genomic position (Mb)") +
        ylab("") +
        scale_color_viridis(name=bquote(log[2]~"(copy number)"),option="plasma") +
        facet_wrap(~Transposon,nrow=2, scales = "free_y") +
        theme_classic(base_family='Helvetica') +
        guides(size=FALSE) +
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
    inner_join(metagenome.metadata)

ancestral.clones.df <- read.csv("../data/one-day-expt-ancestral-sample-metadata.csv") %>%
   dplyr::rename(Ancestor = Sample) %>%
    select(-SampleType) %>%
    mutate(gff_name = paste0(Ancestor, ".gff3")) %>%
    mutate(gff_path = file.path(projdir, "results", "one-day-expt-genome-analysis", gff_name))


## Find chromosomal amplifications in all samples, and annotate with their ancestor gff file.
amps.with.ancestors <- map2_df(mixedpop.input.df$path,
                mixedpop.input.df$Sample,
                find.K12.chromosomal.amplifications) %>%
    ungroup() %>%
    left_join(metagenome.metadata) %>%
    select(-SampleType) %>%
    left_join(ancestral.clones.df)

annotated.amps <- annotate.amplifications(amps.with.ancestors)

## CRITICAL TODO:
## make figures, and polish this analysis further.

parallel.amplified.genes <- annotated.amps %>%
    group_by(gene, locus_tag, start, end) %>%
    summarize(parallel.amplifications = n()) %>%
    arrange(desc(parallel.amplifications)) %>%
    filter(parallel.amplifications > 2)

parallel.amplified.genes2 <- parallel.amplified.genes %>%
    arrange(start)

yaiT.amps <- annotated.amps %>%
    filter(gene == "yaiT")

acrABR.amps <- annotated.amps %>%
    filter(str_detect(gene, "acr"))

test <- amps.with.ancestors %>%
    filter(len>6000)

######################################################################
## Plot the plasmid/chromosome and transposon/chromosome ratio in each sample.
## CRITICAL TODO: JUST get the coverage for the actual synthetic transposon,
## not the entire whole region integrated into the chromosome.
replicon.coverage.df <- map_dfr(.x = mixedpop.input.df$path, .f = coverage.nbinom.from.html) %>%
    inner_join(metagenome.metadata) %>%
    ## I am not examining dispersion or variance at this point.
    select(Sample, mean, replicon, Transposon, Plasmid, Population, Tet)

replicon.coverage.ratio.df <- replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Transposon, Plasmid, Population, Tet) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome),
              plasmids.per.chromosome = (mean_plasmid/mean_chromosome),
              transposons.per.plasmid = (mean_transposon/mean_plasmid)) %>%
    pivot_longer(cols = c(transposons.per.chromosome,plasmids.per.chromosome,transposons.per.plasmid),
                 names_to = "ratio_type", values_to = "ratio")

Tet5.ratio.plot <- replicon.coverage.ratio.df %>%
    filter(Tet == 5) %>%
    ggplot(aes(y = ratio, x = Plasmid, color = ratio_type, shape = Transposon)) +
    geom_point() +
    theme_classic() +
    facet_wrap(ratio_type~Transposon, scales = "free") +
    theme(strip.background = element_blank()) +
    ggtitle("5 ug/mL tetracycline, Day 1") +
    guides(color= "none", shape = "none")

Tet0.ratio.plot <- replicon.coverage.ratio.df %>%
    filter(Tet == 0) %>%
    ggplot(aes(y = ratio, x = Plasmid, color = ratio_type, shape = Transposon)) +
    geom_point() +
    theme_classic() +
    facet_wrap(ratio_type~Transposon, scales = "free") +
    theme(strip.background = element_blank()) +
    ggtitle("0 ug/mL tetracycline, Day 1") +
    guides(color = "none", shape = "none")

ratio.figure.Fig4 <- plot_grid(Tet5.ratio.plot, Tet0.ratio.plot, labels=c('A','B'),nrow=2)
ggsave("../results/one-day-expt-coverage-ratios.pdf", ratio.figure.Fig4)

## let's write out the table too.
write.csv(replicon.coverage.ratio.df, "../results/one-day-expt-plasmid-transposon-coverage-ratios.csv",
          quote=F, row.names=FALSE)
