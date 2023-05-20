##antibiotic-generality-expt-analysis.R by Rohan Maddamsetti.

## S1 Figure A: plot of the change in transposon copy number.
## S1 Figure B: matrix plot 

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

## additional dependencies for matrix plot.
library(patchwork)
library(ggthemes)


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


annotate.amplifications <- function(evolved.amps) {
    evolved.amps %>% split(.$Sample) %>%
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
        xlab("Genomic position (Mb)") +
        ylab("") +
        ## draw a line at attB site for HK022 integration (between torS and torT)
        geom_vline(xintercept=1.056188,color="red") +
        scale_color_viridis(name="copy number",option="plasma") +
        facet_wrap(~Plasmid,nrow=2, scales = "free_y") +
        theme_classic(base_family='Helvetica') +
        guides(size= "none") +
        theme(legend.position="bottom") +
        theme(axis.ticks=element_line(size=0.1))
    return(segmentplot)
}


##############################################################################
## Functions for making the matrix figure.

MakeMutCountMatrix <- function(evolved.muts, show.all=FALSE) {
    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
        ## unite the Antibiotic and AntibioticConcentration columns together.
        unite("Treatment", Antibiotic:AntibioticConcentration, sep="\n", remove = FALSE) %>%
        group_by(Gene, Sample, Transposon, Plasmid, Antibiotic, AntibioticConcentration, Treatment) %>%
        summarize(mutation.count = n()) %>%
        ## This is for sorting mutations.
        mutate(is.MOB = ifelse(str_detect(Gene,"mini"), TRUE, FALSE))
    
    total.muts <- matrix.data %>%
        group_by(Gene) %>%
        summarize(total.mutation.count = sum(mutation.count))
    
    matrix.data <- left_join(matrix.data, total.muts)
    
    if (!show.all) { ## then filter out genes that are only hit in one sample.
        matrix.data <- matrix.data %>%
            filter(total.mutation.count > 1)
    }
    
    ## sort genes by number of mutations in each row, but put all the transposon mutations together.
    ## also check out the alternate sorting method that follows.
    gene.hit.sort <- matrix.data %>%
        group_by(Gene, is.MOB, .drop = FALSE) %>%
        summarize(hits=sum(mutation.count)) %>%
        arrange(desc(is.MOB), desc(hits))
    ## now sort genes.
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))
    ## cast mutation.count into a factor for plotting.
    matrix.data$mutation.count <- factor(matrix.data$mutation.count)
    
    return(matrix.data)
}


MakeMatrixPanel <- function(mdata, treatment, leg=FALSE) {
    panel.data <- filter(mdata,Treatment==treatment)
    fig <- ggplot(panel.data,
                  aes(x=Sample,
                      y=Gene,
                      fill=mutation.count,
                      frame=Treatment)
                  ) +
        geom_tile(color="black",size=0.1) +
        ggtitle(treatment) +
        theme_tufte(base_family='Helvetica') +
        theme(axis.ticks = element_blank(),
              axis.text.x = element_text(size=10,angle=45,hjust=1),
              axis.text.y = element_text(size=10,hjust=1,face="italic"),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
        scale_fill_manual(name="Mutations",
                          values = c("#ffdf00", "#bebada", "#fb8072", "#80b1d3", "#fdb462"))
    
    if (leg == FALSE) {
        fig <- fig + guides(fill = "none")
    }
    return(fig)
}


##############################################################################
## K-12 MG1655 oriC replication origin annotation
## annotated as rep_origin in the genbank file.
## Also see: https://biocyc.org/ECOLI/NEW-IMAGE?type=EXTRAGENIC-SITE&object=G0-10506
## 3,925,744 -> 3,925,975
K12_oriC_START = 3925744
K12_oriC_END = 3925975
K12_oriC_MID = (K12_oriC_START+K12_oriC_END)/2


rotate.K12.chr <- function(my.position) {
    #' function to rotate genome coordinates,
    #' setting oriC at the center of plots
    #' that examine mutation bias over the chromosome.
    ## we want to change coordinates so that c is the new origin.
    GENOME.LENGTH <- 4641652
    midpoint <- GENOME.LENGTH/2
    oriC <- 3925860
    
    if (oriC >= midpoint) {
        L <- oriC - midpoint
        ifelse(my.position > L, my.position - oriC, GENOME.LENGTH - oriC + my.position)
    } else { ## midpoint is greater than new.origin.
        L <- midpoint + oriC
        ifelse(my.position > L, my.position - GENOME.LENGTH - oriC, my.position - oriC)
    }
}

#######################################
## Analysis time!

## assert that we are in the src directory, such that
## proj.dir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("ARG-duplications","src")))
projdir <- file.path("..")

## get metadata for all the evolved population metagenomes.
metagenome.metadata <- read.csv("../data/generality-expts-sample-metadata.csv")

breseq.output.dir <- file.path(projdir, "results/generality-expts-genome-analysis")
all.samples <- list.files(breseq.output.dir,pattern='^RM')
all.sample.paths <- sapply(all.samples, function(x) file.path(breseq.output.dir,x))

all.evolved.pop.input.df <- data.frame(Sample=all.samples, path=all.sample.paths) %>%
    ## IMPORTANT: only include samples in metagenome.metadata.
    right_join(metagenome.metadata)

## Analyze the antibiotic generality expt data.
evolved.pop.input.df <- filter(all.evolved.pop.input.df, Antibiotic != "Tetracycline")

######################################################################
## Plot the plasmid/chromosome and transposon/chromosome ratio in each sample.

## Get the actual coverage for the tetA gene in each transposon.
## This is calculated by get-generality-expts-transposon-coverage.py.
transposon.coverage.file <- file.path(projdir,
                                      "results/generality-expts-genome-analysis",
                                      "generality-expts-transposon-coverage.csv")
transposon.coverage.df <- read.csv(transposon.coverage.file) %>%
    ## let's add metadata and rename columns for compatibility.
    dplyr::rename(mean = TransposonCoverage) %>%
    dplyr::mutate(replicon = "transposon")

## data for the evolved samples.
all.evolved.transposon.coverage.df <- metagenome.metadata %>%
    left_join(transposon.coverage.df) %>%
    ## don't need this when joining to evolved.replicon.coverage.df
    select(-SampleType)

## Analyze the antibiotic generality expt data.
evolved.transposon.coverage.df <- filter(all.evolved.transposon.coverage.df, Antibiotic != "Tetracycline")

evolved.replicon.coverage.df <- map_dfr(.x = evolved.pop.input.df$path, .f = coverage.nbinom.from.html) %>%
    inner_join(metagenome.metadata) %>%
    ## I am not examining dispersion or variance at this point.
    select(Sample, mean, replicon, Transposon, Plasmid, Population, Antibiotic, AntibioticConcentration) %>%
    ## add rows for the transposon coverage data.
    full_join(evolved.transposon.coverage.df) %>%
    ## set NA coverage values to zero.
    mutate(mean=sapply(mean, function(x) ifelse(is.na(x), 0,x)))


evolved.replicon.coverage.ratio.df <- evolved.replicon.coverage.df %>%
    pivot_wider(names_from = replicon, values_from = mean, names_prefix = "mean_") %>%
    group_by(Sample, Transposon, Plasmid, Population, Antibiotic, AntibioticConcentration) %>%
    summarise(transposons.per.chromosome = (mean_transposon/mean_chromosome))
              

generality.ratio.fig.df <- evolved.replicon.coverage.ratio.df %>%
    mutate(`Transposon copy number` = transposons.per.chromosome) %>%
    mutate(Population = as.factor(Population)) %>%
    mutate(Day = ifelse(AntibioticConcentration==0, 0, 1)) %>%
    mutate(Day= factor(Day))

           
S1FigA <-  ggplot(generality.ratio.fig.df,
                  aes(x = Day,
                      y = `Transposon copy number`,
                      color = Antibiotic)) +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    theme(strip.background = element_blank()) +
    ylab("transposons per chromosome") +
    ylim(0.5,2.3) ## to get all points into the plot.


## let's write out the table too.
write.csv(evolved.replicon.coverage.ratio.df, "../results/antibiotic-generality-expt-transposon-coverage-ratios.csv",
          quote=F, row.names=FALSE)

################################################################################
### Figure: make a matrix plot of genes with mutations in two or more clones.
################################################################################

## colorblind-friendly palette.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## This file *just* has the evolved populations (Day 9).
all.pop.clone.labels <- read.csv(
  file.path(projdir,
            "data/generality-expts-sample-metadata.csv"),
  stringsAsFactors=FALSE)

## This is the key data file for the analysis.
all.evolved.mutations <- read.csv(
    file.path(projdir,
              "results/generality-expts-genome-analysis/evolved_mutations.csv"),
    stringsAsFactors=FALSE) %>%
    mutate(Mbp.coordinate=Position/1000000) %>%
    ## update the names of the Plasmid factor for a prettier plot.
    mutate(Plasmid_factor = fct_recode(as.factor(Plasmid),
                                `No plasmid` = "None",
                                p15A = "p15A"))

## let's focus the analysis on the experiments generalizing across antibiotics.
pop.clone.labels <- filter(all.pop.clone.labels, Antibiotic != "Tetracycline")
evolved.mutations <- filter(all.evolved.mutations, Antibiotic != "Tetracycline") %>%
    ## these parallel mutations looks like false positives to me... filter out.
    filter(Gene != "gltP/yjcO")

## make the matrix for the antibiotic generality expt.
non.tet.evolved.mut.count.matrix <- MakeMutCountMatrix(evolved.mutations, show.all=TRUE)

## look at the levels, because each level is one panel.
unique(non.tet.evolved.mut.count.matrix$Treatment)

B90.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "Spectinomycin\n250")
## Remove the gene labels to save space.
B91.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "Kanamycin\n250") +
    theme(axis.text.y=element_blank())

## get the legend from this panel, because this shows all the colors.
B92.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "Carbenicillin\n2000", leg=TRUE) +
    theme(axis.text.y=element_blank())

B95.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "Chloramphenicol\n70") +
    theme(axis.text.y=element_blank())

## get the legend.
Fig.legend <- get_legend(B92.matrix.panel)
## now remove the legend from the panel.
B92.matrix.panel <- B92.matrix.panel + guides(fill = "none")

matrix.panels <-
    B90.matrix.panel +
    B91.matrix.panel +
    B92.matrix.panel +
    B95.matrix.panel +
    Fig.legend +
    plot_layout(nrow = 1)

## hack to label x-axis from comments at: https://github.com/thomasp85/patchwork/issues/150
matrix.panels.grob <- patchwork::patchworkGrob(matrix.panels)
S1FigB <- gridExtra::grid.arrange(matrix.panels.grob, left = "", bottom = "Evolved populations")

#################################################################################
## Now make S1 Figure.
## save the figure.
S1Fig <- cowplot::plot_grid(S1FigA, S1FigB, labels=c("A","B"),nrow=2)
ggsave(filename="../results/S1Fig.pdf", plot=S1Fig, height=9, width=8)

