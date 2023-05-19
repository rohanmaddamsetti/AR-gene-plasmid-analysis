## generality-expts-metagenomics.R by Rohan Maddamsetti.

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggthemes)
library(viridis)
library(ggrepel)

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


## colorblind-friendly palette.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## assert that we are in the src directory, such that
## projdir is the parent of the current directory.
stopifnot(endsWith(getwd(), file.path("ARG-duplications","src")))
projdir <- file.path("..")

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
evolved.mutations <- filter(all.evolved.mutations, Antibiotic != "Tetracycline")

###############################################
## Plot the distribution of measured allele frequencies in each population.
make.allele.freq.histogram <- function(evolved.mutations.df, my.title,annotate=FALSE) {
    p <- ggplot(evolved.mutations.df, aes(x=Frequency)) +
        geom_histogram(bins = 200) +
        theme_classic() +
        ylab("Count") +
        xlab("Allele Frequency") +
        scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1), limits = c(0,1.1)) +
        ggtitle(my.title) +
        facet_wrap(Plasmid~Antibiotic) +
    geom_vline(xintercept=0.10,color="red",linetype="dashed",size=0.2)

    muts.to.label <- filter(evolved.mutations.df, Frequency>0.5)
    if (annotate && nrow(muts.to.label) > 0) {
        p <- p +
            geom_text_repel(
                ## label mutations with > 25% allele frequency
                data= muts.to.label,
                aes(x=Frequency,y=1,label=Gene),
                fontface = "italic",size=1.5,show.legend=FALSE,inherit.aes=FALSE)
        }
    return(p)
}


no.antibiotic.evolved.mutations <- filter(evolved.mutations, AntibioticConcentration==0)
antibiotic.evolved.mutations <- filter(evolved.mutations, AntibioticConcentration>0)

no.antibiotic.freq.panel <- make.allele.freq.histogram(no.antibiotic.evolved.mutations, "No antibiotic populations",TRUE)
antibiotic.freq.panel <- make.allele.freq.histogram(antibiotic.evolved.mutations, "Antibiotic treated populations",TRUE)

## This is a very important plot: what does this distribution say about
## the possibility of false positives? how can I interpret this?
## any theoretical basis in population genetics?

## Idea for an empirical control.
## 1) downsample reads from the treatment without plasmid, and re-run breseq
## to see if false positive mutation calls arise when coverage is ~40X rather than
## 300X.

## IMPORTANT TODO: There seems to be a bug in which "missing data" is removed, but right now I have
## no idea what this is about or what is being removed from the plot. Figure this out!!!
freq.panel <- plot_grid(no.antibiotic.freq.panel, antibiotic.freq.panel, labels = c('A','B'))
freq.panel.output <- "../results/generality-expt-allele-freqs.pdf"
ggsave(freq.panel, file=freq.panel.output,width=10,height=4)

###############################################
## make a stacked bar plot of the kinds of mutations in each treatment.

## This function sums mutations per replicate population.
make.mutation.class.df <- function(evolved.mutations.df) {
    evolved.mutations.df %>%
        ## give nicer names for mutation classes.
        mutate(Mutation=recode(Mutation,
                               MOB = "Mobile element transposition",
                               DEL = "Indel",
                               INS = "Indel",
                               SUB = "Multiple-base substitution",
                               nonsynonymous = "Nonsynonymous",
                               synonymous = "Synonymous",
                               nonsense = "Nonsense",
                               pseudogene = "Pseudogene",
                               intergenic = "Intergenic",
                               )) %>%
        group_by(Sample, Transposon, Plasmid, Antibiotic, Population, Mutation, AntibioticConcentration) %>%
        summarize(Count=n(),WeightedCount = sum(Frequency)) %>%
        ungroup() %>%
        data.frame() %>%
        mutate(Mutation=as.factor(as.character(Mutation)))
}


plot.mutation.summary.stackbar <- function(mutation.class.df, leg=FALSE, weight.by.freq=FALSE) {

    if (weight.by.freq) {
        fig <- ggplot(mutation.class.df, aes(x=Plasmid, y=WeightedCount, fill=Mutation)) +
            ylab("Summed Allele Frequency")
    } else {
        fig <- ggplot(mutation.class.df, aes(x=Plasmid, y=Count, fill=Mutation)) +
            ylab("Count")
    }

    fig <- fig +
        facet_wrap(AntibioticConcentration~Antibiotic) +
        geom_bar(stat='identity') +
        scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +        
        theme_classic(base_family='Helvetica') +
        theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
              axis.text.y=element_text(size=12),
              panel.border=element_blank(),
              strip.background = element_blank(),
              panel.spacing.x=unit(1, "cm"),
              panel.spacing.y=unit(0.5, "cm"))

    if (leg == TRUE) {
        fig <- fig +
            theme(legend.title=element_text(size=8, face="bold"),
                  legend.title.align=0.5,
                  legend.text=element_text(size=8),
                  legend.position="bottom")
    } else {
        fig <- fig + guides(fill = "none")
    }
    
    return(fig)
}

## Now make a Figure.
mutation.class.df <- make.mutation.class.df(evolved.mutations)
generality.expts.mutation.classes.plot <- plot.mutation.summary.stackbar(mutation.class.df, TRUE, FALSE)
mutation.classes.output <- "../results/generality-expts-mutation-classes.pdf"
ggsave(generality.expts.mutation.classes.plot, file=mutation.classes.output,width=6,height=5)


#################################################################################
## analysis of parallel evolution at the same nucleotide.
## discuss numbers and finding in the text (no figure.).
## This could be a Supplementary Table.

bp.parallel.mutations <- evolved.mutations %>% group_by(Position) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.MOB <- filter(bp.parallel.mutations,Mutation=='MOB')
## no parallel DEL muts at bp level.
parallel.INS <- filter(bp.parallel.mutations,Mutation=='INS')
parallel.dN <- filter(bp.parallel.mutations,Mutation=='nonsynonymous')
parallel.dS <- filter(bp.parallel.mutations,Mutation=='synonymous')

## examine parallel evolution at amino acid level (only one case, in robA).
parallel.AA.dN <- evolved.mutations %>% filter(Mutation=='nonsynonymous') %>% group_by(Position) %>% summarize(count=n()) %>% filter(count > 1)
parallel.dN.Table <- filter(evolved.mutations, Position %in% parallel.AA.dN$Position) %>% arrange(Position)

## check parallel evolution for synonymous mutations too.
parallel.dS.Table <- filter(evolved.mutations, Position %in% parallel.dS$Position) %>% arrange(Position)


#################################################################################
## analysis of parallel evolution at the gene level (including intergenic regions).

gene.level.parallel.mutations <- evolved.mutations %>% group_by(Gene) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.genes <- gene.level.parallel.mutations %>%
    select(Gene, count, Plasmid, Transposon, Antibiotic, AntibioticConcentration) %>%
    distinct() %>%
    arrange(desc(count))

#################################################################################
### Figure: make a matrix plot of genes with mutations in two or more clones.
################################################################################
MakeMutCountMatrix <- function(evolved.muts, show.all=FALSE) {
    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
        ## unite the Transposon, Plasmid, Antibiotic, AntibioticConcentration columns together.
        unite("Treatment", Transposon:AntibioticConcentration, sep="\n", remove = FALSE) %>%
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


MakeSummedAlleleFrequencyMatrix <- function (evolved.muts,
                                             allele.freq.threshold = 0.20, ## this parameter only matters if show.all == FALSE.
                                             show.all=FALSE) { ## if TRUE, all mutations are shown, regardless of allele frequency.

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
        ## unite the Transposon, Plasmid, Antibiotic, AntibioticConcentration columns together.
        unite("Treatment", Transposon:AntibioticConcentration, sep="\n", remove = FALSE) %>%
        group_by(Gene, Sample, Transposon, Plasmid, Antibiotic, AntibioticConcentration, Treatment) %>%
        summarize(summed.Allele.Frequency = sum(Frequency)) %>%
        ## This is for sorting mutations.
        mutate(is.MOB = ifelse(str_detect(Gene,"mini"), TRUE, FALSE))

    total.allele.freqs <- matrix.data %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(total.Allele.Frequency = sum(summed.Allele.Frequency))

    matrix.data <- left_join(matrix.data, total.allele.freqs)
    
    ## filter matrix.data for genes that pass the allele frequency threshold,
    ## based on total allele frequency summed across all pops.
    if (!show.all) {
        matrix.data <- matrix.data %>%
            filter(total.Allele.Frequency > allele.freq.threshold)
    }

    ## sort genes by the total allele frequency in each row, but put all the transposon mutations together.
    gene.freq.sort <- matrix.data %>%
        group_by(Gene, is.MOB, .drop = FALSE) %>%
        summarize(totalallelefreq = sum(summed.Allele.Frequency)) %>%
        arrange(desc(is.MOB), desc(totalallelefreq))
    ## sort the genes.
    matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.freq.sort$Gene))

    return(matrix.data)
}


MakeAlleleFreqMatrixPanel <- function(mdata, treatment, use.legend=TRUE) {
    fig <- ggplot(filter(mdata,Treatment==treatment),
                  aes(x=Sample,
                      y=Gene,
                      fill=summed.Allele.Frequency,
                      frame=Treatment)
                  ) +
        geom_tile(color="black",size=0.1) +
        ggtitle(treatment) +
        theme_tufte(base_family='Helvetica') +
        theme(axis.ticks = element_blank(),
              axis.text.x = element_text(size=10,angle=45,hjust=1),
              axis.text.y = element_text(size=10,hjust=1,face="italic"),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              ) +
        scale_y_discrete(drop=FALSE) + ## don't drop missing genes.
        theme(legend.position = "bottom") + ## arrange the legend on the bottom.
        scale_fill_viridis_c(option = "inferno", limits = c(0,1)) ## important: need a uniform scale across samples.
    
    if (use.legend == FALSE) {
        fig <- fig + guides(fill= "none")
    }
    return(fig)
}


##########################
## Make mutation matrices.
## first make transposon generality matrix.
tet.evolved.mutations <- filter(antibiotic.evolved.mutations, Antibiotic == "Tetracycline")
tet.evolved.mut.count.matrix <- MakeMutCountMatrix(tet.evolved.mutations, show.all=FALSE)
## look at the levels, because each level is one panel.
unique(tet.evolved.mut.count.matrix$Treatment)

## make one matrix for the Tet5, multiple transposons experiment.
## make Tet5 panels.
B107.noPlasmid.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                                    "B107\nNone\nTetracycline\n5")
## Remove the gene labels to save space.
B111.noPlasmid.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                                    "B111\nNone\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
B123.noPlasmid.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                                    "B123\nNone\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
B109.noPlasmid.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                                    "B109\nNone\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
B110.noPlasmid.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                                    "B110\nNone\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
## p15A panels.
B107.p15A.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                               "B107\np15A\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
B111.p15A.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                               "B111\np15A\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
B123.p15A.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                               "B123\np15A\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
B109.p15A.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                               "B109\np15A\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
B110.p15A.Tet5.matrix.panel <- MakeMatrixPanel(tet.evolved.mut.count.matrix,
                                               "B110\np15A\nTetracycline\n5") +
    theme(axis.text.y=element_blank())
    
## Using the patchwork library for layout.
transposon.generality.matrix.figure <-
    B107.noPlasmid.Tet5.matrix.panel +
    B111.noPlasmid.Tet5.matrix.panel +
    B123.noPlasmid.Tet5.matrix.panel +
    B109.noPlasmid.Tet5.matrix.panel +
    B110.noPlasmid.Tet5.matrix.panel +
    ## p15A panels.
    B107.p15A.Tet5.matrix.panel +
    B111.p15A.Tet5.matrix.panel +
    B123.p15A.Tet5.matrix.panel +
    B109.p15A.Tet5.matrix.panel +
     B110.p15A.Tet5.matrix.panel +
    plot_layout(nrow = 1)
## save the figure.
ggsave("../results/transposon-generality-mutation-matrix.pdf",
       transposon.generality.matrix.figure, height=4, width=12)


##################################################
## now make the matrix for the antibiotic generality expt.
non.tet.evolved.mutations <- filter(antibiotic.evolved.mutations, Antibiotic != "Tetracycline")
non.tet.evolved.mut.count.matrix <- MakeMutCountMatrix(non.tet.evolved.mutations, show.all=TRUE)

## look at the levels, because each level is one panel.
unique(non.tet.evolved.mut.count.matrix$Treatment)

B90.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "B90\nNone\nSpectinomycin\n250")
## Remove the gene labels to save space.
B91.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "B91\nNone\nKanamycin\n250") +
    theme(axis.text.y=element_blank())
B92.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "B92\nNone\nCarbenicillin\n2000") +
    theme(axis.text.y=element_blank())
B95.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "B95\nNone\nChloramphenicol\n70") +
    theme(axis.text.y=element_blank())

antibiotic.generality.matrix.figure <-
    B90.matrix.panel +
    B91.matrix.panel +
    B92.matrix.panel +
    B95.matrix.panel +
    plot_layout(nrow = 1)

## save the figure.
ggsave("../results/antibiotic-generality-mutation-matrix.pdf",
       antibiotic.generality.matrix.figure, height=6, width=6)

