## one-day-expt-metagenomics.R by Rohan Maddamsetti.

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
pop.clone.labels <- read.csv(
  file.path(projdir,
            "data/one-day-expt-evolved-sample-metadata.csv"),
  stringsAsFactors=FALSE) %>%
    ## remove the pUC samples from the analysis
    filter(Plasmid != "pUC")

## This is the key data file for the analysis.
evolved.mutations <- read.csv(
    file.path(projdir,
              "results/one-day-expt-genome-analysis/evolved_mutations.csv"),
    stringsAsFactors=FALSE) %>%
    mutate(Mbp.coordinate=Position/1000000) %>%
    ## remove the pUC samples from the analysis.
    filter(Plasmid != "pUC") %>%
    ## update the names of the Transposon factor for a prettier plot.
    mutate(Transposon_factor = fct_recode(as.factor(Transposon),
                                   `Tn5+` = "B30",
                                   `Tn5-` = "B59")) %>%
    ## update the names of the Plasmid factor for a prettier plot.
    mutate(Plasmid_factor = fct_recode(as.factor(Plasmid),
                                `No plasmid` = "None",
                                p15A = "p15A")) %>%
    mutate(Tet_factor = fct_recode(as.factor(Tet),
                                   `Tet 0` = "0",
                                   `Tet 5` = "5"))



evolved.MOB <- filter(evolved.mutations, Mutation=='MOB')

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
        facet_grid(Plasmid~.) +
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


Tet0.evolved.mutations <- filter(evolved.mutations, Tet==0)
Tet5.evolved.mutations <- filter(evolved.mutations, Tet==5)

Tet0.freq.panel <- make.allele.freq.histogram(Tet0.evolved.mutations, "Tet 0 populations",TRUE)
Tet5.freq.panel <- make.allele.freq.histogram(Tet5.evolved.mutations, "Tet 5 populations",TRUE)

## This is a very important plot: what does this distribution say about
## the possibility of false positives? how can I interpret this?
## any theoretical basis in population genetics?

## Idea for an empirical control.
## 1) downsample reads from the treatment without plasmid, and re-run breseq
## to see if false positive mutation calls arise when coverage is ~40X rather than
## 300X.

## IMPORTANT TODO: There seems to be a bug in which "missing data" is removed, but right now I have
## no idea what this is about or what is being removed from the plot. Figure this out!!!
freq.panel <- plot_grid(Tet0.freq.panel, Tet5.freq.panel, labels = c('A','B'))
freq.panel.output <- "../results/one-day-expt-allele-freqs.pdf"
ggsave(freq.panel, file=freq.panel.output,width=10,height=4)

###############################################
## Figure 2B: make a stacked bar plot of the kinds of mutations in each treatment.

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
        group_by(Sample, Transposon, Plasmid, Tet, Population, Mutation, Transposon_factor, Plasmid_factor, Tet_factor) %>%
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
        ## show both tetracycline concentrations.
        facet_wrap(Transposon_factor~Tet_factor) +
        geom_bar(stat='identity') +
        scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +        
        theme_classic(base_family='Helvetica') +
        theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
              axis.text.y=element_text(size=12),
              legend.text = element_text(size = 14),
              strip.text = element_text(size = 12),
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


## Write Figure 2B Source Data.
mutation.class.df <- make.mutation.class.df(evolved.mutations)
write.csv(mutation.class.df, "../results/Source-Data/Fig2B-Source-Data.csv", row.names=FALSE, quote=FALSE)
## Now make Figure 2B.
Fig2B <- plot.mutation.summary.stackbar(mutation.class.df, TRUE, FALSE)
fig2B.output <- "../results/Fig2B.pdf"
ggsave(Fig2B, file=fig2B.output,width=6,height=5)


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
    select(Gene, count, Plasmid, Transposon, Tet) %>%
    distinct() %>%
    arrange(desc(count))

#################################################################################
### Figure 2C: make a matrix plot of genes with mutations in two or more clones.
################################################################################
MakeMutCountMatrixFigure <- function(evolved.muts, show.all=FALSE, use.treatment.hit.sort=FALSE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
        ## unite the Transposon_factor, Plasmid_factor, Tet_factor columns together.
        unite("Treatment", Transposon_factor:Tet_factor, sep="\n", remove = FALSE) %>%
        group_by(Gene, Sample, Transposon, Plasmid, Tet, Treatment) %>%
        summarize(mutation.count = n()) %>%
        ## This is for sorting mutations.
        mutate(is.MOB = ifelse(str_detect(Gene,"tetA-Tn5"), TRUE, FALSE))
    
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
    if (use.treatment.hit.sort) {
        ## alternate sorting method: difference in hits between environments,
        ## AKA the (absolute value of the) difference in number of pops with hits
        ## between the None and p15A treatments.
        p15A.hit.count.df <- filter(matrix.data,Plasmid=="p15A") %>%
            group_by(Gene, is.MOB, .drop = FALSE) %>%
            summarize(p15A.hit.count=n())
        
        noPlasmid.hit.count.df <- filter(matrix.data,Plasmid=="None") %>%
            group_by(Gene, is.MOB, .drop = FALSE) %>%
            summarize(noPlasmid.hit.count=n())
        
        treatment.hit.sort <- full_join(p15A.hit.count.df, noPlasmid.hit.count.df) %>%
            mutate(hit.diff = p15A.hit.count - noPlasmid.hit.count) %>%
            arrange(desc(is.MOB), desc(hit.diff))

        ## cast Gene into a factor for plotting.
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(treatment.hit.sort$Gene))
    } else {
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.hit.sort$Gene))
    }
    
    ## cast mutation.count into a factor for plotting.
    matrix.data$mutation.count <- factor(matrix.data$mutation.count)

    make.matrix.panel <- function(mdata, treatment, leg=FALSE) {
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

    
    ## make Tet5 panels.
    ## get the legend from the first panel.
    Inactive.noPlasmid.Tet5.matrix.panel <- make.matrix.panel(matrix.data, "Tn5-\nNo plasmid\nTet 5", leg=TRUE)

    Fig.legend <- get_legend(Inactive.noPlasmid.Tet5.matrix.panel)

    ## now remove the legend from the first panel.
    Inactive.noPlasmid.Tet5.matrix.panel <- Inactive.noPlasmid.Tet5.matrix.panel + guides(fill = "none")
        
    ## Remove the gene labels to save space.
    Inactive.A31.Tet5.matrix.panel <- make.matrix.panel(matrix.data, "Tn5-\np15A\nTet 5") +
        theme(axis.text.y=element_blank())
    
    Active.noPlasmid.Tet5.matrix.panel <- make.matrix.panel(matrix.data,"Tn5+\nNo plasmid\nTet 5") +
        theme(axis.text.y=element_blank())
    Active.A31.Tet5.matrix.panel <- make.matrix.panel(matrix.data, "Tn5+\np15A\nTet 5") +
        theme(axis.text.y=element_blank())
    
    ## Using the patchwork library for layout.
    matrix.panels <-
        Inactive.noPlasmid.Tet5.matrix.panel +
        Inactive.A31.Tet5.matrix.panel +
        Active.noPlasmid.Tet5.matrix.panel +
        Active.A31.Tet5.matrix.panel +
        Fig.legend +
        plot_layout(nrow = 1)

    ## hack to label x-axis from comments at: https://github.com/thomasp85/patchwork/issues/150
    matrix.panels.grob <- patchwork::patchworkGrob(matrix.panels)
    matrix.figure <- gridExtra::grid.arrange(matrix.panels.grob, left = "", bottom = "Evolved populations")
    
    return(matrix.figure)
}


## Use summed allele frequency for the heatmap.
MakeSummedAlleleFrequencyMatrixFigure <- function(evolved.muts,
                                                  allele.freq.threshold = 0.20, ## this parameter only matters if show.all == FALSE.
                                                  show.all=FALSE, ## if TRUE, all mutations are shown, regardless of allele frequency.
                                                  use.treatment.hit.sort=FALSE,
                                                  add.legend = TRUE) {

    ## First, make a mutation matrix for plotting.
    matrix.data <- evolved.muts %>%
        ## unite the Transposon_factor, Plasmid_factor, Tet_factor columns together.
        unite("Treatment", Transposon_factor:Tet_factor, sep="\n", remove = FALSE) %>%
        group_by(Gene, Sample, Transposon, Plasmid, Tet, Treatment) %>%
        summarize(summed.Allele.Frequency = sum(Frequency)) %>%
        ## This is for sorting mutations.
        mutate(is.MOB = ifelse(str_detect(Gene,"tetA-Tn5"), TRUE, FALSE))

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
    ## also check out the alternate sorting method that follows.
    gene.freq.sort <- matrix.data %>%
        group_by(Gene, is.MOB, .drop = FALSE) %>%
        summarize(totalallelefreq = sum(summed.Allele.Frequency)) %>%
        arrange(desc(is.MOB), desc(totalallelefreq))
    
    ## alternative sorting method:
    ## difference in allele frequency between the p15A and noPlasmid treatments..
    p15A.allele.freq.df <- filter(matrix.data, Plasmid=="p15A") %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(p15A.allele.frequency=sum(summed.Allele.Frequency))
    
    noPlasmid.allele.freq.df <- filter(matrix.data, Plasmid=="None") %>%
        group_by(Gene, .drop = FALSE) %>%
        summarize(noPlasmid.allele.frequency = sum(summed.Allele.Frequency))
    
    treatment.freq.sort <- full_join(p15A.allele.freq.df, noPlasmid.allele.freq.df) %>%
        replace_na(list(p15A.allele.frequency = 0, p15A.allele.frequency = 0)) %>%
        mutate(allele.diff = noPlasmid.allele.frequency - noPlasmid.allele.frequency) %>%
        arrange(desc(allele.diff))

    ## sort the genes.
    if (use.treatment.hit.sort) {
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(treatment.freq.sort$Gene))
    } else {
        matrix.data$Gene <- factor(matrix.data$Gene,levels=rev(gene.freq.sort$Gene))
    }

    make.allele.freq.matrix.panel <- function(mdata, treatment, leg=use.legend) {
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
        

        if (leg == FALSE) {
            fig <- fig + guides(fill= "none")
        }
        return(fig)
    }
    

    Inactive.noPlasmid.Tet5.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "Tn5-\nNo plasmid\nTet 5", add.legend)
    ## Remove the gene labels for the additional matrices to save space.
    Inactive.A31.Tet5.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "Tn5-\np15A\nTet 5", add.legend)  +
        theme(axis.text.y=element_blank())
    
    Active.noPlasmid.Tet5.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "Tn5+\nNo plasmid\nTet 5", add.legend)  +
        theme(axis.text.y=element_blank())
    Active.A31.Tet5.matrix.panel <- make.allele.freq.matrix.panel(matrix.data, "Tn5+\np15A\nTet 5", add.legend) +
        theme(axis.text.y=element_blank())

    if (add.legend) {
        ## get the legend from the last panel.
        my.legend <- get_legend(Active.A31.Tet5.matrix.panel)
        ## now remove the legend from the panel.
        Active.A31.Tet5.matrix.panel <- Active.A31.Tet5.matrix.panel + guides(fill = "none")
    }
    
    ## Using the patchwork library for layout.
    matrix.figure <-
        Inactive.noPlasmid.Tet5.matrix.panel +
        Inactive.A31.Tet5.matrix.panel +
        Active.noPlasmid.Tet5.matrix.panel +
        Active.A31.Tet5.matrix.panel +
        plot_layout(nrow = 1) +
        plot_layout(guides = "collect") & theme(legend.position = 'bottom')

    return(matrix.figure)
}

## Write Figure 2C Source Data.
write.csv(Tet5.evolved.mutations, "../results/Source-Data/Fig2C-Source-Data.csv", row.names=FALSE, quote=FALSE)
## Figure 2C.
Fig2C <- MakeMutCountMatrixFigure(Tet5.evolved.mutations, show.all=TRUE, use.treatment.hit.sort=TRUE)
Fig2C.outf <- "../results/Fig2C.pdf"
ggsave(Fig2C.outf, Fig2C, height=4, width=7)

 
