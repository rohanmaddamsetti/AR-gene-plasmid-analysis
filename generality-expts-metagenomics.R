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

#################################################################################
### Figure: make a matrix plot of genes with mutations in two or more clones.
################################################################################
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


##########################
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
B92.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "Carbenicillin\n2000") +
    theme(axis.text.y=element_blank())
B95.matrix.panel <- MakeMatrixPanel(non.tet.evolved.mut.count.matrix,
                                     "Chloramphenicol\n70") +
    theme(axis.text.y=element_blank())

antibiotic.generality.matrix.figure <-
    B90.matrix.panel +
    B91.matrix.panel +
    B92.matrix.panel +
    B95.matrix.panel +
    plot_layout(nrow = 1)

## save the figure.
ggsave("../results/antibiotic-generality-mutation-matrix.pdf",
       antibiotic.generality.matrix.figure, height=6, width=7)

