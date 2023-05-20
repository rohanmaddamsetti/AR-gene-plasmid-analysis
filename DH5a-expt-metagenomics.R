## DH5a-expt-metagenomics.R by Rohan Maddamsetti.

## Here's a list of NEB reference genomes. I use the most recent version of NEB5-alpha as reference.
## https://international.neb.com/tools-and-resources/usage-guidelines/competent-e-coli-genome-sequences-tool

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggthemes)
library(viridis)
library(ggrepel)

## DH5a origin, based on aligning Jeff Barrick's manual annotation of the
## REL606 oriC sequence against NZ_CP017100 using NCBI BLAST.
NZ_CP017100_oriC_START = 3866291
NZ_CP017100_oriC_END = 3866522
NZ_CP017100_oriC_MID = (NZ_CP017100_oriC_START+NZ_CP017100_oriC_END)/2

## GC skew calculations using the webserver at:
## https://genskew.csb.univie.ac.at/webskew
## Also, see: https://skewdb.org/view/?seq=NZ_CP017100.1.
## GCskew_max <- 1494058
## GCskew_min <- 3858886
## GCskew_min is in atpA, so this is not exactly right.
## use Jeff Barrick's annotation.

## IDEA, perhaps for future work: verify GC skew and replication origin correlation,
## by comparing GC skew against actual replication of origin, based on looking
## at wave pattern in sequencing coverage in my Tet50 genome sequencing samples--
## some of these cultures were still in exponential phase when I sequenced them.
## what about for plasmids? see skewDB.

rotate.NEB5alpha.chr <- function(my.position) {
    #' function to rotate REL606 genome coordinates,
    #' setting oriC at the center of plots
    #' that examine mutation bias over the chromosome.
    ## we want to change coordinates so that c is the new origin.
    GENOME.LENGTH <- 4583637
    midpoint <- GENOME.LENGTH/2
    oriC <- 3886105
    
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
            "data/DH5a-genome-sequencing/DH5a-evolved-populations-and-clones.csv"),
  stringsAsFactors=FALSE)

## This is the key data file for the analysis.
evolved.mutations <- read.csv(
    file.path(projdir,
              "results/DH5a-expt-genome-analysis/evolved_mutations.csv"),
    stringsAsFactors=FALSE) %>%
    mutate(Mbp.coordinate=Position/1000000) %>%
        ## remove the pUC samples from the analysis.
    filter(Plasmid != "pUC") %>%
    ## update the symbol used for intergenic regions on the plasmid
    ## so that it is not replaced by "..." in the figures.
    mutate(Gene = str_replace(Gene,"–/KanR", "-/KanR")) %>%
    ## update the names of the Transposon factor for a prettier plot.
    mutate(Transposon_factor = fct_recode(as.factor(Transposon),
                                          `Tn5+ (TetA++)` = "B30",
                                          `Tn5+ (TetA+)` = "B20",
                                          `Tn5- (TetA++)` = "B59")) %>%
    ## update the names of the Plasmid factor for a prettier plot.
    mutate(Plasmid_factor = fct_recode(as.factor(Plasmid),
                                       `No plasmid` = "None",
                                       p15A = "p15A")) %>%
    mutate(Tet_factor = fct_recode(as.factor(Tet),
                                   `Tet 0` = "0",
                                   `Tet 50` = "50"))



#################################################################################
### Figure 1E: make a matrix plot of genes with mutations in two or more clones.
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
        p15A.hit.count.df <- filter(matrix.data,Plasmid=="pUC") %>%
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

    
    ## make Tet50 panels.
    B59.noPlasmid.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "Tn5- (TetA++)\nNo plasmid\nTet 50")
    ## Remove the gene labels to save space.
    B59.A31.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "Tn5- (TetA++)\np15A\nTet 50") +
        theme(axis.text.y=element_blank())
    
    B20.noPlasmid.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "Tn5+ (TetA+)\nNo plasmid\nTet 50") +
        theme(axis.text.y=element_blank())

    B20.A31.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "Tn5+ (TetA+)\np15A\nTet 50") +
        theme(axis.text.y=element_blank())
    
    B30.noPlasmid.Tet50.matrix.panel <- make.matrix.panel(matrix.data,"Tn5+ (TetA++)\nNo plasmid\nTet 50") +
        theme(axis.text.y=element_blank())

    ## get the legend from the last panel, because this shows all the colors.
    B30.A31.Tet50.matrix.panel <- make.matrix.panel(matrix.data, "Tn5+ (TetA++)\np15A\nTet 50", leg=TRUE) +
        theme(axis.text.y=element_blank())

    Fig.legend <- get_legend(B30.A31.Tet50.matrix.panel)

    ## now remove the legend from the last panel.
    B30.A31.Tet50.matrix.panel <- B30.A31.Tet50.matrix.panel + guides(fill = "none")
    
    ## Using the patchwork library for layout.
    matrix.panels <- B59.noPlasmid.Tet50.matrix.panel +
        B20.noPlasmid.Tet50.matrix.panel +
        B30.noPlasmid.Tet50.matrix.panel +
        B59.A31.Tet50.matrix.panel +
        B20.A31.Tet50.matrix.panel +
        B30.A31.Tet50.matrix.panel +
        Fig.legend +
        plot_layout(nrow = 1)

    ## hack to label x-axis from comments at: https://github.com/thomasp85/patchwork/issues/150
    matrix.panels.grob <- patchwork::patchworkGrob(matrix.panels)
    matrix.figure <- gridExtra::grid.arrange(matrix.panels.grob, left = "", bottom = "Evolved populations")
        
    return(matrix.figure)
}

################################################################################
## Now actually make the figure.

## ALL of these MOB insertions are miniTn5-Tet insertions, either into the KanR gene
## on the plasmid, or into chromosomal genes in the no plasmid treatment.
evolved.MOB <- evolved.mutations %>% filter(Mutation == "MOB") 


gene.level.parallel.mutations <- evolved.mutations %>% group_by(Gene) %>%
summarise(count = n()) %>% filter(count>1) %>% inner_join(evolved.mutations)

parallel.genes <- gene.level.parallel.mutations %>%
    select(Gene, count, Plasmid, Transposon, Tet) %>%
    distinct() %>%
    arrange(desc(count))

## examine genes that show parallelism across Tet 0 and Tet 50.
parallel.genes.across.Tet0.and.Tet50 <- parallel.genes %>%
    select(Gene, Tet) %>%
    distinct() %>%
    group_by(Gene) %>%
    summarize(num.tet.conc.found.in = n()) %>%
    filter(num.tet.conc.found.in == 2)

parallel.mutations.across.Tet0.and.Tet50 <- evolved.mutations %>%
    filter(Gene %in% parallel.genes.across.Tet0.and.Tet50$Gene)

## the intergenic mutations on the plasmid occur in the oriC sequences.
## IMPORTANT TODO: analyze these in a separate figure, to relate to the
## evolved changes in plasmid copy number.
plasmid.origin.muts <- parallel.mutations.across.Tet0.and.Tet50 %>%
    filter(Gene == "–/KanR")

## remove these from the other parallel mutations across the 2 treatments.
## the remaining mutations are not that interesting.
## 3 loci: mrcA, narU, rpsA, yeeJ, 2 parallel mutations in each
filtered.parallel.mutations.across.Tet0.and.Tet50 <- parallel.mutations.across.Tet0.and.Tet50 %>%
    filter(Gene != "–/KanR")


## genes that only show parallelism in Tet 0.
parallel.genes.in.Tet0 <- parallel.genes %>%
    select(Gene, Tet) %>%
    filter(Tet == 0) %>%
    distinct() %>%
    filter(!(Gene %in% parallel.genes.across.Tet0.and.Tet50$Gene))
## not a whole lot that is interesting. maybe the high frequency oxyR mutations.
## Off the top of my head, I believe that is a high level I-modulon regulator.
parallel.mutations.in.only.Tet0 <- evolved.mutations %>%
    filter((Gene %in% parallel.genes.in.Tet0$Gene))


## genes that show parallelism in Tet 50, and all MOB (these are only in Tet50 treatment anyway).
parallel.genes.in.Tet50 <- parallel.genes %>%
    select(Gene, Tet) %>%
    filter(Tet == 50) %>%
   distinct()

parallel.mutations.in.Tet50 <- evolved.mutations %>%
    filter(Gene %in% parallel.genes.in.Tet50$Gene)

Fig1E.data <- full_join(evolved.MOB,
                       filter(parallel.mutations.in.Tet50, Allele != "MOB"))

Fig1E <- MakeMutCountMatrixFigure(Fig1E.data,
                                 show.all=TRUE, ## This is needed to show the MOB insertions too.
                                 use.treatment.hit.sort=FALSE)

## Figure 1E in the ARG duplications manuscript.
Fig1E.outf <- "../results/Fig1E.pdf"
ggsave(Fig1E.outf, Fig1E, height=6, width=12)
ggsave("../results/Fig1E.pdf", Fig1E, height=6, width=12)


