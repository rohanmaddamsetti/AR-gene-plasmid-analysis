## qPCR-analysis.R by Rohan Maddamsetti.
## This script analyzes qPCR data using Yi's transposon
## and antibiotic resistance marker system,
## for my first ARG duplication paper.

library(tidyverse)
library(cowplot)


calc.probe.fold.differences <- function(well.df) {
    ## calculate probe fold differences per well.

    ## data analysis using constants calculated from Yi's standard curve calibration.
    T.per.C.constant <- 0.39071847356712
    K.per.C.constant <- 0.58657387456313
    T.per.K.constant <- 0.666102754504912

    C <- filter(well.df, probe == 'Cm')$cycle_at_threshold
    T <- filter(well.df, probe == 'Tet')$cycle_at_threshold
    K <- filter(well.df, probe == 'Kan')$cycle_at_threshold

    T.per.C <- 2^(C - T)/T.per.C.constant
    K.per.C <- 2^(C - K)/K.per.C.constant
    T.per.K <- 2^(K - T)/T.per.K.constant

    ## This is what Yi does on his spreadsheet.
    ## He subtracts 1 in the denominator to account
    ## for the copy of the transposon on the chromosome.
    Yi.transposon.on.plasmid.fraction.calc <- 1 - (K.per.C - T.per.C)/(K.per.C - 1)

    return.df <- data.frame(Well = unique(well.df$Well),
                            Transposon = unique(well.df$Transposon),
                            Plasmid = unique(well.df$Plasmid),
                            Day = unique(well.df$Day),
                            TetConc = unique(well.df$TetConc),
                            Replicate = unique(well.df$Replicate),
                            transposons.per.chromosome = T.per.C,
                            plasmids.per.chromosome = K.per.C,
                            transposons.per.plasmid = T.per.K,
                            Yi.transposon.frac = Yi.transposon.on.plasmid.fraction.calc
                            )
    
    return(return.df)
}

######################################################################

## Day 1 of experiment, using DH5a as strain.
april.14.data <- read.csv("../data/qPCR/2022-04-14_DH5a_Tet5-day1-culture_qPCR.csv")


## Use Yi's calibration curve.
april.14.results <- april.14.data %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))


## Day 2 of experiment, using DH5a as strain.
april.15.data <- read.csv("../data/qPCR/2022-04-15_DH5a_Tet5-day2-culture_qPCR.csv")


april.15.results <- april.15.data %>%
    split(.$Well) %>%
    map_dfr(calc.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))


## let's join the results and make one figure.
april.14.15.results <- rbind(april.14.results,april.15.results)

## Using Yi's calibration produces a more sensible result.
april.14.15.fig <- ggplot(april.14.15.results,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Replicate,
                           shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point() +
    geom_line() +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1,2)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
    guides(color=FALSE) +
    ##geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
    ggtitle("Selection for tetracycline resistance causes transposition-mediated duplications")

ggsave("../results/DH5a-B30-qPCR-2022-4-14-and-15.pdf", april.14.15.fig, width=7, height=3)


