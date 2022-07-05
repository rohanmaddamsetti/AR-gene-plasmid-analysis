## qPCR-analysis.R by Rohan Maddamsetti.
## This script analyzes qPCR data using Yi's transposon
## and antibiotic resistance marker system,
## for my first ARG duplication paper.

library(tidyverse)
library(cowplot)
library(forcats)


calc.all.probe.fold.differences <- function(well.df) {
    ## this is a helper function for calculating probe fold differences
    ## per well.

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


calc.only.Tet.Cm.probe.fold.differences <- function(well.df) {
    ## calculate probe fold differences per well.

    ## data analysis using constants calculated from Yi's standard curve calibration.
    T.per.C.constant <- 0.39071847356712

    C <- filter(well.df, probe == 'Cm')$cycle_at_threshold
    T <- filter(well.df, probe == 'Tet')$cycle_at_threshold
    T.per.C <- 2^(C - T)/T.per.C.constant
    
    return.df <- data.frame(Well = unique(well.df$Well),
                            Transposon = unique(well.df$Transposon),
                            Plasmid = unique(well.df$Plasmid),
                            Day = unique(well.df$Day),
                            TetConc = unique(well.df$TetConc),
                            Replicate = unique(well.df$Replicate),
                            transposons.per.chromosome = T.per.C)
    return(return.df)
}

######################################################################
## Figure 4D of the paper.

## Day 1 of experiment, using DH5a-B30 as strain.
april.14.data <- read.csv("../data/qPCR/2022-04-14_DH5a-B30_Tet5-day1-culture_qPCR.csv")


## Use Yi's calibration curve.
april.14.results <- april.14.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))


## Day 2 of experiment, using DH5a-B30 as strain.
april.15.data <- read.csv("../data/qPCR/2022-04-15_DH5a-B30_Tet5-day2-culture_qPCR.csv")


april.15.results <- april.15.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))


## let's join the results and make one figure.
april.14.15.results <- rbind(april.14.results,april.15.results) %>%
    ## update the names of the Plasmid factor for a prettier plot.
    mutate(Plasmid = fct_recode(as.factor(Plasmid),
                      `No plasmid` = "no_plasmid",
                      p15A = "p15A_plasmid",
                      pUC = "pUC_plasmid"))

## Using Yi's calibration produces a more sensible result.
Fig4D <- ggplot(april.14.15.results,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = Replicate,
                           shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    geom_line(size=0.5) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1,2)) + ## set scale for Days.
    scale_shape_discrete(name = "Tetracycline concentration\n(ug/mL)") +
    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per chromosome")

ggsave("../results/Fig4D.pdf", Fig4D, width=7, height=3)

######################################################################

## Day 1 of experiment, using DH5a-B59 as strain.
DH5a.B59.day1.data <- read.csv("../data/qPCR/2022-05-23_DH5a-B59_Tet5-day1-culture_qPCR.csv")

## Use Yi's calibration curve.
DH5a.B59.day1.results <- DH5a.B59.day1.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))


## Day 2 of experiment, using DH5a-B59 as strain.
DH5a.B59.day2.data <- read.csv("../data/qPCR/2022-05-23_DH5a-B59_Tet5-day2-culture_qPCR.csv")

DH5a.B59.day2.results <- DH5a.B59.day2.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))


## let's join the results and make one figure.
DH5a.B59.results <- rbind(DH5a.B59.day1.results, DH5a.B59.day2.results)

## Using Yi's calibration produces a more sensible result.
DH5a.B59.fig <- ggplot(DH5a.B59.results,
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
    ggtitle("No duplications during tetracycline selection in the absence of transposase")

ggsave("../results/DH5a-B59-qPCR-2022-5-22-fig1.pdf", DH5a.B59.fig, width=7, height=3)



## plot transposons per plasmids
DH5a.B59.fig2 <- ggplot(DH5a.B59.results,
                       aes(x = Day,
                           y = transposons.per.plasmid,
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
    ylab("Transposons per plasmid") +
    ggtitle("Transposons per plasmid DH5a+B59")

ggsave("../results/DH5a-B59-qPCR-2022-5-22-fig2.pdf", DH5a.B59.fig2, width=7, height=3)



DH5a.B59.fig3 <- ggplot(DH5a.B59.results,
                       aes(x = Day,
                           y = plasmids.per.chromosome,
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
    ylab("plasmids per chromosome") +
    ggtitle("plasmids per chromosome DH5a+B59")

ggsave("../results/DH5a-B59-qPCR-2022-5-22-fig3.pdf", DH5a.B59.fig3, width=7, height=3)

