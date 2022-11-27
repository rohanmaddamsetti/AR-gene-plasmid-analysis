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
oldFig4D <- ggplot(april.14.15.results,
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

ggsave("../results/old-Fig4D.pdf", oldFig4D, width=7, height=3)

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
    geom_point(size=3) +
    geom_line(size=0.5) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1,2)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
    guides(color= "none") +
    theme(strip.background = element_blank()) +
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
    geom_point(size=3) +
    geom_line(size=0.5) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1,2)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
    guides(color= "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per plasmid") +
    ggtitle("Transposons per plasmid DH5a+B59")

ggsave("../results/DH5a-B59-qPCR-2022-5-22-fig2.pdf", DH5a.B59.fig2, width=7, height=3)

DH5a.B59.fig3 <- ggplot(DH5a.B59.results,
                       aes(x = Day,
                           y = plasmids.per.chromosome,
                           color = Replicate,
                           shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    geom_line(size=0.5) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1,2)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("plasmids per chromosome") +
    ggtitle("plasmids per chromosome DH5a+B59")

ggsave("../results/DH5a-B59-qPCR-2022-5-22-fig3.pdf", DH5a.B59.fig3, width=7, height=3)

######################################################################
## Experiment done on 9/1/2022.
## 1 day qPCR + whole genome sequencing experiment, using K12-B59 and K12-B30 in parallel.
## IMPORTANT: I used a fresh mix of qPCR probes for this experiment.

K12.one.day.data <- read.csv("../data/qPCR/2022-09-01_K12-B59-B30-Tet0-Tet5-Day-1.csv")

## Use Yi's calibration curve.
K12.one.day.results <- K12.one.day.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))

## Using Yi's calibration produces a more sensible result.
K12.one.day.fig <- ggplot(K12.one.day.results,
                          aes(x = Day,
                              y = transposons.per.chromosome,
                              color = Transposon,
                              shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per chromosome")

ggsave("../results/K12-B59-B30-qPCR-2022-9-01-fig1.pdf", K12.one.day.fig, width=7, height = 3)


## plot transposons per plasmids
K12.one.day.fig2 <- ggplot(K12.one.day.results,
                           aes(x = Day,
                               y = transposons.per.plasmid,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per plasmid")

ggsave("../results/K12-B59-B30-qPCR-2022-9-01-fig2.pdf", K12.one.day.fig, width=7, height = 3)

K12.one.day.fig3 <- ggplot(K12.one.day.results,
                           aes(x = Day,
                               y = plasmids.per.chromosome,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("plasmids per chromosome")

ggsave("../results/K12-B59-B30-qPCR-2022-9-01-fig3.pdf", K12.one.day.fig3, width=7, height = 3)

#############################################################################
## Experiment done on 10/21/2022.
## 2 day qPCR + whole genome sequencing experiment, using K12-B107, K12-B111,
## K12-B123, K12-B134, K12-B142, K12-B143 to generalize over native transposons
## with varying transposition kinetics.
## IMPORTANT: I used the same mix of qPCR probes used on 9/1/2022.

K12.generality.two.day.data <- read.csv("../data/qPCR/2022-10-21_K12-native-transposon-Tet5-Day2-culture_qPCR.csv")

## Use Yi's calibration curve.
K12.generality.two.day.results <- K12.generality.two.day.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))

## Using Yi's calibration produces a more sensible result.
K12.generality.fig <- ggplot(K12.generality.two.day.results,
                          aes(x = Day,
                              y = transposons.per.chromosome,
                              color = Transposon,
                              shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
 ##   guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per chromosome")

ggsave("../results/K12-generality-qPCR-2022-10-21-fig1.pdf", K12.generality.fig, width=7, height = 3)


## plot transposons per plasmids
K12.generality.fig2 <- ggplot(K12.generality.two.day.results,
                           aes(x = Day,
                               y = transposons.per.plasmid,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
 ##   guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per plasmid")

ggsave("../results/K12-generality-qPCR-2022-10-21-fig2.pdf", K12.generality.fig2, width=7, height = 3)

K12.generality.fig3 <- ggplot(K12.generality.two.day.results,
                           aes(x = Day,
                               y = plasmids.per.chromosome,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
##    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("plasmids per chromosome")

ggsave("../results/K12-generality-qPCR-2022-10-21-fig3.pdf", K12.generality.fig3, width=7, height = 3)

#############################################################################
## Experiment done on 11/01/2022.
## 7 day qPCR , using K12-B107, K12-B111,
## K12-B123, K12-B134, K12-B142, K12-B143 to generalize over native transposons
## with varying transposition kinetics.
## IMPORTANT: I used a fresh qPCR mix made on 11/1/2022.

## IMPORTANT: This experiment does NOT have Day 0 data. I could redo this experiment, using cells from my streak
## plate, that would probably work if I wanted to redo this experiment.
## for now, use the Tet0 samples as a comparison.

K12.generality.seven.day.data <- read.csv("../data/qPCR/2022-11-01_K12-native-transposon-Tet5-Day7-culture_qPCR.csv")

## Use Yi's calibration curve.
K12.generality.seven.day.results <- K12.generality.seven.day.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))

## Using Yi's calibration produces a more sensible result.
K12.generality.fig <- ggplot(K12.generality.seven.day.results,
                          aes(x = Day,
                              y = transposons.per.chromosome,
                              color = Transposon,
                              shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
 ##   guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per chromosome")

ggsave("../results/K12-generality-qPCR-2022-11-01-fig1.pdf", K12.generality.fig, width=7, height = 3)


## plot transposons per plasmids
K12.generality.fig2 <- ggplot(K12.generality.seven.day.results,
                           aes(x = Day,
                               y = transposons.per.plasmid,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
 ##   guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per plasmid")

ggsave("../results/K12-generality-qPCR-2022-11-01-fig2.pdf", K12.generality.fig2, width=7, height = 3)

K12.generality.fig3 <- ggplot(K12.generality.seven.day.results,
                           aes(x = Day,
                               y = plasmids.per.chromosome,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
##    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("plasmids per chromosome")

ggsave("../results/K12-generality-qPCR-2022-11-01-fig3.pdf", K12.generality.fig3, width=7, height = 3)

#############################################################################
## Experiment done on 11/13/2022.
## 1 day qPCR, using K12-B107, K12-B111,
## K12-B123, K12-B134, K12-B142, K12-B143 to generalize over native transposons
## with varying transposition kinetics.
## I used a qPCR mix made on 11/1/2022.
## Wells in Row G only have the A31 p15A plasmid, and were grown in 20mL LB+Tet5 from a 200uL bottleneck.
## Wells in Row H were grown in 3mL LB+Tet5 from a 30uL bottleneck.


Nov13.data <- read.csv("../data/qPCR/2022-11-13_K12-native-transposon-Tet5-Day1-culture_qPCR.csv")

## Use Yi's calibration curve.
Nov13.results <- Nov13.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))

## Using Yi's calibration produces a more sensible result.
Nov13.fig <- ggplot(Nov13.results,
                          aes(x = Day,
                              y = transposons.per.chromosome,
                              color = Transposon,
                              shape = TetConc)) +
    facet_wrap(Replicate~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
 ##   guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per chromosome")

ggsave("../results/K12-generality-qPCR-2022-11-13-fig1.pdf", Nov13.fig, width=7, height = 6)


## plot transposons per plasmids
Nov13.fig2 <- ggplot(Nov13.results,
                           aes(x = Day,
                               y = transposons.per.plasmid,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
 ##   guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per plasmid")

ggsave("../results/K12-generality-qPCR-2022-11-13-fig2.pdf", Nov13.fig2, width=7, height = 3)

Nov13.fig3 <- ggplot(Nov13.results,
                           aes(x = Day,
                               y = plasmids.per.chromosome,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
##    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("plasmids per chromosome")

ggsave("../results/K12-generality-qPCR-2022-11-13-fig3.pdf", Nov13.fig3, width=7, height = 3)

#############################################################################
## Experiment done on 11/18/2022.
## 2 day qPCR, using K12-B107, K12-B111,
## K12-B123, K12-B109, K12-B110 to generalize over native E. coli transposons
## with varying transposition kinetics.
## I used a qPCR mix made on 11/1/2022.

Nov18.data <- read.csv("../data/qPCR/2022-11-18_K12-native-transposon-Tet5-Day2-culture_qPCR.csv")

## Use Yi's calibration curve.
Nov18.results <- Nov18.data %>%
    split(.$Well) %>%
    map_dfr(calc.all.probe.fold.differences) %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))

## Using Yi's calibration produces a more sensible result.
Nov18.fig <- ggplot(Nov18.results,
                          aes(x = Day,
                              y = transposons.per.chromosome,
                              color = Transposon,
                              shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
 ##   guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per chromosome")

ggsave("../results/K12-generality-qPCR-2022-11-18-fig1.pdf", Nov18.fig, width=9, height = 3)


## plot transposons per plasmids
Nov18.fig2 <- ggplot(Nov18.results,
                           aes(x = Day,
                               y = transposons.per.plasmid,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
 ##   guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("Transposons per plasmid")

ggsave("../results/K12-generality-qPCR-2022-11-18-fig2.pdf", Nov18.fig2, width=7, height = 3)

Nov18.fig3 <- ggplot(Nov18.results,
                           aes(x = Day,
                               y = plasmids.per.chromosome,
                               color = Transposon,
                               shape = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point(size=3) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=c(0,1)) + ## set scale for Days.
    scale_shape_discrete(name = "tetracycline concentration\n(ug/mL)") +
##    guides(color = "none") +
    theme(strip.background = element_blank()) +
    ylab("plasmids per chromosome")

ggsave("../results/K12-generality-qPCR-2022-11-18-fig3.pdf", Nov18.fig3, width=7, height = 3)
