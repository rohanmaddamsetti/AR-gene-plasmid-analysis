## qPCR-analysis.R by Rohan Maddamsetti.
## This script analyzes qPCR data using Yi's transposon
## and antibiotic resistance marker system,
## for my first ARG duplication paper.

library(tidyverse)
library(cowplot)


calc.probe.fold.differences <- function(well.df) {
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


calc.tet.cm.CNV.with.day0.control <- function(df.with.day0.control) {

    m.cm <- 2^1.055 ## calibrated with May 6 2021 standard curve.
    m.tet <- 2^0.904 ## calibrated with May 6 2021 standard curve.

    control.data <- filter(df.with.day0.control, Day == 0)
    tet.control.data <- filter(control.data, probe == "Tet")
    cm.control.data <- filter(control.data, probe == "Cm")

    ## question: maybe I should be using geometric mean here instead?
    mean.control.tet.Cq <- mean(tet.control.data$cycle_at_threshold)
    mean.control.cm.Cq <- mean(cm.control.data$cycle_at_threshold)
    
    beta <- 2^(m.tet * mean.control.tet.Cq - m.cm * mean.control.cm.Cq)

    helper.func <- function(well.df) {
        ## note: this helper uses the beta variable defined above.
        C <- filter(well.df, probe == 'Cm')$cycle_at_threshold
        T <- filter(well.df, probe == 'Tet')$cycle_at_threshold
        T.per.C <- beta * 2^(m.cm * C - m.tet * T)

        return.df <- data.frame(Well = unique(well.df$Well),
                                Transposon = unique(well.df$Transposon),
                                Plasmid = unique(well.df$Plasmid),
                                Day = unique(well.df$Day),
                                TetConc = unique(well.df$TetConc),
                                Replicate = unique(well.df$Replicate),
                                transposons.per.chromosome = T.per.C)
        return(return.df)
    }

    
    results <- df.with.day0.control %>%
        split(.$Well) %>%
        map_dfr(helper.func)
    return(results)
}

######################################################################

## Day 1 of experiment, using DH5a as strain.
april.14.data <- read.csv("../data/qPCR/2022-04-14_DH5a_Tet5-day1-culture_qPCR.csv")

april.14.results <- april.14.data %>%
    calc.tet.cm.CNV.with.day0.control() %>%
    mutate(Replicate = as.factor(Replicate)) %>%
    mutate(TetConc = as.factor(TetConc))
    

april.14.fig <- ggplot(april.14.results,
                       aes(x = Day,
                           y = transposons.per.chromosome,
                           color = TetConc)) +
    facet_wrap(.~Plasmid, scales="free") +
    geom_point() +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_color_discrete(name = "tetracycline concentration\n(ug/mL)") +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    ylab("Transposons per chromosome") +
    ggtitle("Selection for tetracycline resistance causes transposition-mediated duplications")

ggsave("../results/DH5a-qPCR-day1-2022-4-14.pdf", april.14.fig, width=7, height=3)
