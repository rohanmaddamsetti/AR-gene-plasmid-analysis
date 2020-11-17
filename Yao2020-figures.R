## Yao2020-figures.R by Rohan Maddamsetti.
## This script generates some figures for Yao et al. 2020.

## for dplyr and ggplot2 packages.
library(tidyverse)

## follow directions here to use Calibri font
## and embed in PDF: https://github.com/wch/extrafont/
library(extrafont)


## function for plotting better y-axis labels.
## see solution here for nice scientific notation on axes.
## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
fancy_scientific <- function(x) {
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}


## Make Figure 4B.

## import Yi's data.
Fig4B.data <- read.csv("../data/Yao2020/4B-8figure.csv")

## add columns for error bars and double-check
## mean and sd calculations.
Fig4B.dataframe <- Fig4B.data %>%
    group_by(name, Tet, mean, SD) %>%
    summarize(my.mean = mean(value), my.sd = sd(value)) %>%
    ## error bars are +/- 1 standard deviation.
    mutate(left.error = mean - SD) %>%
    mutate(right.error = mean + SD)

## make the figure object.
Fig4B <- ggplot(data=Fig4B.dataframe, aes(x= Tet,
                                          y = mean,
                                          ymin = left.error,
                                          ymax = right.error)) +
    geom_point(shape = 1, size = 2) +
    geom_line(size = 0.2) +
    geom_errorbar(inherit.aes=TRUE, size = 0.3, width = 0.3) +
    theme_classic() +
    facet_wrap(.~name, nrow=2) + 
    xlab(expression(paste("[Tet] (",mu,"g/mL)"))) +
    ylab(expression('f'[PT]*" frequency")) +
    theme(text=element_text(family="Calibri", size=18)) +
    theme(strip.background = element_blank()) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

## save the figure to file as PDF.
ggsave("../results/Yao2020/Fig4B.pdf", Fig4B,height=3,width=5)
## embed the Calibri font into the PDF itself.
embed_fonts("../results/Yao2020/Fig4B.pdf")

## Make Figure 4C.

## import Yi's data.
Fig4C.data <- read.csv("../data/Yao2020/4C-3figure.csv")

calc.Fig4C.error.bars <- function(Fig4C.df) {

    ## weird equation to calculate error bars, from Yi.
    weird.equation <- function(little.df) {
        antibiotic.row <- filter(little.df,Treatments=="Antibiotic")
        no.antibiotic.row <- filter(little.df,Treatments=="No_Antibiotic")
        SDA <- antibiotic.row$my.sd
        SDL <- no.antibiotic.row$my.sd
        MA <- antibiotic.row$my.mean
        ML <- no.antibiotic.row$my.mean
        weird.eq.for.error <- 0.5*(MA/ML)*sqrt((SDL/ML)^2 + (SDA/MA)^2)
        little.df.with.errorbar <- little.df
        little.df.with.errorbar$error <- weird.eq.for.error
        return(little.df.with.errorbar)
}

    Fig4C.df %>%
        split(list(.$name,.$Tet)) %>%
        map_dfr(.f = weird.equation)
}

calc.Fig4C.freq.values <- function(Fig4C.df) {

    calc.freq <- function(little.df) {
        antibiotic.row <- filter(little.df,Treatments=="Antibiotic")
        no.antibiotic.row <- filter(little.df,Treatments=="No_Antibiotic")
        MA <- antibiotic.row$my.mean
        ML <- no.antibiotic.row$my.mean
        little.df.with.freq <- little.df
        little.df.with.freq$freq <- MA/ML
        return(little.df.with.freq)
}

    Fig4C.df %>%
        split(list(.$name,.$Tet)) %>%
        map_dfr(.f = calc.freq)
}

Fig4C.dataframe <- Fig4C.data %>%
    group_by(name, Tet, Treatments) %>%
    summarize(my.mean = mean(value), my.sd = sd(value)) %>%
    calc.Fig4C.error.bars() %>%
    calc.Fig4C.freq.values() %>%
    select(name, Tet, error, freq) %>%
    distinct() %>%
    ## error bars are +/- 1 the propagated error by the weird equation.
    mutate(left.error = freq - error) %>%
    mutate(right.error = freq + error)

## make the figure object.
Fig4C <- ggplot(data=Fig4C.dataframe, aes(x= Tet,
                                          y = freq,
                                          ymin = left.error,
                                          ymax = right.error)) +
    geom_point(shape = 1, size = 2) +
    geom_line(size = 0.2) +
    geom_errorbar(inherit.aes=TRUE, size = 0.3, width = 0.3) +
    theme_classic() +
    facet_wrap(.~name, nrow=1) + 
    xlab(expression(paste("[Tet] (",mu,"g/mL)"))) +
    ylab(expression('f'[PT]*" frequency")) +
    theme(text=element_text(family="Calibri", size=18)) +
    theme(strip.background = element_blank()) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    scale_y_log10(labels=fancy_scientific)#,
                  ##breaks = scales::extended_breaks(n = 6),
##                       limits = c(0, NA))

## save the figure to file as PDF.
ggsave("../results/Yao2020/Fig4C.pdf", Fig4C,height=3,width=5)
## embed the Calibri font into the PDF itself.
embed_fonts("../results/Yao2020/Fig4C.pdf")
