## Yao2020-figures.R by Rohan Maddamsetti.
## This script generates some figures for Yao et al. 2020.

## for dplyr and ggplot2 packages.
library(tidyverse)

## follow directions here to use Calibri font
## and embed in PDF: https://github.com/wch/extrafont/
library(extrafont)

## import Yi's data.
Fig4B.data <- read.csv("../data/Yao2020/4B-8figure.csv")

## add columns for error bars and double-check
## mean and sd calculations.
Fig4B.dataframe <- Fig4B.data %>%
    group_by(name, Tet,mean,SD) %>%
    summarize(my.mean=mean(value),my.sd=sd(value)) %>%
    ## error bars are +/- 1 standard deviation.
    mutate(left.error = mean - SD) %>%
    mutate(right.error = mean + SD)

## make the figure object.
Fig4B <- ggplot(data=Fig4B.dataframe, aes(x=Tet,
                                   y = mean,
                                   ymin = left.error,
                                   ymax = right.error)) +
    geom_point(shape=1) +
    geom_line() +
    geom_errorbar(inherit.aes=TRUE) +
    theme_classic() +
    facet_wrap(.~name,nrow=2) +
    xlab(expression(paste("[Tet] (",mu,"g/mL)"))) +
    ylab(expression('f'[PT]*" frequency")) +
    theme(text=element_text(family="Calibri", size=18)) +
    theme(strip.background = element_blank())

## save the figure to file as PDF.
ggsave("../results/Yao2020/Fig4B.pdf", Fig4B,height=3,width=5)
## embed the Calibri font into the PDF itself.
embed_fonts("../results/Yao2020/Fig4B.pdf")
                    
