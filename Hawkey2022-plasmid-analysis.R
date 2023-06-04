## Hawkey2022-plasmid-analysis.R by Rohan Maddamsetti.
## analyze the plasmid copy number results made by
## plasmid-copy-number-pipeline.py.

##     Make a scatterplot of plasmid copy numbers against plasmid length,
##     and color dots by presence of ARGs.
##    
##     based on these results, plasmids with ARGs do NOT have high copy number.

##    The Hawkey et al. 2022 paper specifically focuses on ESBL resistance,
##    so let's focus on beta-lactamases, and can compare to other kinds of resistances in these data.


library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)


## antibiotic-specific keywords.
chloramphenicol.keywords <- "chloramphenicol|Chloramphenicol"
tetracycline.keywords <- "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
MLS.keywords <- "macrolide|lincosamide|streptogramin"
multidrug.keywords <- "Multidrug resistance|multidrug resistance|antibiotic resistance"
beta.lactam.keywords <- "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\\S*"
glycopeptide.keywords <- "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
polypeptide.keywords <- "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
diaminopyrimidine.keywords <- "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
sulfonamide.keywords <- "sulfonamide|Sul1|sul1|sulphonamide"
quinolone.keywords <- "quinolone|Quinolone|oxacin|qnr|Qnr"
aminoglycoside.keywords <- "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
macrolide.keywords <- "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
antimicrobial.keywords <- "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\\S*"


antibiotic.keywords <- paste(chloramphenicol.keywords, tetracycline.keywords, MLS.keywords, multidrug.keywords,
    beta.lactam.keywords, glycopeptide.keywords, polypeptide.keywords, diaminopyrimidine.keywords,
    sulfonamide.keywords, quinolone.keywords, aminoglycoside.keywords, macrolide.keywords, antimicrobial.keywords, sep="|")

################################################################################
## get lengths of all the replicons.
replicon.length.data <- read.csv("../results/Hawkey2022_replicon_lengths.csv")

## get ARG copy number data.
ARG.copy.number.data <- read.csv("../results/Hawkey2022_ARG_copy_numbers.csv") %>%
    mutate(beta.lactam.resistance = ifelse(str_detect(product,beta.lactam.keywords), TRUE, FALSE))

beta.lactam.ARGs <- filter(ARG.copy.number.data, beta.lactam.resistance==TRUE)
non.beta.lactam.ARGs <- filter(ARG.copy.number.data, beta.lactam.resistance==FALSE)

chromosome.plasmid.copy.number.data <- read.csv("../results/Hawkey2022_chromosome_plasmid_copy_numbers.csv") %>%
    full_join(replicon.length.data) %>%
    mutate(has.ARG = ifelse(SeqID %in% ARG.copy.number.data$SeqID, TRUE, FALSE)) %>%
    mutate(has.beta.lactamase = ifelse(SeqID %in% beta.lactam.ARGs$SeqID, TRUE, FALSE)) %>%
    ## 0 == no ARG, 1 == has ARG, 2 == has beta-lactamase.
    mutate(ARG.classification = has.ARG + has.beta.lactamase) %>%
    mutate(ARG.classification = as.factor(ARG.classification)) %>%
    ## remove outlier points with very low coverage.
    filter(CopyNumber > 0.5)


## beta-lactamases have higher copy number compared to other ARGs in these strains.
wilcox.test(beta.lactam.ARGs$CopyNumber, non.beta.lactam.ARGs$CopyNumber,alternative="greater")$p.value

mean(beta.lactam.ARGs$CopyNumber)
mean(non.beta.lactam.ARGs$CopyNumber)

median(beta.lactam.ARGs$CopyNumber)
median(non.beta.lactam.ARGs$CopyNumber)

## Plasmids with ARGs actually have lower copy numbers than
## plasmids without ARGs.

plasmid.copy.number.data <- chromosome.plasmid.copy.number.data %>%
    filter(SeqType == "plasmid") %>%
    arrange(CopyNumber)

beta.lactamase.plasmid.data <- plasmid.copy.number.data %>%
    filter(has.beta.lactamase==TRUE)

no.beta.lactamase.plasmid.data <- plasmid.copy.number.data %>%
    filter(has.beta.lactamase == FALSE)

ARG.plasmid.data <- plasmid.copy.number.data %>%
    filter(has.ARG==TRUE)

no.ARG.plasmid.data <- plasmid.copy.number.data %>%
    filter(has.ARG == FALSE)

mean(beta.lactamase.plasmid.data$CopyNumber)
mean(no.beta.lactamase.plasmid.data$CopyNumber)

mean(ARG.plasmid.data$CopyNumber)
mean(no.ARG.plasmid.data$CopyNumber)

median(beta.lactamase.plasmid.data$CopyNumber)
median(no.beta.lactamase.plasmid.data$CopyNumber)

median(ARG.plasmid.data$CopyNumber)
median(no.ARG.plasmid.data$CopyNumber)




## This will go into the Supplement after polishing.
plasmid.copy.number.plot <- ggplot(plasmid.copy.number.data,
                                   aes(x=log10(replicon_length),y=log10(CopyNumber),
                                       color=ARG.classification)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray")

plasmid.copy.number.plot

