library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(janitor)

#CHANGE TAXON $NUMBER IF 79 taxon tree or 117 taxon tree
NUMBER=79

setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/ORTHOFINDER/Results_Leca${NUMBER}T/")

#load in data
Orthogroups <- read_tsv("Orthogroups.GeneCount.tsv")
head(Orthogroups)

#filter for rows where all columns are either 0 or 1 and then where the total number of genes is at least half the number of taxa
#ncol(Orthogroups) - 2 gives the number of taxa as each taxon is a column but there's two extra columns with orthogroup and total which need to be removed
#can change the 0.75 to higher proportion to reduce amount of missing data
Orthogroups %>% filter(across(ends_with(c("renamed","proteins")), ~ .x < 2), Total >= (ncol(Orthogroups) -2) * 0.5) -> single_copy_orthologues_50
nrow(single_copy_orthologues_50)

#save orthogroup column as a text file
write.table(single_copy_orthologues_50$Orthogroup,"Orthologues_Leca${NUMBER}T_50p.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
