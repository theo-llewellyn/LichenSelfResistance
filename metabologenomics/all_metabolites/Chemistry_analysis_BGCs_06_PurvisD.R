###################################################
# Code to analyse chemical data alongside BGCsc   #
# paired dataset of 83 genomes and metabolomes    #
# SCRIPT 6: Fritz Purvis' D                       #
# Author: Theo Llewellyn                          #
###################################################

library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(ape)
library(picante)
library(vegan)
library(cowplot)
library(ggdendro)
library(phytools)
library(reshape2)
library(gplots)
library(caper)
library(tidyselect)
library(janitor)

setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/")
set.seed(12343434)

##############
#            #
#  BGC DATA  #
#            #
##############

#load in data
bgcs <- read_tsv("BIGSCAPE/BigScape_c0.45_clansoff_Leca82T/Network_Annotations_Full.tsv")
bgcs %>% filter(!grepl('BGC', BGC)) -> filtered_bgcs

#load in data
pks1_clusters <- read_tsv("BIGSCAPE/BigScape_c0.45_clansoff_Leca82T/PKSI_clustering_c0.45.tsv")
nrps_clusters <- read_tsv("BIGSCAPE/BigScape_c0.45_clansoff_Leca82T/NRPS_clustering_c0.45.tsv")
pks_nrp_hybrids <- read_tsv("BIGSCAPE/BigScape_c0.45_clansoff_Leca82T/PKS-NRP_Hybrids_clustering_c0.45.tsv")
pks_other_clusters <- read_tsv("BIGSCAPE/BigScape_c0.45_clansoff_Leca82T/PKSother_clustering_c0.45.tsv")
Terpene_clusters <- read_tsv("BIGSCAPE/BigScape_c0.45_clansoff_Leca82T/Terpene_clustering_c0.45.tsv")
Other_clusters <- read_tsv("BIGSCAPE/BigScape_c0.45_clansoff_Leca82T/Others_clustering_c0.45.tsv")

#change colname
pks1_clusters <- rename(pks1_clusters, BGC = `#BGC Name`, Family.Number = `Family Number`)
nrps_clusters <- rename(nrps_clusters, BGC = `#BGC Name`, Family.Number = `Family Number`)
pks_nrp_hybrids <- rename(pks_nrp_hybrids, BGC = `#BGC Name`, Family.Number = `Family Number`)
pks_other_clusters <- rename(pks_other_clusters, BGC = `#BGC Name`, Family.Number = `Family Number`)
Terpene_clusters <- rename(Terpene_clusters, BGC = `#BGC Name`, Family.Number = `Family Number`)
Other_clusters <- rename(Other_clusters, BGC = `#BGC Name`, Family.Number = `Family Number`)

table_list <- list(pks1_clusters, nrps_clusters, pks_nrp_hybrids, pks_other_clusters, Terpene_clusters, Other_clusters)

#detach plyr after so it doesnt interfere with dplyr
library(plyr)
all_bgc_clusters <- join_all(table_list, type = "full")
detach("package:plyr", unload = TRUE)

#link tbls by BGC column, the nrows drops as it removes any clusters from MiBig
all_bgcs_linked <- inner_join(x = all_bgc_clusters, y = filtered_bgcs)

#change BGC name so it just has the accession code
all_bgcs_linked$BGC <- gsub(pattern = "_scaffold_.*", replacement = "", x = all_bgcs_linked$BGC)

#make a dissimilarity format matrix with taxa as rows and BGCFs a columns
all_bgcs_linked[,c(1,2)] %>% add_column(presence = 1) %>% 
  pivot_wider(names_from = Family.Number,
              values_from = presence,
              values_fn = list(presence = mean),
              values_fill = 0) %>% 
  column_to_rownames(var="BGC") -> dissimilarity_matrix


###############
#             #
#  Chem DATA  #
#             #
###############

#Damiens MolNotator compounds
all_compounds <- read_csv("CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/MolNotator_Damien/node_table.csv")

###Repeat formatting for Negative ion
#removes any ions present in blank samples, merge apothecia and thallus rows
negative_ions_pres_abs <- negative_ions[,c(14:121)] %>% 
  filter(`blk3.mzML` == 0 & `blk.mzML` == 0 & `Blk7.mzML` == 0 & `Blk8.mzML` == 0 & `Blk2.mzML` == 0 & `Blk5.mzML` == 0 & `Blk4.mzML` == 0 & `Blk6.mzML` == 0) %>% 
  dplyr::select(-c(`blk3.mzML`,`blk.mzML`,`Blk7.mzML`,`Blk8.mzML`,`Blk2.mzML`,`Blk5.mzML`,`Blk4.mzML`,`Blk6.mzML`)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA.mzML`, `LIQ109XAAUT.mzML`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA.mzML`, `LIQ146XSPT.mzML`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA.mzML`, `LIQ164TEHYT.mzML`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA.mzML`, `LIQ165TEEXT.mzML`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA.mzML`, `LIQ166TFLAVT.mzML`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA.mzML`, `LIQ240SEVILT.mzML`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA.mzML`, `LIQ69CCART.mzML`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA.mzML`, `LIQ72TVILT.mzML`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA.mzML`, `LIQ78THCRT.mzML`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA.mzML`, `LIQ84UMVET.mzML`)))) %>% 
  dplyr::select(-c(2,3,6,7,11:16,66,67,79,80,83,84,89,90,95,96)) %>% 
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))

# edit column names and change typos in accession codes, also remove any samples not in genome database
colnames(negative_ions_pres_abs) <- gsub(".mzML","", colnames(negative_ions_pres_abs))
negative_ions_pres_abs_filtered <- negative_ions_pres_abs %>% 
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG =  LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

rotated_negative_ions_pres_abs_filtered <- negative_ions_pres_abs_filtered %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>%
  as.data.frame()

rownames(rotated_negative_ions_pres_abs_filtered) <- rotated_negative_ions_pres_abs_filtered$name
rotated_negative_ions_pres_abs_filtered <- rotated_negative_ions_pres_abs_filtered[,-1]
rowSums(rotated_negative_ions_pres_abs_filtered)

#check same length
nrow(dissimilarity_matrix) == nrow(rotated_negative_ions_pres_abs_filtered)


# FORMAT MolNotatori compounds

###Filter for computed neutrals, removes blanks, merge apothecia and thallus, edit accession typos and convert to presence absence, also remove any empty rows
molecules_pres_abs <- all_compounds %>% 
  #remove blanks and only keep neutrals
  filter(status_universal == "neutral" & Blk2 == 0 & Blk4 == 0 & Blk5 == 0 & Blk6 == 0 & Blk7 == 0 & Blk8 == 0 & blk == 0 & blk3 == 0) %>% 
  dplyr::select(-c(1:24,125:145)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))) %>% 
  dplyr::select(-c(2,3,6,7,11:16,66,67,79,80,83,84,89,90,95,96)) %>% 
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG =  LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) %>% 
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))

rotated_molecules_pres_abs <- molecules_pres_abs %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>%
  as.data.frame()

#set taxa as rownames and remove first column and any rows with a total of 0 as cant do dissimilarity on empty rows
rownames(rotated_molecules_pres_abs) <- rotated_molecules_pres_abs$name
rotated_molecules_pres_abs <- rotated_molecules_pres_abs[,-1]
rowSums(rotated_molecules_pres_abs)

#check same length
nrow(dissimilarity_matrix) == nrow(rotated_molecules_pres_abs)
#make new dissimilarity matrix removing empty rows in molecules
row.names.remove <- c("LIQ179PYOC","LIQ230CAOCH","LIQ231CALA")
dissimilarity_matrix_filtered_mol <- dissimilarity_matrix[!(row.names(dissimilarity_matrix) %in% row.names.remove),]

#check same length
nrow(dissimilarity_matrix_filtered_mol) == nrow(rotated_molecules_pres_abs)


##### TREE DATA
Leca82T_tree <- read.tree("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TREE_INFERENCE/IQTree_Leca82T_50p_tAL/concat_OF_tAl_renamed_rooted.treefile")
Leca82T_tree$tip.label <- gsub("'","",Leca82T_tree$tip.label)
Leca82T_tree$node.label <- NULL

##MOLECULES
#subset tree to contain just taxa in molecule matrix and order molecule matrix to be same as tree
presence_absence_ordered_mol <- match.phylo.data(phy = Leca82T_tree, data = rotated_molecules_pres_abs)
presence_absence_ordered_bgcf <- match.phylo.data(phy = Leca82T_tree, data = dissimilarity_matrix_filtered_mol)


#Fritz and Purvis' D to see whether the anthraquinone BGCs or SMs show phylogenetic signal
#need to use the alternative version of phylo.d for it to work in a loop
source("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/04_R_SCRIPTS/Random_useful_functions.R")

#remove BGCs that are absent in all samples
bgc_matrix_named <- presence_absence_ordered_bgcf$data %>% dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
bgc_matrix_named <- cbind(taxa = rownames(bgc_matrix_named), bgc_matrix_named)
write_csv(bgc_matrix_named, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/bgc_matrix_presenceabsence.csv")
write.tree(presence_absence_ordered_bgcf$phy, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/IQTree_Leca80T_SpeciesTree.tre")

#this for loop runs through each column (BGC) and calculates the Purvis D value to test if its has phylogenetic signal
set.seed(15091995)
PhyloD_results_dataframe <- data.frame()
for(i in colnames(bgc_matrix_named)[-1]){
  #print BGC
  print(paste("now analysing BGC",i))
  #runs D calculation
  PhyloD_results <- phylo.d2(phy = presence_absence_ordered_bgcf$phy, 
                             data = bgc_matrix_named, 
                             binvar=i,
                             names.col = taxa)
  #saves BGC name, D-value, p-value whether different from 1, p-value whether different from 0
  PhyloD_output <- c(PhyloD_results$binvar,PhyloD_results$DEstimate[[1]],PhyloD_results$Pval1,PhyloD_results$Pval0)
  PhyloD_results_dataframe <- rbind(PhyloD_results_dataframe, PhyloD_output)
}

#assign distribution
for(i in 1:nrow(PhyloD_results_dataframe)){
  #Brownian = if D is less than 1, pval is dif to 1 and pval is NOT dif to 0
  if(PhyloD_results_dataframe[i,4]>0.05){
    PhyloD_results_dataframe[i,5] <- "Brownian"
  }
  #Random = if D is greater than 0, pval is NOT dif to 1 and pval is dif to 0
  else if(PhyloD_results_dataframe[i,2]>0 & PhyloD_results_dataframe[i,3]>0.05 & PhyloD_results_dataframe[i,4]<=0.05){
    PhyloD_results_dataframe[i,5] <- "Random"
  }
  else{
    PhyloD_results_dataframe[i,5] <- "Other"
  }
}

colnames(PhyloD_results_dataframe) <- c("BGC","D","Pval1","Pval0","Distribution")
write_csv(PhyloD_results_dataframe, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/PhyloD_results_dataframe.csv")
#calculate proportion of BGCs which cant reject null that significantly different from Brownian
length(which(PhyloD_results_dataframe$Pval0 > 0.05))/length(PhyloD_results_dataframe$Pval0)
# 0.929195

#calculate proportion of BGCs which cant reject null that significantly different from random
length(which(PhyloD_results_dataframe$Distribution == "Random"))/length(PhyloD_results_dataframe$Distribution)
# 0.06886518

#show the results for a subset of those BGCs that were important on the pcoa plot for the three subfamilies
spp.scrs.calo <- read_csv("CHAPTER3_METABOLOGENOMICS/significant_BGCs_caloplacoideae_Leca82T.csv")
spp.scrs.telo <- read_csv("CHAPTER3_METABOLOGENOMICS/significant_BGCs_teloschistoideae_Leca82T.csv")
spp.scrs.xantho <- read_csv("CHAPTER3_METABOLOGENOMICS/significant_BGCs_xanthorioideae_Leca82T.csv")

PhyloD_results_dataframe %>% filter(BGC %in% spp.scrs.calo$BGC)
PhyloD_results_dataframe %>% filter(BGC %in% spp.scrs.telo$BGC)
PhyloD_results_dataframe %>% filter(BGC %in% spp.scrs.xantho$BGC)
#same but for anthraquinone ones
PhyloD_results_dataframe %>% filter(BGC %in% c("2452","2448","4606","3986","4586","3063","2518","2527"))

#same for metabolite data
mol_matrix_named <- presence_absence_ordered_mol$data %>% dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
mol_matrix_named$taxa <- rownames(mol_matrix_named)

PhyloD_results_mol_dataframe <- data.frame()
for(i in colnames(mol_matrix_named)[1:length(colnames(mol_matrix_named))-1]){
  print(paste("now analysing metabolite",i))
  PhyloD_results <- phylo.d2(phy = presence_absence_ordered_bgcf$phy, 
                             data = mol_matrix_named, 
                             binvar=i,
                             names.col = taxa)
  PhyloD_output <- c(PhyloD_results$binvar,PhyloD_results$DEstimate[[1]],PhyloD_results$Pval1,PhyloD_results$Pval0)
  PhyloD_results_mol_dataframe <- rbind(PhyloD_results_mol_dataframe, PhyloD_output)
}

#assign distribution
for(i in 1:nrow(PhyloD_results_mol_dataframe)){
  #Brownian = if D is less than 1, pval is dif to 1 and pval is NOT dif to 0
  if(PhyloD_results_mol_dataframe[i,4]>0.05){
    PhyloD_results_mol_dataframe[i,5] <- "Brownian"
  }
  #Random = if D is greater than 0, pval is NOT dif to 1 and pval is dif to 0
  else if(PhyloD_results_mol_dataframe[i,2]>0 & PhyloD_results_mol_dataframe[i,3]>0.05 & PhyloD_results_mol_dataframe[i,4]<=0.05){
    PhyloD_results_mol_dataframe[i,5] <- "Random"
  }
  else{
    PhyloD_results_mol_dataframe[i,5] <- "Other"
  }
}

colnames(PhyloD_results_mol_dataframe) <- c("SM","D","Pval1","Pval0","Distribution")
write_csv(PhyloD_results_dataframe, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/PhyloD_results_mol_dataframe.csv")
#calculate proportion of BGCs which show phylogenetic signal. How many are NOT significantly different from a D value under Brownian motion of evolution across the tree i.e. p-value0 is greater than 0.05
length(which(PhyloD_results_mol_dataframe$Pval0 > 0.05))/length(PhyloD_results_mol_dataframe$Pval0)
# 0.7203166
length(which(PhyloD_results_mol_dataframe$Distribution == "Random"))/length(PhyloD_results_mol_dataframe$Distribution)
# 0.2242744

ggplot(PhyloD_results_mol_dataframe, aes(x = as.numeric(D))) +
  geom_histogram(bins = 50, aes(fill = "SM"), alpha = 0.5) +
  #geom_density(aes(y=..count../10), colour = "steelblue") +
  geom_histogram(data = PhyloD_results_dataframe, bins = 50, aes(fill = "BGC"), alpha = 0.5) +
  #geom_density(data = PhyloD_results_dataframe, aes(y=..count../10),colour = "orange") +
  scale_fill_manual(name='Dataset',breaks=c("SM","BGC"),
                    values=c("SM"='steelblue', "BGC" = "orange")) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5)



#SAME BUT FOR TELOSCHISTACEAE AND USING MCMCTREE PHYLOGENY
##### TREE DATA
load("PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/phy.mean.post.GBM.filt.OLD.RData")
phy.mean.post.GBM.filt.OLD$tip.label <- gsub('L6R42','LIQ80XSP', gsub('LQ337', 'LIQ153CAOA', gsub('LQ341', 'LIQ74CAUR', gsub('LQ352', 'LIQ145TFLA',gsub('L6R38', 'LIQ146XSP',gsub('L6R39', 'LIQ85CALI',gsub('LQ348', 'LIQ92SELA_2',gsub('LQ345', 'LIQ72TVIL',gsub('L6R36', 'LIQ73XASTE-2',gsub('LQ339', 'LIQ75XAME-2',gsub('LQ351', 'LIQ93LETR',gsub('LQ338', 'LIQ69CCAR',gsub('L6R37', 'LIQ70TSP',gsub('L6R41', 'LIQ106LELE',gsub('LQ346', 'LIQ82CAATT',gsub('LQ349', 'LIQ109XAAU_2',gsub('L6R35', 'LIQ143CAAG',gsub('LQ343', 'LIQ76CEHR',gsub('LQ350', 'LIQ78TCHR-2',gsub('L6R40', 'LIQ81XMZF-2',gsub('LQ342', 'LIQ94LETR-2',gsub('LQ344', 'LIQ71TLAC',gsub('Xanthoria_parietina', 'Xanpa2',gsub('LQ340', 'LIQ79DIDIA',gsub('LQ347', 'LIQ84UMVE',phy.mean.post.GBM.filt.OLD$tip.label)))))))))))))))))))))))))

##MOLECULES
#subset tree to contain just taxa in molecule matrix and order molecule matrix to be same as tree
presence_absence_ordered_mol <- match.phylo.data(phy = phy.mean.post.GBM.filt.OLD, data = rotated_molecules_pres_abs)
presence_absence_ordered_bgcf <- match.phylo.data(phy = phy.mean.post.GBM.filt.OLD, data = dissimilarity_matrix_filtered_mol)

#remove BGCs that are absent in all samples
bgc_matrix_named <- presence_absence_ordered_bgcf$data %>% dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
bgc_matrix_named <- cbind(taxa = rownames(bgc_matrix_named), bgc_matrix_named)
write_csv(bgc_matrix_named, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/bgc_matrix_presenceabsence_Teloschistaceae.csv")
write.tree(presence_absence_ordered_bgcf$phy, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/IQTree_Telos69T_SpeciesTree.tre")


#this for loop runs through each column (BGC) and calculates the Purvis D value to test if its has phylogenetic signal
PhyloD_results_dataframe_Telos <- data.frame()
for(i in colnames(bgc_matrix_named)[-1]){
  #print BGC
  print(paste("now analysing BGC",i))
  #runs D calculation
  PhyloD_results <- phylo.d2(phy = presence_absence_ordered_bgcf$phy, 
                             data = bgc_matrix_named, 
                             binvar=i,
                             names.col = taxa)
  #saves BGC name, D-value, p-value whether different from 1, p-value whether different from 0
  PhyloD_output <- c(PhyloD_results$binvar,PhyloD_results$DEstimate[[1]],PhyloD_results$Pval1,PhyloD_results$Pval0)
  PhyloD_results_dataframe_Telos <- rbind(PhyloD_results_dataframe_Telos, PhyloD_output)
}

#assign distribution
for(i in 1:nrow(PhyloD_results_dataframe_Telos)){
  #Brownian = if D is less than 1, pval is dif to 1 and pval is NOT dif to 0
  if(PhyloD_results_dataframe_Telos[i,4]>0.05){
    PhyloD_results_dataframe_Telos[i,5] <- "Brownian"
  }
  #Random = if D is greater than 0, pval is NOT dif to 1 and pval is dif to 0
  else if(PhyloD_results_dataframe_Telos[i,2]>0 & PhyloD_results_dataframe_Telos[i,3]>0.05 & PhyloD_results_dataframe_Telos[i,4]<=0.05){
    PhyloD_results_dataframe_Telos[i,5] <- "Random"
  }
  else{
    PhyloD_results_dataframe_Telos[i,5] <- "Other"
  }
}

colnames(PhyloD_results_dataframe_Telos) <- c("BGC","D","Pval1","Pval0","Distribution")
write_csv(PhyloD_results_dataframe_Telos, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/PhyloD_results_dataframe_Telos_Teloschistaceae.csv")
#calculate proportion of BGCs which cant reject null that significantly different from Brownian
length(which(PhyloD_results_dataframe_Telos$Pval0 > 0.05))/length(PhyloD_results_dataframe_Telos$Pval0)
# 0.9301649
#calculate proportion of BGCs which cant reject null that significantly different from random
length(which(PhyloD_results_dataframe_Telos$Distribution == "Random"))/length(PhyloD_results_dataframe_Telos$Distribution)
# 0.06595538

#show the results for a subset of those BGCs that were important on the pcoa plot for the three subfamilies
PhyloD_results_dataframe_Telos %>% filter(BGC %in% spp.scrs.calo$BGC)
PhyloD_results_dataframe_Telos %>% filter(BGC %in% spp.scrs.telo$BGC)
PhyloD_results_dataframe_Telos %>% filter(BGC %in% spp.scrs.xantho$BGC)
#same but for anthraquinone ones
PhyloD_results_dataframe_Telos %>% filter(BGC %in% c("2452","2448","4606","3986","4586","3063","2518","2527"))

#same for metabolite data
mol_matrix_named <- presence_absence_ordered_mol$data %>% dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
mol_matrix_named$taxa <- rownames(mol_matrix_named)

PhyloD_results_mol_dataframe_Telos <- data.frame()
for(i in colnames(mol_matrix_named)[1:length(colnames(mol_matrix_named))-1]){
  print(paste("now analysing metabolite",i))
  PhyloD_results <- phylo.d2(phy = presence_absence_ordered_bgcf$phy, 
                             data = mol_matrix_named, 
                             binvar=i,
                             names.col = taxa)
  PhyloD_output <- c(PhyloD_results$binvar,PhyloD_results$DEstimate[[1]],PhyloD_results$Pval1,PhyloD_results$Pval0)
  PhyloD_results_mol_dataframe_Telos <- rbind(PhyloD_results_mol_dataframe_Telos, PhyloD_output)
}

#assign distribution to values
for(i in 1:nrow(PhyloD_results_mol_dataframe_Telos)){
  #Brownian = if D is less than 1, pval is dif to 1 and pval is NOT dif to 0
  if(PhyloD_results_mol_dataframe_Telos[i,4]>0.05){
    PhyloD_results_mol_dataframe_Telos[i,5] <- "Brownian"
  }
  #Random = if D is greater than 0, pval is NOT dif to 1 and pval is dif to 0
  else if(PhyloD_results_mol_dataframe_Telos[i,2]>0 & PhyloD_results_mol_dataframe_Telos[i,3]>0.05 & PhyloD_results_mol_dataframe_Telos[i,4]<=0.05){
    PhyloD_results_mol_dataframe_Telos[i,5] <- "Random"
  }
  else{
    PhyloD_results_mol_dataframe_Telos[i,5] <- "Other"
  }
}


colnames(PhyloD_results_mol_dataframe_Telos) <- c("SM","D","Pval1","Pval0", "Distribution")
write_csv(PhyloD_results_dataframe, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/PhyloD_results_mol_dataframe_Telos_Teloschistaceae.csv")
#calculate proportion of BGCs which show phylogenetic signal. How many are NOT significantly different from a D value under Brownian motion of evolution across the tree i.e. p-value0 is greater than 0.05
length(which(PhyloD_results_mol_dataframe_Telos$Pval0 > 0.05))/length(PhyloD_results_mol_dataframe_Telos$Pval0)
# 0.7229551
#proportion are random
length(which(PhyloD_results_mol_dataframe_Telos$Distribution == "Random"))/length(PhyloD_results_mol_dataframe_Telos$Distribution)
# 0.2189974

ggplot(PhyloD_results_mol_dataframe_Telos, aes(x = as.numeric(D))) +
  geom_histogram(bins = 50, aes(fill = "SM"), alpha = 0.5) +
  #geom_density(aes(y=..count../10), colour = "steelblue") +
  geom_histogram(data = PhyloD_results_dataframe, bins = 50, aes(fill = "BGC"), alpha = 0.5) +
  #geom_density(data = PhyloD_results_dataframe, aes(y=..count../10),colour = "orange") +
  scale_fill_manual(name='Dataset',breaks=c("SM","BGC"),
                    values=c("SM"='steelblue', "BGC" = "orange")) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5)



#plot them!
(Dvals_Telos <- ggplot(PhyloD_results_mol_dataframe_Telos, aes(y = as.numeric(D))) +
  geom_violin( aes(x = "SM")) +
  geom_jitter(aes(col = Distribution, x = "SM"), alpha = 0.5, width = 0.2, size = 1) +
  geom_boxplot( aes(x = "SM"), width = 0.1, outlier.shape = NA) +
  #BGC data
  geom_violin(data = PhyloD_results_dataframe_Telos, aes(x = "BGC")) +
  geom_jitter(data = PhyloD_results_dataframe_Telos, aes(col = Distribution, x = "BGC"), alpha = 0.5, width = 0.2, size = 1) +
  geom_boxplot(data = PhyloD_results_dataframe_Telos, aes(x = "BGC"), width = 0.1, outlier.shape = NA) +
  scale_color_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Brownian","NA","Random")) +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  theme(#aspect.ratio=1,
        legend.position = "none") +
  ylab("Fritz & Purvis D") +
  xlab("") +
  ylim(-6,4.5))

(Dvals_all <- ggplot(PhyloD_results_mol_dataframe, aes(y = as.numeric(D))) +
  geom_violin( aes(x = "SM")) +
  geom_jitter(aes(col = Distribution, x = "SM"), alpha = 0.5, width = 0.2, size = 1) +
  geom_boxplot( aes(x = "SM"), width = 0.1, outlier.shape = NA) +
  #BGC data
  geom_violin(data = PhyloD_results_dataframe, aes(x = "BGC")) +
  geom_jitter(data = PhyloD_results_dataframe,aes(col = Distribution, x = "BGC"), alpha = 0.5, width = 0.2, size = 1) +
  geom_boxplot(data = PhyloD_results_dataframe, aes(x = "BGC"), width = 0.1, outlier.shape = NA) +
  scale_color_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Brownian","NA","Random")) +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  theme(#aspect.ratio=1,
        legend.position = c(0.89,0.88), legend.background = element_rect(size = 30)) +
  ylab("Fritz & Purvis D") +
  xlab("") +
  ylim(-6,4.5))

png("FIGURES/FritzPurvisD_BGCvsSM_boxplot_all_Teloschistaceae_brownian.png",  res = 300, width = 3500, height = 2000)
plot_grid(Dvals_all, Dvals_Telos, ncol=2, labels = c("(a)","(b)"))
dev.off()

(Dvals_all_nolegend <- ggplot(PhyloD_results_mol_dataframe, aes(y = as.numeric(D))) +
    geom_violin( aes(x = "SM")) +
    geom_jitter(aes(col = Distribution, x = "SM"), alpha = 0.5, width = 0.2, size = 1) +
    geom_boxplot( aes(x = "SM"), width = 0.1, outlier.shape = NA) +
    #BGC data
    geom_violin(data = PhyloD_results_dataframe, aes(x = "BGC")) +
    geom_jitter(data = PhyloD_results_dataframe,aes(col = Distribution, x = "BGC"), alpha = 0.5, width = 0.2, size = 1) +
    geom_boxplot(data = PhyloD_results_dataframe, aes(x = "BGC"), width = 0.1, outlier.shape = NA) +
    scale_color_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Brownian","NA","Random")) +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    theme(legend.position = "none") +
    ylab("Fritz & Purvis D") +
    xlab("") +
    ylim(-6,4.5))

Dvals_all_prop <- ggplot(PhyloD_results_dataframe,aes(x = "BGC",fill = Distribution)) + 
  geom_bar(position = "fill") +
  geom_text(aes(label = ..count..), stat = "count", position = position_fill(.5)) +
  geom_bar(data = PhyloD_results_mol_dataframe, aes(x = "SM"),position = "fill") +
  geom_text(data = PhyloD_results_mol_dataframe, aes(x = "SM", label = ..count..), stat = "count", position = position_fill(.5)) +
  scale_fill_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Brownian","NA","Random")) +
  ylab("proportion") +
  xlab("")

png("FIGURES/FritzPurvisD_BGCvsSM_boxplot_all_proportions.png",  res = 300, width = 3500, height = 2000)
plot_grid(Dvals_all_nolegend, Dvals_all_prop, ncol=2, labels = c("(a)","(b)"))
dev.off()
