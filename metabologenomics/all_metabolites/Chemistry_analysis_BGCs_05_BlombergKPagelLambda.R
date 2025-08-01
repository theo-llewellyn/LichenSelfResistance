###################################################
# Code to analyse chemical data alongside BGCsc   #
# paired dataset of 83 genomes and metabolomes    #
# SCRIPT 5: Blombergs K & Pagels Lambda           #
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

#mzmine ions
positive_ions <- read_csv("CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/Positive_MS2filter20_quant.csv")
negative_ions <- read_csv("CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/Negative_MS2filter20_quant.csv")
#Damiens MolNotator compounds
all_compounds <- read_csv("CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/MolNotator_Damien/node_table.csv")
#GNPS clusters
positive_gnps <- read_csv("CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/KewGardens_Telos_POS_cos0.55_minmatfrag4_metadata_ms2filter20.csv")
negative_gnps <- read_csv("CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/KewGardens_Telos_NEG_cos0.55_minmatfrag4_metadata_ms2filter20.csv")


# FORMAT MZMine ions

#reformat column names to make filtering a bit easier
colnames(positive_ions) <- gsub(" Peak area", "", x = colnames(positive_ions))
colnames(negative_ions) <- gsub(" Peak area", "", x = colnames(negative_ions))

#removes any ions present in blank samples, merge apothecia and thallus rows
positive_ions_pres_abs <- positive_ions[,c(14:121)] %>% 
  filter(`blk3-1.mzML` == 0 & `blk.mzML` == 0 & `Blk7.mzML` == 0 & `Blk8.mzML` == 0 & `Blk2.mzML` == 0 & `Blk5.mzML` == 0 & `Blk4.mzML` == 0 & `Blk6.mzML` == 0) %>% 
  dplyr::select(-c(`blk3-1.mzML`,`blk.mzML`,`Blk7.mzML`,`Blk8.mzML`,`Blk2.mzML`,`Blk5.mzML`,`Blk4.mzML`,`Blk6.mzML`)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA.mzML`, `LIQ109XAAUT.mzML`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA.mzML`, `LIQ146XSPT.mzML`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA.mzML`, `LIQ164TEHYT.mzML`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA.mzML`, `LIQ165TEEXT.mzML`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA.mzML`, `LIQ166TFLAVT.mzML`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA.mzML`, `LIQ240SEVILT.mzML`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA.mzML`, `LIQ69CCART.mzML`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA.mzML`, `LIQ72TVILT.mzML`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA.mzML`, `LIQ78THCRT.mzML`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA.mzML`, `LIQ84UMVET.mzML`)))) %>% 
  dplyr::select(-c(2,3,6,7,11:16,66,67,79,80,83,84,89,90,95,96)) %>% 
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))

# edit column names and change typos in accession codes, also remove any samples not in genome database
colnames(positive_ions_pres_abs) <- gsub(".mzML","", colnames(positive_ions_pres_abs))
positive_ions_pres_abs_filtered <- positive_ions_pres_abs %>% 
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG =  LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) %>% 
  mutate_if(is.numeric, ~1 * (. > 0))


rotated_positive_ions_pres_abs_filtered <- positive_ions_pres_abs_filtered %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>%
  as.data.frame()

rownames(rotated_positive_ions_pres_abs_filtered) <- rotated_positive_ions_pres_abs_filtered$name
rotated_positive_ions_pres_abs_filtered <- rotated_positive_ions_pres_abs_filtered[,-1]
rowSums(rotated_positive_ions_pres_abs_filtered)

#check same length
nrow(dissimilarity_matrix) == nrow(rotated_positive_ions_pres_abs_filtered)

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


##FORMAT GNPS DATA
#make column names more reader friendly
colnames(positive_gnps) <- gsub("GNPSGROUP:","", colnames(positive_gnps))  %>% gsub(".mzML","", .)
positive_gnps %>%
  #only keep the subcluster numbers and the column that shows which accessions its found in
  dplyr::select(c(34,53,91:190)) %>%
  #remove clusters present in blanks
  filter(BLANK == 0) %>% 
  #remove blank column
  dplyr::select(-c(2)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  #merge columns for apothecia and thallus of same sample
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))) %>%
  #remove thallus and apothecia columns to leave a single column per accession
  dplyr::select(-c(3,4,7,8,12:17,67,68,80,81,84,85,90,91,96,97)) %>%
  #correct typos
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) -> positive_gnps_full

#aggregate the subnetworks
positive_gnps_full %>%
  filter(componentindex != -1) %>% 
  aggregate(. ~  componentindex, data = ., sum) %>%
  tibble() %>%
  mutate_at(vars(-('componentindex')), ~1 * (. > 0)) %>%
  column_to_rownames('componentindex') -> positive_gnps_aggregates

#isolate the singleton nodes
positive_gnps_full %>%
  filter(componentindex == -1) %>%
  dplyr::select(-1) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) -> positive_gnps_singletons

#merge the aggregates and singletons
positive_gnps_agg_sing <- rbind(positive_gnps_aggregates,positive_gnps_singletons)

positive_gnps_agg_sing %>% 
  rownames_to_column() %>%
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>%
  as.data.frame() -> positive_gnps_pres_abs

rownames(positive_gnps_pres_abs) <- positive_gnps_pres_abs$name
#convert to binary
positive_gnps_pres_abs <- positive_gnps_pres_abs[,-1]

rowSums(positive_gnps_pres_abs)
#remove any gnps clusters absent in all remaining genomes
positive_gnps_pres_abs <- positive_gnps_pres_abs[,colSums(positive_gnps_pres_abs) > 0]


#reformat GNPS data, remove blank samples, convert to presence-absence matrix, set accessions as rownames, remove accession column and then convert values to 0 or 1
colnames(negative_gnps) <- gsub("GNPSGROUP:","", colnames(negative_gnps))  %>% gsub(".mzML","", .)
negative_gnps_full <- negative_gnps %>%
  #only keep the subcluster numbers and the column that shows which accessions its found in
  dplyr::select(c(34,53,91:190)) %>%
  #remove clusters present in blanks
  filter(BLANK == 0) %>% 
  #remove blank column
  dplyr::select(-c(2)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))) %>%
  #remove thallus and apothecia columns to leave a single column per accession
  dplyr::select(-c(3,4,7,8,12:17,67,68,80,81,84,85,90,91,96,97)) %>%
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP))

#aggregate the subnetworks
negative_gnps_full %>%
  filter(componentindex != -1) %>% 
  aggregate(. ~  componentindex, data = ., sum) %>%
  tibble() %>%
  mutate_at(vars(-('componentindex')), ~1 * (. > 0)) %>%
  column_to_rownames('componentindex') -> negative_gnps_aggregates

#isolate the singleton nodes
negative_gnps_full %>%
  filter(componentindex == -1) %>%
  dplyr::select(-1) %>%
  mutate_if(is.numeric, ~1 * (. > 0)) -> negative_gnps_singletons

#merge the aggregates and singletons
negative_gnps_agg_sing <- rbind(negative_gnps_aggregates,negative_gnps_singletons)

negative_gnps_agg_sing %>% 
  rownames_to_column() %>%
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>%
  as.data.frame() -> negative_gnps_pres_abs

rownames(negative_gnps_pres_abs) <- negative_gnps_pres_abs$name
#convert to binary
negative_gnps_pres_abs <- negative_gnps_pres_abs[,-1]

rowSums(negative_gnps_pres_abs)
#remove any gnps clusters absent in all remaining genomes
negative_gnps_pres_abs <- negative_gnps_pres_abs[,colSums(negative_gnps_pres_abs) > 0]
#remove any genomes absent in all remaining gnps clusters
negative_gnps_pres_abs <- negative_gnps_pres_abs[rowSums(negative_gnps_pres_abs) > 0,]


#check same length
nrow(dissimilarity_matrix) == nrow(negative_gnps_pres_abs)

#merge
molspersample <- as.data.frame(colSums(molecules_pres_abs))
molspersample$accession <- rownames(molspersample)
pospersample <- as.data.frame(rowSums(rotated_positive_ions_pres_abs_filtered))
pospersample$accession <- rownames(pospersample)
negpersample <- as.data.frame(rowSums(rotated_negative_ions_pres_abs_filtered))
negpersample$accession <- rownames(negpersample)
posGNPSpersample <- as.data.frame(rowSums(positive_gnps_pres_abs))
posGNPSpersample$accession <- rownames(posGNPSpersample)
negGNPSpersample <- as.data.frame(rowSums(negative_gnps_pres_abs))
negGNPSpersample$accession <- rownames(negGNPSpersample)
bgcspersample<- as.data.frame(rowSums(dissimilarity_matrix_filtered_mol))
bgcspersample$accession <- rownames(bgcspersample)

merged_tables <- inner_join(x = molspersample, y = bgcspersample)
merged_tables <- inner_join(x = merged_tables, y = pospersample)
merged_tables <- inner_join(x = merged_tables, y = negpersample)
merged_tables <- inner_join(x = merged_tables, y = posGNPSpersample)
merged_tables <- inner_join(x = merged_tables, y = negGNPSpersample)
colnames(merged_tables) <- c("molecules","accession","bgcs","pos_ions","neg_ions","pos_GNPS","neg_GNPS")

#
melted_merged_tables <- merged_tables %>% melt



##### TREE DATA
Leca82T_tree <- read.tree("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TREE_INFERENCE/IQTree_Leca82T_50p_tAL/concat_OF_tAl_renamed_rooted.treefile")
Leca82T_tree$tip.label <- gsub("'","",Leca82T_tree$tip.label)
Leca82T_tree$node.label <- NULL

##MOLECULES
#subset tree to contain just taxa in molecule matrix and order molecule matrix to be same as tree
presence_absence_ordered_mol <- match.phylo.data(phy = Leca82T_tree, data = rotated_molecules_pres_abs)
presence_absence_ordered_bgcf <- match.phylo.data(phy = Leca82T_tree, data = dissimilarity_matrix_filtered_mol)


#Blomberg's K and Pagel's lambda to test effect of phylogeny on no. BGCs and SMs. It does a stastical test to calculate the probability that these trait values could be random (K or lambda =0). It shuffles the trait values on the tree and then recalculates K/lambda 10,000 times and then compares to observed K/lambda. It also tests whether observed values for K/lambda are less than expected under Brownian evolution. For K we do this by generating 10,000 datasets under brownian model on the Leca tree and then calculate K. We then count proportion of times simulated values were lower than observed K. For lambda we compute the likelihood that lambda = 1 and then compare the maximum likelihood estimate of lambda from our data. We compare using LRT with 1df.

library(geiger)

traitVector <- c()
KVector <- c()
KpVector <- c()
KpBMVector <- c()
lambdaVector <- c()
lambdapVector <- c()
lambda_Pval_BM_Vector <- c()

#simulate datasets
nullX <- fastBM(presence_absence_ordered_bgcf$phy, nsim = 10000)
#test phylogenetic signal
nullK <- apply(nullX, 2, phylosig, tree = presence_absence_ordered_bgcf$phy)

for(i in c("bgcs","molecules","pos_ions","neg_ions","pos_GNPS","neg_GNPS")){
  #make vector of trait values named by accession 
  trait <- setNames(merged_tables[[i]], merged_tables$accession)
  #calculate blombergs K and test difference to random (K=0)
  K <- phylosig(tree = presence_absence_ordered_bgcf$phy, x = trait, method = 'K', test = T, nsim = 10000)
  #calculate Pagels lambda and test difference to random 
  lambda <- phylosig(tree = presence_absence_ordered_bgcf$phy, x = trait, method = 'lambda', test = T)
  #make a vector with name of trait
  traitVector <- c(traitVector,i)
  #save K and lamda values and the probability theyre different from 0
  KVector <- c(KVector,K$K)
  KpVector <- c(KpVector,K$P)
  KpBMVector <- c(KpBMVector,mean(nullK<=K$K))
  lambdaVector <- c(lambdaVector,lambda$lambda)
  lambdapVector <- c(lambdapVector,lambda$P)
  #do LRT of lambda = 1 ie Brownian evolution
  LR <- -2*(lambda$lik(1)-lambda$logL)
  #calculate p-value of observed being different to Brownian
  lambda_Pval_BM <- pchisq(LR, df = 1, lower.tail = FALSE)
  lambda_Pval_BM_Vector <- c(lambda_Pval_BM_Vector,lambda_Pval_BM)
}

K_lambda_results <- data.frame(traitVector,round(KVector,4),round(KpVector,4),round(KpBMVector,4),round(lambdaVector,4),round(lambdapVector,4),round(lambda_Pval_BM_Vector,4))
write_csv(K_lambda_results,"~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/K_lambda_metabologenomics_Leca82T.csv")


#compare fit of three models for each trait
aic_results <- c()
for(i in c("bgcs","molecules","pos_ions","neg_ions","pos_GNPS","neg_GNPS")){
  #make vector of trait values named by accession 
  trait <- setNames(merged_tables[[i]], merged_tables$accession)
  #fit BM model
  fit_BM <- fitContinuous(presence_absence_ordered_bgcf$phy, trait, model = "BM")
  #fit early burst model
  fit_EB <- fitContinuous(presence_absence_ordered_bgcf$phy, trait, model = "EB")
  #fit ornstein-uhlenbeck model
  fit_OU <- fitContinuous(presence_absence_ordered_bgcf$phy, trait, model = "OU", 
                          bounds = list(alpha=c(0,100)))
  #see which model fits best
  aic <- setNames(c(AIC(fit_BM),AIC(fit_EB),AIC(fit_OU)),c("BM","EB","OU"))
  aic_weighted <- setNames(aic.w(aic),c("BM_w","EB_w","OU_w"))
  aic_results <- rbind(aic_results, c(aic, aic_weighted))
}
aic_results <- as.data.frame(aic_results)
aic_results$trait <- c("bgcs","molecules","pos_ions","neg_ions","pos_GNPS","neg_GNPS")

write_csv(aic_results, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/AIC_models_metabologenomics_Leca82T.csv")



#SAME AS ABOVE BUT JUST FOR TELOSCHISTACEAE AND USING MCMCTREE DATED PHYLOGENY


load("PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/phy.mean.post.GBM.filt.OLD.RData")
phy.mean.post.GBM.filt.OLD$tip.label <- gsub('L6R42','LIQ80XSP', gsub('LQ337', 'LIQ153CAOA', gsub('LQ341', 'LIQ74CAUR', gsub('LQ352', 'LIQ145TFLA',gsub('L6R38', 'LIQ146XSP',gsub('L6R39', 'LIQ85CALI',gsub('LQ348', 'LIQ92SELA_2',gsub('LQ345', 'LIQ72TVIL',gsub('L6R36', 'LIQ73XASTE-2',gsub('LQ339', 'LIQ75XAME-2',gsub('LQ351', 'LIQ93LETR',gsub('LQ338', 'LIQ69CCAR',gsub('L6R37', 'LIQ70TSP',gsub('L6R41', 'LIQ106LELE',gsub('LQ346', 'LIQ82CAATT',gsub('LQ349', 'LIQ109XAAU_2',gsub('L6R35', 'LIQ143CAAG',gsub('LQ343', 'LIQ76CEHR',gsub('LQ350', 'LIQ78TCHR-2',gsub('L6R40', 'LIQ81XMZF-2',gsub('LQ342', 'LIQ94LETR-2',gsub('LQ344', 'LIQ71TLAC',gsub('Xanthoria_parietina', 'Xanpa2',gsub('LQ340', 'LIQ79DIDIA',gsub('LQ347', 'LIQ84UMVE',phy.mean.post.GBM.filt.OLD$tip.label)))))))))))))))))))))))))


##MOLECULES
#subset tree to contain just taxa in molecule matrix and order molecule matrix to be same as tree
presence_absence_ordered_mol <- match.phylo.data(phy = phy.mean.post.GBM.filt.OLD, data = rotated_molecules_pres_abs)
presence_absence_ordered_bgcf <- match.phylo.data(phy = phy.mean.post.GBM.filt.OLD, data = dissimilarity_matrix_filtered_mol)


#Blomberg's K and Pagel's lambda to test effect of phylogeny on no. BGCs and SMs.

traitVector <- c()
KVector <- c()
KpVector <- c()
KpBMVector <- c()
lambdaVector <- c()
lambdapVector <- c()
lambda_Pval_BM_Vector <- c()

#simulate datasets
nullX <- fastBM(presence_absence_ordered_bgcf$phy, nsim = 10000)
#test phylogenetic signal
nullK <- apply(nullX, 2, phylosig, tree = presence_absence_ordered_bgcf$phy)

for(i in c("bgcs","molecules","pos_ions","neg_ions","pos_GNPS","neg_GNPS")){
  #make vector of trait values named by accession 
  trait <- setNames(merged_tables[[i]], merged_tables$accession)
  #calculate blombergs K and test difference to random (K=0)
  K <- phylosig(tree = presence_absence_ordered_bgcf$phy, x = trait, method = 'K', test = T, nsim = 10000)
  #calculate Pagels lambda and test difference to random 
  lambda <- phylosig(tree = presence_absence_ordered_bgcf$phy, x = trait, method = 'lambda', test = T)
  #make a vector with name of trait
  traitVector <- c(traitVector,i)
  #save K and lamda values and the probability theyre different from 0
  KVector <- c(KVector,K$K)
  KpVector <- c(KpVector,K$P)
  KpBMVector <- c(KpBMVector,mean(nullK<=K$K))
  lambdaVector <- c(lambdaVector,lambda$lambda)
  lambdapVector <- c(lambdapVector,lambda$P)
  #do LRT of lambda = 1 ie Brownian evolution
  LR <- -2*(lambda$lik(1)-lambda$logL)
  #calculate p-value of observed being different to Brownian
  lambda_Pval_BM <- pchisq(LR, df = 1, lower.tail = FALSE)
  lambda_Pval_BM_Vector <- c(lambda_Pval_BM_Vector,lambda_Pval_BM)
}

K_lambda_results <- data.frame(traitVector,round(KVector,4),round(KpVector,4),round(KpBMVector,4),round(lambdaVector,4),round(lambdapVector,4),round(lambda_Pval_BM_Vector,4))
write_csv(K_lambda_results,"~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/K_lambda_metabologenomics_Teloschistaceae.csv")


#compare fit of three models for each trait
aic_results <- c()
for(i in c("bgcs","molecules","pos_ions","neg_ions","pos_GNPS","neg_GNPS")){
  #make vector of trait values named by accession 
  trait <- setNames(merged_tables[[i]], merged_tables$accession)
  #fit BM model
  fit_BM <- fitContinuous(presence_absence_ordered_bgcf$phy, trait, model = "BM")
  #fit early burst model
  fit_EB <- fitContinuous(presence_absence_ordered_bgcf$phy, trait, model = "EB")
  #fit ornstein-uhlenbeck model
  fit_OU <- fitContinuous(presence_absence_ordered_bgcf$phy, trait, model = "OU", 
                          bounds = list(alpha=c(0,100)))
  #see which model fits best
  aic <- setNames(c(AIC(fit_BM),AIC(fit_EB),AIC(fit_OU)),c("BM","EB","OU"))
  aic_weighted <- setNames(aic.w(aic),c("BM_w","EB_w","OU_w"))
  aic_results <- rbind(aic_results, c(aic, aic_weighted))
}
aic_results <- as.data.frame(aic_results)
aic_results$trait <- c("bgcs","molecules","pos_ions","neg_ions","pos_GNPS","neg_GNPS")

write_csv(aic_results, "~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/CHEMISTRY/AIC_models_metabologenomics_Teloschistaceae.csv")
