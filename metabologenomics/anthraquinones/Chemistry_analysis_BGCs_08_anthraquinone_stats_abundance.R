###################################################
# Code to analyse chemical data alongside BGCsc   #
# paired dataset of 83 genomes and metabolomes    #
# SCRIPT 8: Anthraquinone abundance               #
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
library(ggtree)
library(aplot)

setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/")

###############
#             #
#  Chem DATA  #
#             #
###############

#FORMAT GNPS DATA
positive_gnps <- read_csv('CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/Filtered_Positive_MS2filter20_SIRIUS_summary.tsv/GNPS_table_wSIRIUS_compoundnames.csv')

#make column names more reader friendly
colnames(positive_gnps) <- gsub("GNPSGROUP:","", colnames(positive_gnps))  %>% gsub(".mzML","", .)
positive_gnps %>%
  #only keep the subcluster numbers and the column that shows which accessions its found in
  dplyr::select(c(265,291,5,109:208,71)) %>%
  #remove clusters present in blanks
  filter(BLANK == 0) %>% 
  #remove blank column
  dplyr::select(-c(BLANK)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  #merge columns for apothecia and thallus of same sample
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`)))/2, LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`)))/2, LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`)))/2, LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`)))/2, LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`)))/2, LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`)))/2, LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`)))/2, LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`)))/2, LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`)))/2, LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))/2) %>%
  #remove thallus and apothecia columns to leave a single column per accession
  dplyr::select(-c(3,4,7,8,12:17,67,68,80,81,84,85,90,91,96,97)) %>%
  #correct typos
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  #remove samples without annotated compound
  filter_at(vars(`Analog:Compound_Name`, structure_name),any_vars(!is.na(.))) %>%
  #remove samples not in genome dataset
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) -> positive_gnps_full

#merge GNPS and SIRIUS names
positive_gnps_full$compound_name <- apply(positive_gnps_full[,c(1,75)], 1, function(x) x[!is.na(x)][1])
#format compound names
positive_gnps_full$compound_name <- gsub("Massbank:F.* ","", positive_gnps_full$compound_name)  %>% gsub("usnic","Usnic", .)

rotated_GNPS_abundance <- positive_gnps_full %>%
  #just anthraquinones
  slice(c(1,12,17,24,34,38,49,56,57,76,81,101)) %>%
  dplyr::select(-c(1,74,75)) %>% 
  #sum nodes with same compound ID
  group_by(compound_name) %>%
  summarise_all(.funs = sum, na.rm = TRUE)

#NEGATIVE
negative_gnps <- read_csv('CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/Filtered_negative_MS2filter20_SIRIUS_summary.tsv/GNPS_table_wSIRIUS_compoundnames.csv')

#make column names more reader friendly
colnames(negative_gnps) <- gsub("GNPSGROUP:","", colnames(negative_gnps))  %>% gsub(".mzML","", .)
negative_gnps %>%
  #only keep the subcluster numbers and the column that shows which accessions its found in
  dplyr::select(c(264,290,5,109:208,71)) %>%
  #remove clusters present in blanks
  filter(BLANK == 0) %>% 
  #remove blank column
  dplyr::select(-c(BLANK)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  #merge columns for apothecia and thallus of same sample
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`)))/2, LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`)))/2, LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`)))/2, LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`)))/2, LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`)))/2, LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`)))/2, LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`)))/2, LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`)))/2, LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`)))/2, LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))/2) %>%
  #remove thallus and apothecia columns to leave a single column per accession
  dplyr::select(-c(3,4,7,8,12:17,67,68,80,81,84,85,90,91,96,97)) %>%
  #correct typos
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  #remove samples without annotated compound
  filter_at(vars(`Analog:Compound_Name`, structure_name),any_vars(!is.na(.))) %>%
  #remove samples not in genome dataset
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) -> negative_gnps_full

#merge GNPS and SIRIUS names
negative_gnps_full$compound_name <- apply(negative_gnps_full[,c(1,75)], 1, function(x) x[!is.na(x)][1])
#format compound names and synonyms
negative_gnps_full$compound_name <- gsub("Canarion","Canarione", negative_gnps_full$compound_name) %>% gsub("2-chloro-1,3,8-trihydroxy-6-\\(hydroxymethyl\\)anthracene-9,10-dione","7-Chlorocitreosein",.)

rotated_GNPS_abundance_neg <- negative_gnps_full %>%
  #just anthraquinones
  slice(c(3,9,16,24,27,35,42,43,44,51,53,60,69,70)) %>% 
  dplyr::select(-c(1,74,75)) %>%
  #sum nodes with same compound ID
  {
    bind_rows(
      slice(., -c(3, 14)),
      slice(., c(3, 14)) %>%
        summarise(
          across(where(is.numeric), sum),
          across(where(is.character), ~ "7-Chloroemodin")
        )
    )
  } %>%
  {
    bind_rows(
      slice(., -c(1, 11)),
      slice(., c(1, 11)) %>%
        summarise(
          across(where(is.numeric), sum),
          across(where(is.character), ~ "7-Chlorocitreosein")
        )
    )
  }
  

#merge positive and negative, if the same compound is found in both take the average of the two as they are pseudoreplicates
merged_df <- bind_rows(rotated_GNPS_abundance, rotated_GNPS_abundance_neg) %>%
  group_by(compound_name) %>%
  summarise_all(.funs = mean, na.rm = TRUE) %>%
  pivot_longer(-compound_name) %>% 
  pivot_wider(names_from=compound_name, values_from=value) %>%
  as.data.frame()

#set taxa as rownames and remove first column and any rows with a total of 0 as cant do dissimilarity on empty rows
rownames(merged_df) <- merged_df$name
merged_df <- merged_df[,-1]

#ALPHA DIVERSITY
#########################
#make vector assigning the rownames to which order they are in
metadata <- read_csv("CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/metadata_Leca82T.csv")

#function to extract variable and make a list of colours
taxon_colour_function <- function(colours,taxa_list, variable){
  taxon_vector <- c()
  #for each taxon on the list of taxa
  for(taxon in taxa_list){
    #find which group its in at a particular taxonomic rank
    group <- subset(metadata, accession == taxon)[,variable][[1]]
    taxon_vector <- c(taxon_vector,group)
  }
  taxon_vector <- as.factor(taxon_vector)
  #make a vector of colours using the colours given and the taxa list
  colour_vector <- colours[taxon_vector]
  #make a vector given the levels of the taxa
  level_vector <- levels(taxon_vector)
  #combine into data frame to be able to access
  output <- data.frame(taxon_vector,colour_vector)
}

col_vec1 <- c("#5d6941","black","#d99627","#96adcb","#d1d7b3","pink","#a74207", "#666f6e")
colours_alpha_subfam <- taxon_colour_function(colours = col_vec1,
                                              taxa_list = rownames(merged_df),
                                              variable = "subfamily")
#generate variable of whether Teloschistaceae or not
colours_alpha_subfam$is_Telos <- ifelse(colours_alpha_subfam$taxon_vector == 'Caloplacoideae'|
                                        colours_alpha_subfam$taxon_vector == 'Teloschistoideae'|
                                        colours_alpha_subfam$taxon_vector == 'Xanthorioideae', 
                                        'Teloschistaceae', 'Outgroup') %>% replace_na("Outgroup")
colours_alpha_subfam$Telos_colours <- ifelse(colours_alpha_subfam$is_Telos == 'Teloschistaceae', "#a74207", 'grey')
colours_alpha_subfam$taxa <- rownames(merged_df)

####
anth_richness <- tibble(rownames(merged_df),diversity(merged_df, index = "shannon"),colours_alpha_subfam)
colnames(anth_richness) <- c('taxon','anthraquinone_richness','subfamily','colour','is_Telos')


#check if normal
shapiro.test(anth_richness$anthraquinone_richness[anth_richness$is_Telos == 'Teloschistaceae'])
shapiro.test(anth_richness$anthraquinone_richness[anth_richness$is_Telos == 'Outgroup'])
#if homoskedasticity
library(rstatix)
levene_test(anthraquinone_richness ~ is_Telos, data = anth_richness)
#Welch t-test for unequal variance
t.test(anthraquinone_richness ~ is_Telos, data = anth_richness, 
       alternative ="less",
       var.equal = FALSE)

#anthraquinones not more diverse in Telos


#test whether each anthraquinone is over-represented compared to other compounds in Telos
merged_df %>%
  rownames_to_column() -> merged_df_tibble
#make column whether telos or not with telos as the reference category level
merged_df_tibble$is_Telos <- colours_alpha_subfam$is_Telos
merged_df_tibble$is_Telos <- factor(merged_df_tibble$is_Telos)
merged_df_tibble$is_Telos <- relevel(merged_df_tibble$is_Telos, ref = 'Teloschistaceae')


#for every anthraquinone
compound_cols <- setdiff(colnames(merged_df_tibble), c("rowname", "is_Telos"))
merged_df_tibble$total_abundance <- rowSums(merged_df_tibble[, compound_cols])
results <- list()

# Loop through each compound and calculate its ratio
for (compound in compound_cols) {
  
  # Compute the ratio of the compound to total abundance
  merged_df_tibble[[paste0(compound, "_ratio")]] <- merged_df_tibble[[compound]] / merged_df_tibble$total_abundance
  
  # Perform a t-test to check if the compound's ratio differs between groups
  test_result <- t.test(merged_df_tibble[[paste0(compound, "_ratio")]] ~ merged_df_tibble$is_Telos,
                        alternative = 'greater')
  
  # Store the result (p-value) in the results list
  results[[compound]] <- test_result$p.value
}

# Convert the results into a data frame
results_df <- data.frame(
  Compound = names(results),
  P_Value = unlist(results),
  stringsAsFactors = FALSE
)

# View results
print(results_df)

# Optional: Adjust p-values for multiple comparisons using p.adjust
results_df$Adjusted_P_Value <- p.adjust(results_df$P_Value, method = "BH")

# View adjusted results
print(results_df) %>% format(., scientific = FALSE)
write_csv(results_df, "CHEMISTRY/anthraquinone_over-represented_results.csv")
