#################################################################
# Code to do hOUwie analysis and save the results               #
# Author: Theo Llewellyn                                        #
#################################################################

setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/")

library(corHMM)
library(data.table)
library(expm)
library(parallel)
library(tibble)
library(tidyverse)
library(tidyselect)
library(dplyr)
library(phytools)
library(OUwie)
library(ggrepel)
library(vegan)

#read in anth BGC gene presence absence data (tip states)
anth_genes <- read_csv("CLINKER/anth_BGCs_pres_abs_plus_outgroups.csv")
anth_genes <- anth_genes %>%
  mutate_at(vars(!taxon), factor)


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

####
anth_richness <- tibble(rownames(merged_df), diversity(merged_df, index = "shannon"))
anth_pres_abs <- merged_df %>%
  rowwise() %>%
  mutate(non_zero_count = sum(c_across(everything()) != 0)) %>%
  ungroup() %>%
  select(non_zero_count)
anth_richness <- cbind(anth_richness,anth_pres_abs)
colnames(anth_richness) <- c('taxon','anthraquinone_richness_shannon','anthraquinone_richness_pres_abs')


###
# Hypothesis 1: anthraquinone number evolution varies between Telos and non-Telos 
###

#First we test to see if anth num evolution differs in Telos and non-Telos
##### TREE DATA
#read in tree data
Leca82T_tree <- read.tree("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TREE_INFERENCE/IQTree_Leca82T_50p_tAL/concat_OF_tAl_renamed_rooted.treefile")
Leca82T_tree$tip.label <- gsub("'","",Leca82T_tree$tip.label)
Leca82T_tree$node.label <- NULL
#Plot the tree and the internal nodes to highlight the selective regimes:
nod.lab<-(Ntip(Leca82T_tree)+1):(Ntip(Leca82T_tree)+Leca82T_tree$Nnode)
Leca82T_tree$node.label<-nod.lab
plot(Leca82T_tree,show.node.label=1)
#label the Teloschistales nodes as 1 and outgroups as 2
simmap.tr<-paintSubTree(Leca82T_tree,node=83,state="1")
simmap.tr<-paintSubTree(simmap.tr,node=99,state="2")
plot(simmap.tr)

#make a dataset with the tip label, the clade number (1 or 2) and the number of anthraquinones
anth_richness$Reg <- c(rep("2",7),rep("1",3),rep("2",7),'1',rep("2",5),rep("1",2),rep("2",3),rep("1",6),rep("2",30),"1",rep("2",5),rep("1",2),rep("2",9),"1")

anth_richness$anth_num <- anth_richness$anthraquinone_richness_pres_abs/10
#reorder so matches format for OUwie
anth_richness_ouwie <- anth_richness[,c(1,4,5)]

#fit the BM (ignores the mapped regimes)
fitBM <- OUwie(simmap.tr,anth_richness_ouwie,model="BM1", simmap.tree=TRUE,root.station=TRUE)
#fit the BMS (allow each of the mapped regimes to evolve with a different value of the Brownian rate parameter)
fitBMS <- OUwie(simmap.tr,anth_richness_ouwie,model="BMS", simmap.tree=TRUE,root.station=FALSE)
#fit the OUM where the stable phenotypic optimum value (theta) is allowed to vary between clades
# OU associated with natural selection towards a specific value
fitOUM <- OUwie(simmap.tr,anth_richness_ouwie,model="OUM", simmap.tree=TRUE,root.station=FALSE)

#compare fit
aic<-setNames(c(fitBM$AIC,fitBMS$AIC,fitOUM$AIC), c("BM1","BMS","OUM"))
aic
aic.w(aic)


#same but for Shannon diversity
anth_richness_ouwie_shannon <- anth_richness[,c(1,4,2)]

#fit the BM (ignores the mapped regimes)
fitBM_shannon <- OUwie(simmap.tr,anth_richness_ouwie_shannon,model="BM1", simmap.tree=TRUE,root.station=TRUE)
#fit the BMS (allow each of the mapped regimes to evolve with a different value of the Brownian rate parameter)
fitBMS_shannon <- OUwie(simmap.tr,anth_richness_ouwie_shannon,model="BMS", simmap.tree=TRUE,root.station=FALSE)
#fit the OUM where the stable phenotypic optimum value (theta) is allowed to vary between clades
# OU associated with natural selection towards a specific value
fitOUM_shannon <- OUwie(simmap.tr,anth_richness_ouwie_shannon,model="OUM", simmap.tree=TRUE,root.station=FALSE)

#compare fit
aic_shannon<-setNames(c(fitBM_shannon$AIC,fitBMS_shannon$AIC,fitOUM_shannon$AIC), c("BM1","BMS","OUM"))
aic_shannon
aic.w(aic_shannon)


###
# Hypothesis 2: anthraquinone number evolution varies depending on presence-absence of genes in BGC
###

#Now for each gene in the anth BGC see if correlates with changes in anth number
genes <- colnames(anth_genes[,2:11])
aic_results <- c()
model_fit_results <- c()
for(gene in genes){
  print(gene)
  gene_simmap_tree <- read.simmap(paste("ASR_out/anth_jointASR_out_",gene,"_Leca82T_simmap.tre", sep = ""), format = "phylip")
  gene_data <- read_csv(paste("ASR_out/anth_jointASR_out_",gene,"_Leca82T_simmap.csv", sep = ""))
  gene_howie_data <- merge(gene_data, anth_richness_ouwie[,c(1,3)], by = "taxon")
  
  #fit the BM (ignores the mapped regimes)
  fitBM <- OUwie(gene_simmap_tree,gene_howie_data,model="BM1", simmap.tree=TRUE,root.station=TRUE)
  #fit the BMS (allow each of the mapped regimes to evolve with a different value of the Brownian rate parameter)
  fitBMS <- OUwie(gene_simmap_tree,gene_howie_data,model="BMS", simmap.tree=TRUE,root.station=FALSE)
  #fit the OUM where the stable phenotypic optimum value (theta) is allowed to vary between clades
  # OU associated with natural selection towards a specific value
  fitOUM <- OUwie(gene_simmap_tree,gene_howie_data,model="OUM", simmap.tree=TRUE,root.station=FALSE)
  
  #compare fit
  aic<-setNames(c(fitBM$AIC,fitBMS$AIC,fitOUM$AIC), c("BM1","BMS","OUM"))
  #highest weighted AIC is best
  aicw <- aic.w(aic)
  aic <- setNames(c(gene, aicw),c("gene","BM1","BMS","OUM"))
  aic_results <- rbind(aic_results, aic)
  
  #houwie
  print(paste("Starting hOUwie for ",gene," with BMV model ",sep = ""))
  pp_bmv <- hOUwie(gene_simmap_tree, gene_howie_data, rate.cat = 1, discrete_model = "ARD", 
                   continuous_model = "BMV")
  print(paste("Starting hOUwie for ",gene," with OU1 model ",sep = ""))
  pp_ou1 <- hOUwie(gene_simmap_tree, gene_howie_data, rate.cat = 1, discrete_model = "ARD", 
                   continuous_model = "OU1")
  print(paste("Starting hOUwie for ",gene," with BM1 model ",sep = ""))
  pp_bm1 <- hOUwie(gene_simmap_tree, gene_howie_data, rate.cat = 1, discrete_model = "ARD", 
                   continuous_model = "BM1")
  print(paste("Starting hOUwie for ",gene," with OUM model ",sep = ""))
  pp_oum <- hOUwie(gene_simmap_tree, gene_howie_data, rate.cat = 1, discrete_model = "ARD",
                   continuous_model = "OUM")
  print(paste("Starting hOUwie for ",gene," with BMV cid model ",sep = ""))
  pp_bmv_cid <- hOUwie(gene_simmap_tree, gene_howie_data, rate.cat = 2, discrete_model = "ARD", null.model = TRUE,
                       continuous_model = "BMV")
  print(paste("Starting hOUwie for ",gene," with OUM cid model ",sep = ""))
  pp_oum_cid <- hOUwie(gene_simmap_tree, gene_howie_data, rate.cat = 2, discrete_model = "ARD", null.model = TRUE,
                       continuous_model = "OUM")
  #compare models
  print(paste("Comparing models",gene,sep = ""))
  model_set <- list(bm1_fit = pp_bm1, ou1_fit = pp_ou1, bmv_fit = pp_bmv, oum_fit = pp_oum,
                    bmv_cid_fit = pp_bmv_cid, oum_cid_fit = pp_oum_cid)
  model_fit <- cbind(getModelTable(model_set, type = "BIC"),gene = rep(gene,6)) %>% rownames_to_column("Model")
  model_fit_results <- rbind(model_fit_results,model_fit)
  print(paste(gene," Done!",sep = ""))
}


write_csv(as.data.frame(aic_results), "OUWIE/OUwie_outgroup_AIC_results_2025.csv")
write_csv(as.data.frame(model_fit_results), "OUWIE/hOUwie_outgroup_model_fit_results_100simmaps_2025.csv")
