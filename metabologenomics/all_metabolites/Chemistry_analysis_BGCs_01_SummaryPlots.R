###################################################
# Code to analyse chemical data alongside BGCsc   #
# paired dataset of 83 genomes and metabolomes    #
# SCRIPT 1: SUMMARY PLOTS                         #
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


#read in tree data
Leca82T_tree <- read.tree("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TREE_INFERENCE/IQTree_Leca82T_50p_tAL/concat_OF_tAl_renamed_rooted.treefile")
Leca82T_tree$tip.label <- gsub("'","",Leca82T_tree$tip.label)
Leca82T_tree$node.label <- NULL
Leca79T_tree <- drop.tip(Leca82T_tree, c("LIQ179PYOC","LIQ230CAOCH","LIQ231CALA"))
library(phylogram)
Leca79T_dendrogram <- as.dendrogram.phylo(Leca79T_tree)

#order compounds by how many samples found in
rotated_molecules_pres_abs_ordered <- rotated_molecules_pres_abs[,order(colSums(rotated_molecules_pres_abs,na.rm=TRUE), decreasing = TRUE)]
#order taxa by phylogeny
rotated_molecules_pres_abs_ordered <- rotated_molecules_pres_abs_ordered[rev(Leca79T_tree$tip.label),]
#remove singleton compounds
rotated_molecules_pres_abs_ordered <- rotated_molecules_pres_abs_ordered[,colSums(rotated_molecules_pres_abs_ordered) > 1]

png("FIGURES/Molecules_heatmap_hclust_ordered.png",  res = 300, width = 2000, height = 2000)
#heatmap clustered by presence absence patterns of metabolites
molecules_heatmap_clustered <- heatmap.2(data.matrix(rotated_molecules_pres_abs_ordered), trace = "none",
          cexRow = .4,
          scale = "none",
          xlab = "Molnotator computed neutrals",
          dendrogram='column', Rowv=TRUE, Colv=TRUE,
          labCol = FALSE,
          col = c("grey30","orange"),
          #margins = c(2,3),
          offsetRow = -0.1,
          #Rowv = NULL,
          key = FALSE,
          lhei=c(.3,2), lwid = c(0.1,2))
dev.off()

hclust_dendrogram <- molecules_heatmap_clustered$rowDendrogram

png("FIGURES/Molecules_heatmap_phylogeny_ordered.png",  res = 300, width = 2000, height = 2000)
#heatmap ordered by phylogeny
heatmap.2(data.matrix(rotated_molecules_pres_abs_ordered), trace = "none",
          scale = "none",
          cexRow = .4,
          margins = c(2,3),
          offsetRow = -0.1,
          xlab = "Molnotator computed neutrals",
          labCol = FALSE,
          col = c("grey30","orange"),
          Rowv = NULL,
          key = FALSE,
          lhei=c(.3,2), lwid = c(0.1,2))
dev.off()

#compare hclust of taxa due to metabolite presence absence with phylogeny
library(dendextend)
png("FIGURES/Leca82Tphylogeny_Molecules_hclustdendrogram_tanglegram.png",  res = 300, width = 2000, height = 2000)
tanglegram(hclust_dendrogram, as.cladogram(Leca79T_dendrogram), common_subtrees_color_branches = TRUE, common_subtrees_color_lines = TRUE, highlight_branches_lwd = FALSE, axes = FALSE, main_left = "hclust", main_right = "phylogeny", lwd = 1, lab.cex = 0.5)
dev.off()

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
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) -> positive_gnps_full

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

#make a df with the component indexes and the compound names and removes -1 component indexes
pos_table1 <- positive_gnps[,c(34,35)] %>% filter(componentindex != -1)
#make dataframe of column names of datamatrix
pos_table2 <- tibble(componentindex = as.numeric(colnames(positive_gnps_pres_abs)))
#match dataframe to see which component indexes have matched compounds
pos_table3 <- full_join(x = pos_table1, y = pos_table2) %>% 
  arrange(componentindex) %>%
  distinct(.keep_all = TRUE) %>%
  drop_na() 

#these are the positions in the column names that have matched compounds
pos_col_ids <- which(colnames(positive_gnps_pres_abs) %in% unique(pos_table3$componentindex))
col_vec1 <- c("#5d6941","black","#d99627","#96adcb","#d1d7b3","pink","#a74207", "#666f6e","orchid")

#all anotated columns black
colCols_pos <- rep("white",length(colnames(positive_gnps_pres_abs)))
for(i in 1:length(colCols_pos)){
  if(is.element(i, pos_col_ids)){
    colCols_pos[i] <- "black"
  }
}

#each compound a different colour
count <- 0
for(i in pos_col_ids){
  count <- count + 1
  colCols_pos[i] <- col_vec1[count]
}

positive_gnps_pres_abs_ordered <- positive_gnps_pres_abs[,order(colSums(positive_gnps_pres_abs,na.rm=TRUE), decreasing = TRUE)]

png("FIGURES/GNPS_pos_heatmap_1.png", res = 300, width = 3500, height = 3500)
heatmap.2(data.matrix(positive_gnps_pres_abs), trace = "none",
          scale = "none",
          ColSideColors=colCols_pos,
          xlab = "+ve ion GNPS clusters",
          labCol = FALSE,
          col = c("grey30","orange"),
          keysize = .7,
          density.info="none",
          key.title = "",
          key.xlab = NA,
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("absent", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("present", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
legend("left", title = "compounds",legend=c("physcion","Variolaric acid","BIS(2-ETHYLHEXYL)PHTHALATE","Pheophytin a","Caloploicin","emodin","Pheophorbide A","Secalonic acid","Usnic acid"), 
       fill=c("#5d6941","black","#d99627","#96adcb","#d1d7b3","pink","#a74207", "#666f6e","orchid"), cex=0.8, box.lty=0)
dev.off()




#make table of just GNPS clusters that have annotated compounds
positive_gnps_annotated <- positive_gnps[!is.na(positive_gnps$Compound_Name),] %>% 
  dplyr::select(c(34,35,53,91:190)) %>%
  filter(componentindex != -1 & BLANK == 0) %>% 
  dplyr::select(-c(1,3)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  #merge columns for apothecia and thallus of same sample
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))) %>%
  #remove thallus and apothecia columns to leave a single column per accession
  dplyr::select(-c(3,4,7,8,12:17,67,68,80,81,84,85,90,91,96,97)) %>%
  #correct typos
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  #remove taxa not in BGC dataset
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) %>%
  #merge subclusters so there is one row per subcluster rather than per node of network (ion)
  aggregate(. ~  Compound_Name, data = ., sum) %>% 
  tibble() %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

positive_gnps_annotated$anthraquinone <- c("Y","Y","Y","Y","N","N","N","N","N","N")

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
#make new dissimilarity matrix removing empty rows in molecules
row.names.remove <- c("LIQ230CAOCH","LIQ231CALA")
dissimilarity_matrix_filtered_gnps_neg <- dissimilarity_matrix[!(row.names(dissimilarity_matrix) %in% row.names.remove),]
#check same length
nrow(dissimilarity_matrix_filtered_gnps_neg) == nrow(negative_gnps_pres_abs)

#make a df with the component indexes and the compound names and removes -1 component indexes
neg_table1 <- negative_gnps[,c(34,35)] %>% filter(componentindex != -1)
#make dataframe of column names of datamatrix
neg_table2 <- tibble(componentindex = as.numeric(colnames(negative_gnps_pres_abs)))
#match dataframe to see which component indexes have matched compounds
neg_table3 <- full_join(x = neg_table1, y = neg_table2) %>% 
  arrange(componentindex) %>%
  distinct(.keep_all = TRUE) %>%
  drop_na() 

#these are the positions in the column names that have matched compounds
neg_col_ids <- which(colnames(negative_gnps_pres_abs) %in% unique(neg_table3$componentindex))
#for loop that makes a list of colours for the component ids, if the id is one of the matched compounds we will change the colour to orange
colCols_neg <- rep("white",length(colnames(negative_gnps_pres_abs)))
for(i in 1:length(colCols_neg)){
  if(is.element(i, neg_col_ids)){
    colCols_neg[i] <- "black"
  }
}
#check right length
length(colCols_neg) == length(colnames(negative_gnps_pres_abs))

#heatmap
png("FIGURES/GNPS_neg_heatmap.png", res = 300, width = 3500, height = 3500)
heatmap.2(data.matrix(negative_gnps_pres_abs),
          trace = "none",
          scale = "none",
          ColSideColors=colCols_neg,
          xlab = "-ve ion GNPS clusters",
          labCol = FALSE,
          col = c("grey30","orange"),
          keysize = 0.75,
          density.info="none",
          key.title = "",
          key.xlab = NA,
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("absent", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("present", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
dev.off()


#some accession have no gnps clusters so should be removed from dataset
# LIQ179PYOC LIQ230CAOCH LIQ231CALA
negative_gnps_pres_abs <- negative_gnps_pres_abs[-c(9,44,45),]

#make table of just GNPS clusters that have annotated compounds
negative_gnps_annotated <- negative_gnps[!is.na(negative_gnps$Compound_Name),] %>% 
  dplyr::select(c(34,35,53,91:190)) %>%
  filter(componentindex != -1 & BLANK == 0) %>% 
  dplyr::select(-c(1,3)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))) %>%
  dplyr::select(-c(3,4,7,8,12:17,67,68,80,81,84,85,90,91,96,97)) %>%
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) %>% 
  aggregate(. ~  Compound_Name, data = ., sum) %>% 
  tibble() %>% 
  mutate_if(is.numeric, ~1 * (. > 0))
negative_gnps_annotated$anthraquinone <- c("Y","N","N","Y","N","N","Y","N","N","N","N","N","N")

#heatmap
both_gnps <- full_join(x = positive_gnps_annotated, y = negative_gnps_annotated) %>% distinct(Compound_Name, .keep_all = TRUE)

both_gnps$Compound_Name <- gsub(pattern = "2-chloro.*anthracene-9,10-dione",
                                replacement = "7-Chlorocitreosein",
                                both_gnps$Compound_Name) %>%
  gsub(pattern = "Massbank:.*[0-9] ",
       replacement = "")

both_gnps <- column_to_rownames(both_gnps, 'Compound_Name')

#merge the two usnic acid rows as they are the same compound just different spellings as from two databases
merged_row <- both_gnps[8,-83] + both_gnps[9,-83]
merged_row[merged_row==2] <- 1
both_gnps <- add_row(both_gnps,merged_row)
#remove the anthraquinone variable and the old usnic acid rows
both_gnps <- both_gnps[-c(8,9),-83]
rownames(both_gnps)[rownames(both_gnps) == "...22"] <- "Usnic acid"

png("FIGURES/GNPS_Compounds_presabs_heatmap.png", width = 210, height = 210, res = 300, units = "mm")
heatmap.2(t(data.matrix(both_gnps)), trace = "none",
          scale = "none",
          keysize = 1,
          margins = c(10,10),
          col = c("grey30","orange"),
          density.info="none",
          key.title = "",
          key.xlab = NA,
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("absent", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("present", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
dev.off()

#PCOA of anthraquinones
rotated_anthraquinones <- both_gnps[c(2:4,9,11,14),] %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>%
  column_to_rownames('name') %>%
  as.data.frame()

rotated_anthraquinones <- rotated_anthraquinones[rowSums(rotated_anthraquinones[])>0,]

#read in taxonomy info
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
anth_colours_orders <- taxon_colour_function(colours = col_vec1,
                                             taxa_list = rownames(rotated_anthraquinones),
                                             variable = "order")

col_vec2 <- c("#5d6941","black","#d99627","#96adcb")
anth_colours_subfamily <- taxon_colour_function(colours = col_vec2,
                                                taxa_list = rownames(rotated_anthraquinones),
                                                variable = "subfamily")

#plot pcoa
pco_anthraquinones <- wcmdscale(vegdist(rotated_anthraquinones[c(1:7,9:21,24:26,31:58,60:75),], method = "jaccard"), eig = TRUE)
ordipointlabel(pco_anthraquinones,
               display = "sites",
               scaling = "symmetric",
               col = anth_colours_subfamily$colour_vector[c(1:7,9:21,24:26,31:58,60:75)],
               pch = 19,
               xlab = "PCoA axis 1",
               ylab = "PCoA axis 2",
               cex = 0.5,
               xlim=c(-.6,.6))
ordiellipse(pco_anthraquinones$points, groups = na.omit(anth_colours_subfamily$taxon_vector),
            kind = "ehull", col = col_vec2,
            scaling = "symmetric", lwd = 2, lty = 2)
legend("topright",
       legend = unique(anth_colours_orders$taxon_vector),
       col = unique(anth_colours_orders$colour_vector),
       pch = 19)

#summary histograms to see how many taxa each metabolite was found in
molecules_barchart <- ggplot(as.data.frame(colSums(rotated_molecules_pres_abs)), aes(`colSums(rotated_molecules_pres_abs)`)) +
  geom_bar() +
  xlab("# samples") +
  ylab("Molecule counts")
posions_barchart <- ggplot(as.data.frame(colSums(rotated_positive_ions_pres_abs_filtered)), aes(`colSums(rotated_positive_ions_pres_abs_filtered)`)) +
  geom_bar() +
  xlab("# samples") +
  ylab("+ve ion counts")
negions_barchart <- ggplot(as.data.frame(colSums(rotated_negative_ions_pres_abs_filtered)), aes(`colSums(rotated_negative_ions_pres_abs_filtered)`)) +
  geom_bar() +
  xlab("# samples") +
  ylab("-ve ion counts")
posgnps_barchart <- ggplot(as.data.frame(colSums(positive_gnps_pres_abs)), aes(`colSums(positive_gnps_pres_abs)`)) +
  geom_bar() +
  xlab("# samples") +
  ylab("+ve GNPS families counts")
neggnps_barchart <- ggplot(as.data.frame(colSums(negative_gnps_pres_abs)), aes(`colSums(negative_gnps_pres_abs)`)) +
  geom_bar() +
  xlab("# samples") +
  ylab("-ve GNPS families counts")

png("FIGURES/Summary_barcharts_metabolomes.png",height = 1654, width = 2339)
plot_grid(molecules_barchart,posions_barchart,negions_barchart,posgnps_barchart,neggnps_barchart)
dev.off()


#plot of molecule number vs gcf number
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

#plot the bgc,molecule etc values against the phylogenetic tree
library(ggnewscale)
library(scales)
bgc_mol_heatmap <- ggplot(merged_tables, aes(y = accession)) +
  geom_tile(aes(x=1, fill = bgcs)) +
  scale_fill_gradient(high = "black", low = "white") +
  #scale_fill_viridis(option = "A") +
  new_scale_fill() +
  geom_tile(aes(x = 2, fill = molecules)) +
  scale_fill_gradient(high = "darkorange", low = "white") +
  #scale_fill_viridis(option = "D") +
  new_scale_fill() +
  geom_tile(aes(x = 3, fill=pos_GNPS)) +
  scale_fill_gradient(high = "darkgreen", low = "white") +
  #scale_fill_viridis(option = "B") +
  new_scale_fill() +
  geom_tile(aes(x = 4, fill=neg_GNPS)) +
  scale_fill_gradient(high = "darkblue", low = "white") +
  #scale_fill_viridis(option = "C") + 
  coord_equal() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        axis.title = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90))

png("FIGURES/number_bgcs_mols_Leca83T.png",  res = 300, width = 2480, height = 2480)
bgc_mol_heatmap %>% insert_left(g, width = 10)
dev.off()


png("FIGURES/No_moelcules_per_sample_Leca83T.png",  res = 300, width = 3500, height = 2000)
plot0 <- ggplot(merged_tables, aes(x = molecules, y = accession)) +
  geom_bar(stat="identity")
dev.off()

plot1 <- ggplot(merged_tables, aes(x=bgcs, y = molecules)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
plot2 <- ggplot(merged_tables, aes(x=bgcs, y = pos_ions)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
plot3 <- ggplot(merged_tables, aes(x=bgcs, y = neg_ions)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
plot4 <- ggplot(merged_tables, aes(x=bgcs, y = pos_GNPS)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
plot5 <- ggplot(merged_tables, aes(x=bgcs, y = neg_GNPS)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)


png("FIGURES/BGC_vs_moelcules_pos_neg_Leca83T.png",  res = 300, width = 2000, height = 3500)
plot_grid(plot0, plot1, plot2, plot3, plot4, plot5, ncol = 2, labels = c("","molecules","+ve ions","-ve ions","+ve GNPS","-ve GNPS"))
dev.off()

lm_BGC_mol <- lm(molecules ~ bgcs, merged_tables)
summary(lm_BGC_mol)

cor(merged_tables$bgcs, merged_tables$molecules)
cor(merged_tables$bgcs, merged_tables$pos_ions)
cor(merged_tables$bgcs, merged_tables$neg_ions)
cor(merged_tables$bgcs, merged_tables$pos_GNPS)
cor(merged_tables$bgcs, merged_tables$neg_GNPS)
#0.3411398, 0.2809938, 0.4083145, 0.3091673, 0.3819285


### ALTERNATIVE GNPS ANNOTATED COMPOUNDS HEATMAP
positive_gnps_annotated <- positive_gnps[!is.na(positive_gnps$Compound_Name),] %>% 
  dplyr::select(c(34,35,53,91:190)) %>%
  filter(BLANK == 0) %>% 
  dplyr::select(-c(1,3)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  #merge columns for apothecia and thallus of same sample
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))) %>%
  #remove thallus and apothecia columns to leave a single column per accession
  dplyr::select(-c(3,4,7,8,12:17,67,68,80,81,84,85,90,91,96,97)) %>%
  #correct typos
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  #remove taxa not in BGC dataset
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) %>%
  #merge subclusters so there is one row per subcluster rather than per node of network (ion)
  aggregate(. ~  Compound_Name, data = ., sum) %>% 
  tibble() %>% 
  mutate_if(is.numeric, ~1 * (. > 0))


negative_gnps_annotated <- negative_gnps[!is.na(negative_gnps$Compound_Name),] %>% 
  dplyr::select(c(34,35,53,91:190)) %>%
  filter(BLANK == 0) %>% 
  dplyr::select(-c(1,3)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`))), LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`))), LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`))), LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`))), LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`))), LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`))), LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`))), LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`))), LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`))), LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))) %>%
  dplyr::select(-c(3,4,7,8,12:17,67,68,80,81,84,85,90,91,96,97)) %>%
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG = LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) %>% 
  aggregate(. ~  Compound_Name, data = ., sum) %>% 
  tibble() %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

negative_gnps_annotated$Compound_Name <- gsub(pattern = "Caloploicin",
                                              replacement = "Caloploicin1",
                                              negative_gnps_annotated$Compound_Name) %>%
  gsub(pattern = "Variolaric acid",
       replacement = "Variolaric acid1")

#heatmap
both_gnps <- full_join(x = positive_gnps_annotated, y = negative_gnps_annotated)

#remove Asperphenamate_130119
both_gnps <- both_gnps[-1,]

both_gnps$Compound_Name <- gsub(pattern = "2-chloro.*anthracene-9,10-dione",
                                replacement = "7-Chlorocitreosein",
                                both_gnps$Compound_Name) %>%
  gsub(pattern = "Massbank:.*[0-9] ",
       replacement = "")

both_gnps <- column_to_rownames(both_gnps, 'Compound_Name')

#merge the two usnic acid rows as they are the same compound just different spellings as from two databases
merged_row <- c(both_gnps[9,] + both_gnps[10,])
merged_row[merged_row==2] <- 1
merged_row1 <- c(both_gnps[1,] + both_gnps[14,])
merged_row1[merged_row1==2] <- 1
merged_row2 <- c(both_gnps[3,] + both_gnps[19,])
merged_row2[merged_row2==2] <- 1
merged_row3 <- c(both_gnps[11,] + both_gnps[27,])
merged_row3[merged_row3==2] <- 1

both_gnps <- rbind(both_gnps,merged_row,merged_row1,merged_row2,merged_row3)
#remove the unmerged variables and rename merged
both_gnps <- both_gnps[-c(9,10,1,14,3,19,11,27),]
rownames(both_gnps)[rownames(both_gnps) == "1"] <- "Usnic acid"
rownames(both_gnps)[rownames(both_gnps) == "11"] <- "Caloploicin"
rownames(both_gnps)[rownames(both_gnps) == "12"] <- "Isokaempferide"
rownames(both_gnps)[rownames(both_gnps) == "13"] <- "Variolaric acid"
rownames(both_gnps)[rownames(both_gnps) == "MEHP|2-(((2-Ethylhexyl)oxy)carbonyl)benzoic acid|2-(2-ethylhexoxycarbonyl)benzoic acid"] <- "2EOCBA"

rotated_gnps_annotated <- t(data.matrix(both_gnps))
rotated_gnps_annotated <- rotated_gnps_annotated[rev(Leca79T_tree$tip.label),]

png("FIGURES/GNPS_annotated_heatmap_phylogeny_ordered.png",  res = 300, width = 2000, height = 2000)
heatmap.2(rotated_gnps_annotated,
          trace = "none",
          scale = "none",
          cexRow = .4,
          Rowv = NULL,
          Colv=TRUE,
          dendrogram = "none",
          margins = c(10,3),
          offsetRow = -0.1,
          offsetCol = -0.1,
          col = c("grey30","orange"),
          key = FALSE,
          lhei=c(.3,4), lwid = c(0.1,2))
dev.off()


rotated_gnps_annotated <- rotated_gnps_annotated[rowSums(rotated_gnps_annotated[])>0,]
col_vec1 <- c("#5d6941","black","#d99627","#96adcb","#d1d7b3","pink","#a74207", "#666f6e")
gnps_colours_orders <- taxon_colour_function(colours = col_vec1,
                                             taxa_list = rownames(rotated_gnps_annotated),
                                             variable = "order")

pco_gnps_annotated <- wcmdscale(vegdist(rotated_gnps_annotated, method = "jaccard"), eig = TRUE)
ordipointlabel(pco_gnps_annotated,
               display = "sites",
               scaling = "symmetric",
               col = gnps_colours_orders$colour_vector,
               pch = 19,
               xlab = "PCoA axis 1",
               ylab = "PCoA axis 2",
               cex = 0.5,
               xlim=c(-.6,.6))
