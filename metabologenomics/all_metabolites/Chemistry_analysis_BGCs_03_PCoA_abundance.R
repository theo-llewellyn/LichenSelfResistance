###################################################
# Code to analyse chemical data alongside BGCsc   #
# paired dataset of 83 genomes and metabolomes    #
# SCRIPT 3: PCoA                                  #
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

#Damiens MolNotator compounds
all_compounds <- read_csv("CHEMISTRY/OneDrive_2022-08-16/Kew Gardens/MolNotator_Damien/node_table.csv")
# FORMAT MolNotatori compounds
###Filter for computed neutrals, removes blanks, merge apothecia and thallus, edit accession typos and remove any empty rows
molecules_abundance <- all_compounds %>% 
  #remove blanks and only keep neutrals
  filter(status_universal == "neutral" & Blk2 == 0 & Blk4 == 0 & Blk5 == 0 & Blk6 == 0 & Blk7 == 0 & Blk8 == 0 & blk == 0 & blk3 == 0) %>% 
  dplyr::select(-c(1:24,125:145)) %>% 
  dplyr::select(sort(peek_vars())) %>% 
  #take mean of thallus and apothecia from same samples
  mutate(LIQ109XAAU = rowSums(across(c(`LIQ109XAAUA`, `LIQ109XAAUT`)))/2, LIQ146XSP = rowSums(across(c(`LIQ146XSPA`, `LIQ146XSPT`)))/2, LIQ164TEHY = rowSums(across(c(`LIQ164TEHYA`, `LIQ164TEHYT`)))/2, LIQ165TEEX = rowSums(across(c(`LIQ165TEEXA`, `LIQ165TEEXT`)))/2, LIQ166TFLAV = rowSums(across(c(`LIQ166TFLAVA`, `LIQ166TFLAVT`)))/2, LIQ240SEVIL = rowSums(across(c(`LIQ240SEVILA`, `LIQ240SEVILT`)))/2, LIQ69CCAR = rowSums(across(c(`LIQ69CCARA`, `LIQ69CCART`)))/2, LIQ72TVIL = rowSums(across(c(`LIQ72TVILA`, `LIQ72TVILT`)))/2, LIQ78THCR = rowSums(across(c(`LIQ78THCRA`, `LIQ78THCRT`)))/2, LIQ84UMVE = rowSums(across(c(`LIQ84UMVEA`, `LIQ84UMVET`)))/2) %>% 
  #remove unmerged samples
  dplyr::select(-c(2,3,6,7,11:16,66,67,79,80,83,84,89,90,95,96)) %>% 
  rename(LIQ106LELE = LIQ106LELEd, LIQ143CAAG =  LIQ143CAAT, LIQ147TEX = LIQ147TEXA, LIQ189CABR = LIQ189CASP, LIQ194TCHR = LIQ194THCR, LIQ246FLCI = LIQ246FLLId, `LIQ73XASTE-2` = LIQ73XASTE, `LIQ75XAME-2` = LIQ75XAME, `LIQ81XMZF-2`=LIQ81XZMF, LIQ92SELA_2 = LIQ92SELAT, `LIQ94LETR-2`=LIQ94LETR, LIQ109XAAU_2=LIQ109XAAU, `LIQ78TCHR-2`=LIQ78THCR) %>% 
  dplyr::select(-c(LIQ154CACI,LIQ181ARBA,LIQ206HESP,LIQ208TUSP,LIQ229XACA,LIQ241ATPY,LIQ247SOCH,LIQ244PYSP)) %>%
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))

rotated_molecules_abundance <- molecules_abundance %>%
  #normalise by total abundance for that sample
  mutate(across(everything(), ~ .x / sum(.x))) %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) %>%
  as.data.frame()

#set taxa as rownames and remove first column and any rows with a total of 0 as cant do dissimilarity on empty rows
rownames(rotated_molecules_abundance) <- rotated_molecules_abundance$name
rotated_molecules_abundance <- rotated_molecules_abundance[,-1]
rowSums(rotated_molecules_abundance)

#check same length
nrow(dissimilarity_matrix) == nrow(rotated_molecules_abundance)
#make new dissimilarity matrix removing empty rows in molecules
row.names.remove <- c("LIQ179PYOC","LIQ230CAOCH","LIQ231CALA")
dissimilarity_matrix_filtered_mol <- dissimilarity_matrix[!(row.names(dissimilarity_matrix) %in% row.names.remove),]
#check same length
nrow(dissimilarity_matrix_filtered_mol) == nrow(rotated_molecules_abundance)



##### TREE DATA
Leca82T_tree <- read.tree("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TREE_INFERENCE/IQTree_Leca82T_50p_tAL/concat_OF_tAl_renamed_rooted.treefile")
Leca82T_tree$tip.label <- gsub("'","",Leca82T_tree$tip.label)
Leca82T_tree$node.label <- NULL

##MOLECULES
#subset tree to contain just taxa in molecule matrix and order molecule matrix to be same as tree
abundance_ordered_mol <- match.phylo.data(phy = Leca82T_tree, data = rotated_molecules_abundance)
abundance_ordered_bgcf <- match.phylo.data(phy = Leca82T_tree, data = dissimilarity_matrix_filtered_mol)


#PcoA 
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


#make vector of colours
col_vec1 <- c("#5d6941","black","#d99627","#96adcb","#d1d7b3","pink","#a74207", "#666f6e")
bgcs_colours_orders <- taxon_colour_function(colours = col_vec1,
                                             taxa_list = rownames(abundance_ordered_bgcf$data),
                                             variable = "order")

col_vec2 <- c("#5d6941","black","#d99627","#96adcb")
bgcs_colours_subfamily <- taxon_colour_function(colours = col_vec2,
                                                taxa_list = rownames(abundance_ordered_bgcf$data),
                                                variable = "subfamily")

bgcs_colours_substrate <- taxon_colour_function(colours = col_vec2,
                                                taxa_list = rownames(abundance_ordered_bgcf$data),
                                                variable = "substrate")

bgcs_colours_growthform <- taxon_colour_function(colours = col_vec2,
                                                 taxa_list = rownames(abundance_ordered_bgcf$data),
                                                 variable = "growthform")

bgcs_colours_anthphenotype <- taxon_colour_function(colours = col_vec2,
                                                    taxa_list = rownames(abundance_ordered_bgcf$data),
                                                    variable = "anthphenotype")

library(RColorBrewer)
col_vec3 <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(metadata$genus)))
bgcs_colours_genus <- taxon_colour_function(colours = col_vec3,
                                            taxa_list = rownames(abundance_ordered_bgcf$data),
                                            variable = "genus")

#pcoa of bgc and molecules data
par(mfrow=c(1,2))
pco_bgc <- wcmdscale(vegdist(abundance_ordered_bgcf$data, method = "jaccard"), eig = TRUE)
plot(round(eigenvals(pco_bgc),3))
pco_mol <- wcmdscale(vegdist(abundance_ordered_mol$data, method = "bray"), eig = TRUE)
plot(round(eigenvals(pco_mol),3))

#just Teloschistales
pco_bgc_Telos <- wcmdscale(vegdist(abundance_ordered_bgcf$data[13:79,], method = "jaccard"), eig = TRUE)
pco_mol_Telos <- wcmdscale(vegdist(abundance_ordered_mol$data[13:79,], method = "bray"), eig = TRUE)

#fit environmental variables to ordination. Can then add these to ordiplots by just typing plot(fit_env_mol) after
fit_env_mol <- envfit(ord=pco_mol_Telos, env=metadata[c(14:47,50:82),2:6])
fit_env_bgc <- envfit(ord=pco_bgc_Telos, env=metadata[c(14:47,50:82),2:6])

#as we cant access loadings in PCoA due to information loss when converting to dissimilarity, we can fit the BGCs as environmental variables for the biplot to see which ones correlate with the dispersal of points from PCoA
efit <- envfit(pco_bgc_Telos, abundance_ordered_bgcf$data[13:79,])


png("FIGURES/BGC_LCMS_Leca82T_PCoA_labelled_abundance.png",  res = 300, width = 3500, height = 2000)
par(mfrow=c(1,2))
ordipointlabel(pco_bgc,
               display = "sites",
               scaling = "symmetric",
               col = bgcs_colours_orders$colour_vector,
               pch = 19,
               xlab = "PCoA axis 1",
               ylab = "PCoA axis 2",
               cex = 0.5,
               xlim=c(-.65,.4),
               ylim=c(-.4,.4))
abline(h = 0, v = 0, lty = 2)
ordiellipse(pco_bgc$points, groups = na.omit(bgcs_colours_orders$taxon_vector),
            kind = "ehull", col = col_vec1,
            scaling = "symmetric", lwd = 2, lty = 2)
legend("bottomleft",
       legend = unique(bgcs_colours_orders$taxon_vector),
       col = unique(bgcs_colours_orders$colour_vector),
       pch = 19)
ordipointlabel(pco_mol,
               display = "sites",
               scaling = "symmetric",
               col = bgcs_colours_orders$colour_vector,
               pch = 19,
               xlab = "PCoA axis 1",
               ylab = "PCoA axis 2",
               cex = 0.5,
               xlim=c(-.5,.5),
               ylim=c(-.4,.4))
abline(h = 0, v = 0, lty = 2)
ordiellipse(pco_mol$points, groups = na.omit(bgcs_colours_orders$taxon_vector),
            kind = "ehull", col = col_vec1,
            scaling = "symmetric", lwd = 2, lty = 2)
dev.off()


png("FIGURES/BGC_LCMS_Leca82T_Telos_PCoA_labelled_abundance.png",  res = 300, width = 3500, height = 2000)
par(mfrow=c(1,2))
ordipointlabel(pco_bgc_Telos,
               display = "sites",
               scaling = "symmetric",
               pch = 19,
               xlab = "PCoA axis 1",
               ylab = "PCoA axis 2",
               cex = 0.5,
               col = bgcs_colours_subfamily$colour_vector[13:79],
               xlim=c(-.4,.65))
abline(h = 0, v = 0, lty = 2)
ordiellipse(pco_bgc_Telos$points, groups = na.omit(bgcs_colours_subfamily$taxon_vector),
            kind = "ehull", col = col_vec2,
            scaling = "symmetric", lwd = 2, lty = 2)
legend("topright",
       legend = unique(bgcs_colours_subfamily$taxon_vector),
       col = unique(bgcs_colours_subfamily$colour_vector),
       pch = 19)
ordipointlabel(pco_mol_Telos,
               display = "sites",
               scaling = "symmetric",
               pch = 19,
               xlab = "PCoA axis 1",
               ylab = "PCoA axis 2",
               cex = 0.5,
               col = bgcs_colours_subfamily$colour_vector[13:79],
               xlim=c(-.55,.5))
ordiellipse(pco_mol_Telos$points, groups = na.omit(bgcs_colours_subfamily$taxon_vector),
            kind = "ehull", col = col_vec2,
            scaling = "symmetric", lwd = 2, lty = 2)
abline(h = 0, v = 0, lty = 2)
dev.off()
