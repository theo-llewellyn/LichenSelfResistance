###################################################
# Code to perform ancestral state reconstruction  #
# of anthraquinone BGCs for 83 genomes            #
# SCRIPT 7: ASR						              #
# Author: Theo Llewellyn                          #
###################################################

#This differs to original analysis as it includes updated codificaiton of anth BGCs in outgroups

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

##### TREE DATA
Leca82T_tree <- read.tree("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TREE_INFERENCE/IQTree_Leca82T_50p_tAL/concat_OF_tAl_renamed_rooted.treefile")
Leca82T_tree$tip.label <- gsub("'","",Leca82T_tree$tip.label)
Leca82T_tree$node.label <- NULL
#remove taxa without anth. BGC at all
Leca82T_tree <- drop.tip(Leca82T_tree, c("LIQ143CAAG","LIQ179PYOC","LIQ187CASP","LIQ192CASP","LIQ200XMHS","LIQ218CACI","LIQ231CALA","LIQ235IKAU","LIQ242VAFU","LIQ243FLLI","LIQ249FLNA","LIQ230CAOCH"))


#Ancestral state reconstruction of BGCF presence absence
library(corHMM)
#read in anthraquinone gene trait matrix
anth_genes <- read_csv("CLINKER/anth_BGCs_pres_abs_plus_outgroups.csv")

anth_genes <- anth_genes %>% 
  #select(-BGCF) %>% 
  filter(Organism != "LIQ230CAOCH") %>% 
  mutate_at(c('EthD','PKSI','MBL','ABC', 'MethylT', 'GammaG', "AAperm", "Oreduct", 'KinaseP','DUF'), funs(str_replace(., "0", "2_absent"))) %>%
  mutate_at(c('EthD','PKSI','MBL','ABC', 'MethylT', 'GammaG', "AAperm", "Oreduct", 'KinaseP','DUF'), funs(str_replace(., "1", "0_clustered"))) %>%
  mutate_at(c('EthD','PKSI','MBL','ABC', 'MethylT', 'GammaG', "AAperm", "Oreduct", 'KinaseP','DUF'), funs(str_replace(., "2", "1_unclustered"))) %>%
  mutate_at(vars(!Organism), factor)

#prepare dataframe with all the characters
anth_trait_allgenes <- anth_genes %>% column_to_rownames('Organism')
colors_all<-setNames(replicate(ncol(anth_trait_allgenes),setNames(c("white","black","grey"),c("2_absent","0_clustered","1_unclustered")),
                               simplify=FALSE),colnames(anth_trait_allgenes))
plotTree.datamatrix(Leca82T_tree, anth_trait_allgenes,
                    fsize = 0.45, colors = colors_all)


#choose model
aic_results <- c()
for(gene in colnames(anth_trait_allgenes)){
  anth_trait <- anth_trait_allgenes[[gene]]
  names(anth_trait) <- rownames(anth_trait_allgenes)
  fitER <- fitMk(Leca82T_tree, anth_trait, model = "ER")
  fitARD <- fitMk(Leca82T_tree, anth_trait, model = "ARD")
  aic <- c(gene, AIC(fitER),AIC(fitARD))
  aic <- setNames(c(gene,AIC(fitER),AIC(fitARD)),c("gene","ER","ARD"))
  aic_results <- rbind(aic_results, aic)
}


#for loop which goes through each gene and does joint and marginal likelihood ASR then plots on tree
genes <- colnames(anth_genes[,2:11])
for(gene in genes){
  pdf(paste("ASR_out/anth_ASR_out_",gene,"Leca82T_new.pdf", sep = ""))
  print(gene)
  #make dataframe
  anth_trait <- setNames(as.factor(anth_genes[[gene]]), anth_genes$Organism)
  #make a list with the colours to assign to each trait value
  cols <- setNames(c("black","grey","white"), levels(anth_trait))
  #need to convert to dataframe for corHMM
  anth_gene_data <- data.frame(taxon = names(anth_trait), gene = as.numeric(anth_trait)-1)
  #do ASR_out
  fit_marginal <- corHMM(Leca82T_tree, anth_gene_data, node.states = "marginal",
                         rate.cat = 1, model = "ARD")
  fit_joint <- corHMM(Leca82T_tree, anth_gene_data, node.states = "joint",
                      rate.cat = 1, model = "ARD")
  #save joint states as simmap tree with mapped internal node values and then save to file
  simmap <- makeSimmap(Leca82T_tree, anth_gene_data, fit_joint$solution, rate.cat = 1)
  write.simmap(simmap[[1]], paste("ASR_out/anth_jointASR_out_",gene,"_Leca82T_simmap.tre", sep = ""))
  # save joint states table
  write_csv(fit_joint$data, paste("ASR_out/anth_jointASR_out_",gene,"_Leca82T_simmap.csv", sep = ""))
  #plot marginal
  plotTree.datamatrix(Leca82T_tree, as.data.frame(anth_trait),
                      header = FALSE, fsize = 0.45, colors = cols)
  nodelabels(pie = fit_marginal$states, piecol = cols, cex = 0.25)
  legend("topleft", legend = levels(anth_trait), pt.cex = 1.5, pt.bg = cols, bty = "n", cex = 0.8, pch = 22)
  title(gene, line = -1)
  #plot joint
  plotTree.datamatrix(Leca82T_tree, as.data.frame(anth_trait),
                      header = FALSE, fsize = 0.45, colors = cols)
  nodelabels(pie = to.matrix(levels(anth_trait)[fit_marginal$phy$node.label], levels(anth_trait)), piecol = cols, cex = 0.25)
  legend("topleft", legend = levels(anth_trait), pt.cex = 1.5, pt.bg = cols, bty = "n", cex = 0.8, pch = 22)
  title(gene, line = -1)
  dev.off()
}
