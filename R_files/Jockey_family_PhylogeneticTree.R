#################
#################
#### Phylogenetic tree of Consensus G2/Jockey-3 and other Jockey family elements
#################
#################


####
# load libraries
####

library(ape)
library(phytools)
library(geiger)
library(evobiR)
options(scipen=999)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtree)


#### load tree

jock.tree <- ape::read.tree(file="/Users/lucashemmer/Documents/Cecile_simulans/RAxML_bipartitions.alignment_jockey_family_Final.automre")


#### root on true outgroup
jock.tree1 <- root.phylo(jock.tree, "Jockey-1_Dere", resolve.root = TRUE)

#### drop the outgroup
jock.tree1 <- drop.tip(jock.tree1,"Jockey-1_Dmel")

#### drop redudant fragment
jock.tree1 <- drop.tip(jock.tree1,"Jockey-3_Dyak_fragment")


#### assigning labels
jock.tree1$tip.label[jock.tree1$tip.label=="Jockey-3_Dsec_Clade2"] <- "Jockey-3_Dsec_B"
jock.tree1$tip.label[jock.tree1$tip.label=="Jockey-3_Dsim_Clade2"] <- "Jockey-3_Dsim_B"
jock.tree1$tip.label[jock.tree1$tip.label=="Jockey-3_Dmau_Clade2"] <- "Jockey-3_Dmau_B"
jock.tree1$tip.label[jock.tree1$tip.label=="Jockey-3_Dsec_mel-like"] <- "Jockey-3_Dsec_A"
jock.tree1$tip.label[jock.tree1$tip.label=="Jockey-3_Dsim_mel-like"] <- "Jockey-3_Dsim_A"
jock.tree1$tip.label[jock.tree1$tip.label=="Jockey-3_Dmau_mel-like"] <- "Jockey-3_Dmau_A"
jock.tree1$tip.label[jock.tree1$tip.label=="Dyak|CM028598_2|22413314-22416862|+|"] <- "Jockey-3_Dyak_fragment1"
jock.tree1$tip.label[jock.tree1$tip.label=="Dyak|CM028598_2|23203166-23205537|+|"] <- "Jockey-3_Dyak_fragment2"
jock.tree1$tip.label[jock.tree1$tip.label=="Dyak|CM028601_2|24301987-24303493|-|"] <- "Jockey-3_Dyak_fragment3"

#### rotate for clarity
jock.tree2 <- ape::rotate(jock.tree1, c(17))


############
#### Figure
############

pdf(paste("~/Documents/RAxML_bipartitions.alignment_jockey_family_ggtree.pdf",sep = ""), width=7,height=5)
ggtree(jock.tree2) +  
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 75 & as.numeric(label) >= 50,  color = "red"), size = 3,) + 
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 75 & as.numeric(label) < 95, color = "yellow"), size = 3) + 
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 95, color = "forestgreen"), size = 3 ) +
  xlim_tree(1.5) + geom_tiplab() + geom_treescale(x = 0.2, y = 9, width = 0.1) +
  scale_color_identity(guide = 'legend', aesthetics = "color",
                       breaks=c("> 50% Support", "> 75% Support", "> 95% Support"),
                       labels=c('> 50% Support'='red', '> 75% Support'='yellow', '> 95% Support'='forestgreen')) +
  ggforce::geom_mark_circle(x= 0.09, y = 13, expand=0.009, fill = "red", alpha=1, color = "red") +
  ggforce::geom_mark_circle(x= 0.09, y = 12, expand=0.009, fill = "yellow", alpha=1, color = "yellow") +
  ggforce::geom_mark_circle(x= 0.09, y = 11, expand=0.009, fill = "forestgreen", alpha=1, color = "forestgreen") +
  annotate("text", x = 0.25, y = 13, label = "> 50% Support") +
  annotate("text", x = 0.25, y = 12, label = "> 75% Support") + 
  annotate("text", x = 0.25, y = 11, label = "> 95% Support")
dev.off()


