#################
#################
#### Phylogenetic tree of all G2/Jockey-3 insertions
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

jock.tree <- read.tree(file="/Users/lucashemmer/Documents/Cecile_simulans/RAxML_bipartitions.alignment_Jockey-3_melsimyak_400_ORF2_mafft_Jockey-1_Dse.automre")


##### removing tips with high divergence, probably from large indels in sequences

tip.lab.to.remove <- c("Dsec|Y_Contig177|0-433|+|", "Dere|utg000005l|2900967-2902314|+|", "Dmau|Y_Contig132|12513-14183|+|")
jock.tree1 <- drop.tip(jock.tree,tip.lab.to.remove)

## rename

jock.tree1$tip.label[jock.tree1$tip.label=="Jockey-1_DSe_072020_-_ORF2"] <- "Jockey-1_DSe"


#### color dataframe 

tip.lab <- jock.tree1$tip.label
species <- rep(NA,length(tip.lab))
dat <- cbind.data.frame(tip.lab,species)
dat$species[grepl("ISO1",tip.lab)] <- "melanogaster"
dat$species[grepl("Dsim",tip.lab)] <- "simulans"
dat$species[grepl("Dmau",tip.lab)] <- "mauritiana"
dat$species[grepl("Dsec",tip.lab)] <- "sechellia"
dat$species[grepl("Dyak",tip.lab)] <- "yakuba"
dat$species[grepl("Dere",tip.lab)] <- "erecta"
dat$species[grepl("Jockey-1_DSe",tip.lab)] <- "Jockey-1_DSe"


dat$col <- NA
dat$col[dat$species == "melanogaster"] <- "#3288BD"
dat$col[dat$species == "simulans"] <- "#5E4FA2"
dat$col[dat$species == "mauritiana"] <- "#66C2A5"
dat$col[dat$species == "sechellia"] <- "#FEE08B"
dat$col[dat$species == "yakuba"] <- "#F46D43"
dat$col[dat$species == "erecta"] <- "#9E0142"
dat$col[grepl("Jockey-1_DSe",tip.lab)] <- "black"


#### finding common edges to color internally if two branches share the same species

mel <- which(dat$species=="melanogaster")
sim <- which(dat$species=="simulans")
sec <- which(dat$species=="sechellia")
mau <- which(dat$species=="mauritiana")
yak <- which(dat$species=="yakuba")
ere <- which(dat$species=="erecta")
jock <- which(dat$species=="Jockey-1_DSe")
mel.edge <- which.edge(jock.tree1, mel)
sim.edge <- which.edge(jock.tree1, sim)
sec.edge <- which.edge(jock.tree1, sec)
mau.edge <- which.edge(jock.tree1, mau)
yak.edge <- which.edge(jock.tree1, yak)
ere.edge <- which.edge(jock.tree1, ere)
jock.edge <- which.edge(jock.tree1, jock)

#Reduce(intersect, list(mel.edge,sim.edge,sec.edge,mau.edge,yak.edge))
tst <- c(mel.edge,sim.edge,sec.edge,mau.edge,yak.edge,ere.edge, jock.edge)
tst1 <- tst[duplicated(tst)]

clcolr <- rep("black", dim(jock.tree$edge)[1]) 
clcolr[mel.edge] <- "#3288BD"
clcolr[sim.edge] <- "#5E4FA2"
clcolr[sec.edge] <- "#FEE08B"
clcolr[mau.edge] <- "#66C2A5"
clcolr[yak.edge] <- "orange" #"#F46D43"
clcolr[ere.edge] <- "#9E0142"
clcolr[jock.edge] <- "black"
clcolr[tst1] <- "black"


############
#### Assigning centromere or non-centromere 
############

numb.vect <- as.data.frame(jock.tree1$tip.label)
colnames(numb.vect) <- "te.name"
numb.vect <- numb.vect %>% separate(col = te.name, into = c("species","contig","span","direction"), sep="\\|", remove=F) %>%
  separate(col = span, into=c("start","end"), sep = "-", remove = F)

numb.vect$start <- as.numeric(numb.vect$start)
numb.vect$end <- as.numeric(numb.vect$end)

numb.vect <- numb.vect %>% 
  dplyr::mutate(centromere = dplyr::case_when( (contig == "Contig79" & start >= 6279 & end <= 49160) ~ "centromere",
                                                (contig == "tig00057289" & start >= 5435 & end <= 8023) ~ "centromere",
                                                (contig == "3R_5" & start >= 19074 & end <= 69826) ~ "centromere",
                                                (contig == "Contig119" & start >= 29780 & end <= 68785) ~ "centromere",
                                                (contig == "Y_Contig26" & start >= 559 & end <= 99154) ~ "centromere",
                                                (contig == "Contig_137" & start >= 27920 & end <= 118193) ~ "centromere",
                                                (contig == "Contig_14" & start >= 101481 & end <= 193347) ~ "centromere",
                                                (contig == "Contig_2" & start >= 960566 & end <= 1121734) ~ "centromere",
                                                (contig == "Contig_124") ~ "centromere",
                                                (contig == "Y_Contig_135" & start >= 1872124 & end <= 2107756) ~ "centromere",
                                                (contig == "Contig147" & start >= 0 & end <= 69783) ~ "centromere",
                                                (contig == "Contig6" & start >= 277419 & end <= 427757) ~ "centromere",
                                                (contig == "Contig51" & start >= 0 & end <= 96603) ~ "centromere",
                                                (contig == "Contig134" ) ~ "centromere",
                                                (contig == "Y_scaffold3" & start >= 2273714 & end <= 2665043) ~ "centromere",
                                                (contig == "Contig183" & start >= 654208 & end <= 831291) ~ "centromere",
                                                (contig == "Contig46" & start >= 22969 & end <= 77085) ~ "centromere",
                                                (contig == "Contig188" & start >= 181986 & end <= 491408) ~ "centromere",
                                                (contig == "Contig40" & start >= 25236 & end <= 89288) ~ "centromere",
                                                (contig == "Y_scaffold1" & start >= 857067 & end <= 1016359) ~ "centromere",
                                                TRUE ~ "other") )

#### assign shape

numb.vect$pch <- NA
#numb.vect$pch[numb.vect$centromere=="centromere"] <- 4
numb.vect$pch[numb.vect$centromere=="centromere"] <- 8



############
#### Figure
############

pdf(paste("~/Documents/Jockey-3_melsimyak_400_ORF2_mafft_Jockey-1_Dse_unrooted_cent.pdf",sep = ""),width=8,height=8)
plot(jock.tree1, show.tip.label=F, edge.width = 2, edge.color=clcolr,no.margin=T,type="unrooted",label.offset = 0.001, rotate.tree = -45)
tiplabels(cex=0.5, pch=numb.vect$pch, col="black")
box("outer")
legend("topright", legend=c(expression(italic("D. melanogaster")),expression(italic("D. simulans")),expression(italic("D. mauritiana")),expression(italic("D. sechellia")),expression(italic("D. yakuba")),expression(italic("D. erecta"))), col=c("#3288BD","#5E4FA2","#66C2A5","#FEE08B","#F46D43", "#9E0142"), cex=2, pch=c(19,19,19,19), bty='n',lty=NULL)
add.scale.bar(cex = 1, font = 2, col = "black")
text(0.058, 0.35, "Jockey-1_Dse", cex = 1.5)
legend(x = 0, y = 0.15, legend=c("Centromere\nCandidate"), col=c("black"), cex=2, pch=c(4), bty='n',lty=NULL)
dev.off()

