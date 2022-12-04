##############################################################################################
# SUPPLEMENTARY CODE ACCOMPANYING THE PHD THESIS ENTITLED:
##############################################################################################
# "A molecular investigation of stem-galling Tetramesa Walker (Hymenoptera: Eurytomidae) 
# on African grasses: an application to biological control" 
##############################################################################################
# by Clarke van Steenderen, December 2022 / January 2023
# Department of Zoology and Entomology, the Centre for Biological Control
# Rhodes University, Grahamstown/Makhanda, 6139, South Africa
##############################################################################################

# **************************************************
# **************************************************

# CONTENTS OF THIS SCRIPT:
# (Use Ctrl + F to navigate)

###################################################################################################
# CHAPTER 2
###################################################################################################
# PART 1: PHYLOGENY PLOTTING
# PART 2: GMYC SPECIES DELIMITATION TESTS
# PART 3: P-DISTANCES
# PART 4: AMOVAS FOR THE GOLDEN-SHOULDERED TETRAMESA ON ERAGROSTIS CURVULA
# PART 5: MANTEL TESTS: GENETIC DISTANCE VS GEOGRAPHIC DISTANCE

###################################################################################################
# CHAPTER 3
###################################################################################################
# PART 6: COPHYLOGENETIC ANALYSES

###################################################################################################
# CHAPTER 4
###################################################################################################
# PART 7: SPEDE-SAMPLER

###################################################################################################

###################################################################################################
# MISCELLANEOUS
###################################################################################################
# PART 8: MAPS:
# (1) TETRAMESA SPECIES DESCRIPTIONS GLOBALLY
# (2) GRASS DISTRIBUTIONS IN NATIVE AND INVADED RANGES

# **************************************************
# **************************************************

##############################################################################################
# REQUIRED PACKAGES
##############################################################################################

library(ape)
library(pegas)
library(vegan)
library(splits)
library(rlist)
library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(magrittr)
library(dplyr)
library(strap)
#remotes::install_github("fmichonneau/phyloch")
library(phyloch)
library(magrittr)
library(janitor)
library(ggrepel)
library(gridExtra)
library(readxl)
library(tidyverse)
library(tidyr)
library(splitstackshape)
library(cowplot)

##############################################################################################
##############################################################################################
##############################################################################################
# PART 1: PHYLOGENY PLOTTING
##############################################################################################
##############################################################################################
##############################################################################################

##############################################################################################
# Upload the sample sheet of sequence information
##############################################################################################

# extract host plant information
host.plants = readxl::read_xlsx("source_modifiers.xlsx", sheet = "host_plants")
host.plants = janitor::clean_names(host.plants)
# extract the full sample sheet
source.mods = readxl::read_xlsx("source_modifiers.xlsx", sheet = "cleaned_tetramesa_clarke")
source.mods = janitor::clean_names(source.mods)

str(host.plants)
categories = colnames(host.plants) ; categories
head(host.plants)

# change the NA values to zeroes, and numeric values to characters
host.plants[is.na(host.plants)] <- 0
host.plants[] <- lapply(host.plants, as.character)
# change the layout of the data to make it compatible with the upcoming functions
host.plants = tibble::column_to_rownames(host.plants, var = "sample_id")

###############
# 28S
###############

##############################################################################################
# Read in the 28S phylogenetic trees from MrBayes and IQTREE
##############################################################################################

iqtree_28S = treeio::read.mrbayes("phylogenies/28S_iqtree.nex")
mrbayes_28S = treeio::read.mrbayes("phylogenies/28S_mrbayes.nex")

##############################################################################################
# Generate the ML phylogeny
##############################################################################################

iqtree_plot_28S = iqtree_28S %>% ggtree::ggtree(., color="black", layout = "rectangular", lwd=0.1) + 
  geom_tiplab(size=1, color="black", font = 1.5) +
  geom_label2(aes(subset=!isTip, label=round(as.numeric(boots),2)), size=2, color="black", alpha=0, label.size = 0, nudge_x = -0.005, nudge_y = 1) +  # add node numbers  +
  geom_treescale(fontsize=4, linesize=1)

iqtree_plot_28S

##############################################################################################
# Generate the Bayesian phylogeny
##############################################################################################

mrbayes_plot_28S = mrbayes_28S %>% ggtree::ggtree(., color="black", layout = "rectangular", lwd=0.1) + 
  geom_tiplab(size=1, color="black", font = 1.5) +
  geom_label2(aes(subset=!isTip, label=round(as.numeric(prob),2)), size=2, color="black", alpha=0, label.size = 0, nudge_x = -0.005, nudge_y = 1) +  # add node numbers  +
  geom_treescale(fontsize=4, linesize=1)

mrbayes_plot_28S 

##############################################################################################
# Add a coloured matrix to the right of the phylogeny, indicating which host plant each sample
# is associated with
##############################################################################################

heatmap_28S =  ggtree::gheatmap(mrbayes_plot_28S, host.plants, 
                                color="black",
                                offset=0.002, width=0.2, colnames_position = "top", 
                                colnames_angle = 90, font.size=1,  
                                colnames_offset_y = 2) + 
  scale_fill_manual(values=c("white", "darkblue")) +
  theme(legend.position="bottom") 

plot(heatmap_28S)

###############
# COI
###############

##############################################################################################
# Read in the COI phylogenetic trees from MrBayes and IQTREE
##############################################################################################

COI_iqtree = treeio::read.mrbayes("phylogenies/COI_iqtree.nex")
COI_mrbayes = treeio::read.mrbayes("phylogenies/COI_mrbayes.nex")

##############################################################################################
# Generate the ML phylogeny
##############################################################################################

COI_iqtree_plot = COI_iqtree %>% ggtree::ggtree(., color="black", layout = "rectangular", lwd=0.1) + 
  geom_tiplab(size=1, color="black", font = 1.5) +
  geom_label2(aes(subset=!isTip, label=round(as.numeric(boots),2)), size=2, color="black", alpha=0, label.size = 0, nudge_x = -0.005, nudge_y = 1) +  # add node numbers  +
  geom_treescale(fontsize=4, linesize=1)

COI_iqtree_plot

##############################################################################################
# Generate the Bayesian phylogeny
##############################################################################################

COI_mrbayes_plot = COI_mrbayes %>% ggtree::ggtree(., color="black", layout = "rectangular", lwd=0.1) + 
  geom_tiplab(size=1, color="black", font = 1.5) +
  geom_label2(aes(subset=!isTip, label=round(as.numeric(prob),2)), size=2, color="black", alpha=0, label.size = 0, nudge_x = -0.005, nudge_y = 1) +  # add node numbers  +
  geom_treescale(fontsize=4, linesize=1)

COI_mrbayes_plot

##############################################################################################
# Add a coloured matrix to the right of the phylogeny, indicating which host plant each sample
# is associated with
##############################################################################################

COI_heatmap =  ggtree::gheatmap(COI_mrbayes_plot, host.plants, 
                                color="black",
                                offset=0.002, width=0.2, colnames_position = "top", 
                                colnames_angle = 90, font.size=1,  
                                colnames_offset_y = 2) + 
  scale_fill_manual(values=c("white", "darkblue")) +
  theme(legend.position="bottom") 

plot(COI_heatmap)

# This is code to remove the outgroups from the tree, and write the result to a tree file. This is used later for species delimitation tests
mrbayes_COI_noout = ape::drop.tip(COI_mrbayes, c("AY317231", "AY317230", "AY317232") )
ape::write.tree(mrbayes_COI_noout, "mrbayes_COI_noout.nwk", digits = 0)

##############################################################################################
# BEAST DATED PHYLOGENY - 28S gene using a standard molecular clock rate, with no 
# fossil calibrations
##############################################################################################

beast_28S = phyloch::read.beast("phylogenies/28S_beast_mcc")
beast_28S$root.time = 373.1

strap::geoscalePhylo(tree=ladderize(beast_28S,right=FALSE), units=c("Period", "Epoch"),
                     boxes="Epoch", cex.tip=0.15, cex.age=0.5, cex.ts=0.7, label.offset=0.5, 
                     lwd=1, width=0.6, ts.col = F, tick.scale = 10)

HPDbars(beast_28S, lwd = 0.6, col = "cornflowerblue")

#LTT PLOTS (Lineage through time)
ltt.coplot(beast_28S, show.tip.label = FALSE)
# or using the ape equivalent
ape::ltt.plot(beast_28S, backward = T, log="y", col = "black", lwd = 3)

##############################################################################################
##############################################################################################
##############################################################################################
# PART 2: GMYC SPECIES DELIMITATION TESTS
##############################################################################################
##############################################################################################
##############################################################################################

####################
# 28S
####################

gmyc_28S_beasttree = ape::read.tree("species_delimitation/28S_beast_gmyc.nwk")
# remove outgroups
gmyc_28S_beasttree_no_outgroup = drop.tip(gmyc_28S_beasttree, c("AY317177", "AY317165", "AY317156", "AY317161") )
gmyc.28S = splits::gmyc(gmyc_28S_beasttree_no_outgroup, quiet = F, method = "single")

# GMYC summary statistics

splits::summary.gmyc(gmyc.28S)
splits::spec.list(gmyc.28S)
support = splits::gmyc.support(gmyc.28S, p = 0.95) ; support
support[support==0] <- NA
support = support[support >= 0.85] ; support

plot(gmyc.28S$tree, cex = 0.2)
nodelabels(round(support, 2), cex=0.6, frame="n", adj=1)

splits::confset.gmyc(gmyc.28S)
splits::plot.gmyc(gmyc.28S)

####################
# COI
####################

gmyc_COI_beasttree = ape::read.tree("species_delimitation/COI_beast_gmyc.nwk")
# remove outgroups
gmyc_COI_beasttree_no_outgroup = drop.tip(gmyc_COI_beasttree, c("AY317231", "AY317230", "AY317232") )
gmyc.COI = splits::gmyc(gmyc_COI_beasttree_no_outgroup, quiet = F, method = "single")

splits::summary.gmyc(gmyc.COI)
splits::spec.list(gmyc.COI)
support = splits::gmyc.support(gmyc.COI, p = 0.95) ; support
support[support==0] <- NA
support = support[support >= 0.85] ; support

plot(gmyc.COI$tree, cex = 0.2)
nodelabels(round(support, 2), cex=0.6, frame="n", adj=1)

splits::confset.gmyc(gmyc.COI)
splits::plot.gmyc(gmyc.COI)

##############################################################################################
##############################################################################################
##############################################################################################
# PART 3: P-DISTANCES
##############################################################################################
##############################################################################################
##############################################################################################

####################
# 28S
####################

##########################
# INTERSPECIFIC DISTANCES
##########################

# interspecific distances, obtained from MEGA
inter_pdists_28S = read.csv("p_distances/28S_interspecific_dists.csv", row.names = 1)
# reformat the data
inter_pdists_28S = reshape2::melt(as.matrix(inter_pdists_28S), na.rm = TRUE) %>% janitor::clean_names()
# some quick summary statistics for the ingroup data only
inter_pdists_28S_ingroup = dplyr::filter(inter_pdists_28S, var1 != "OUTGROUP", var2 != "OUTGROUP")
summary(inter_pdists_28S_ingroup$value)

################################
# Plot a heatmap of p-distances
################################
                       
#######################
# including the outgroup
#######################

ggplot(inter_pdists_28S, aes(x = var1, y=var2, fill=value, label=round(value*100,3)))+
  geom_tile()+
  geom_text(size=3,color="white") +
  xlab("") +
  ylab("") +
  #scale_fill_gradient(low = 'lightblue', high = 'darkblue') +
  scale_fill_gradient2(low = "darkblue", mid = "lightblue", high = "darkred", midpoint = 0.125) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

###########################################################
# excluding the outgroup, and colouring by distance range
###########################################################

ggplot(inter_pdists_28S, aes(x = var1, y=var2, fill=value, label=round(value*100,3)))+
  #geom_tile(fill = ifelse(p_dists_28S$value > 0.035 && p_dists_28S$value < 0.05, "yellow", NA), color = "black", alpha = 0.5) +
  geom_tile(fill = ifelse(inter_pdists_28S$value > 0.07, "blue", NA), color = "black", alpha = 0.5) +
  geom_tile(fill = ifelse(inter_pdists_28S$value < 0.035, "red", NA), color = "black", alpha = 0.5) +
  geom_text(size=3,color="black") +
  xlab("") +
  ylab("") +
  #scale_fill_gradient(low = 'blue', high = 'red') +
  #scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.125) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

##########################
# INTRASPECIFIC DISTANCES
##########################

# Intraspecific distances from MEGA
intra_pdists_28S = read.csv("p_distances/28S_intraspecific_dists.csv")

intra_pdists_28S$p_dist = as.numeric(intra_pdists_28S$p_dist)
intra_pdists_28S = na.omit(intra_pdists_28S)

##############################################################################
# Extract only the golden and non-golden-shouldered clades of interest
##############################################################################
# GOLDEN-SHOULDERED

# intraspecific
GS_only_28S_intra = intra_pdists_28S %>% 
  dplyr::filter(group == "GS rigidior" | group == "GS" | group == "GS sporob" | group == "ngs hirta" | group == "romana")
# interspecific
GS_only_28S_inter = inter_pdists_28S %>% 
  dplyr::filter(var1 == "agayanus" | var1 == "GS" | var1 == "GS_rigidior" | var1 == "GS_sporob" | var1 == "ngs_hirta" | var1 == "romana" | var1 == "t_bambusae")

#NON-GOLDEN SHOULDERED

NGS_only_28S_intra = intra_pdists_28S %>%
  dplyr::filter(group == "NGS"| group == "NGS1"| group == "NGS2"| group == "NGS3"| group == "NGS4" | group == "ngs_SPY" | group == "ngs_hirta2")

NGS_only_28S_inter= inter_pdists_28S %>%
  dplyr::filter(var1 == "NGS" | var1 == "NGS1" | var1 == "NGS2" | var1 == "NGS3" | var1 == "NGS4")

# set colours for A. gayanus and H. hirta
colors_gayanus_28S = ifelse(GS_only_28S_inter$var2 == "agayanus" | GS_only_28S_inter$var2 =="ngs_hirta",
                            "agayanus/hhir", "other")

##############################################################################
# BARCODE GAP PLOT USING INTRA AND INTERSPECIFIC P DISTANCES
##############################################################################

GS_gap_plot = ggplot() +
  # intraspecific GS
  geom_histogram(data = GS_only_28S_intra, aes(x = p_dist), bins = 35, fill = "orange", alpha = 1) +
  # interspecific GS
  geom_histogram(data = GS_only_28S_inter, aes(x = value, fill = colors_gayanus_28S), bins = 35, alpha = 0.75 ) +
  scale_fill_manual(values = c("black", "lightblue"))+
  scale_x_continuous(breaks = seq(0, 0.15, by = 0.01)) +
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 14) ) +
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 14) ) +
  ylab("Frequency") + 
  xlab("P-distance (K2P)") +
  ggtitle("28S intra and interspecific distances") ; GS_gap_plot

##############################################################################
# THRESHOLD EXPERIMENTATION - FINDING GENETIC BOUNDARIES
##############################################################################

# Select a range of threshold values to start with
THRESH1 = 0.01
THRESH2 = 0.015
THRESH3 = 0.02
THRESH4 = 0.025

groups_28S = levels(inter_pdists_28S$var2)

# Have a look at which groups (via pairwise comparisons) are conspecific when a particular p-distance threshold is set
for(i in 1:length(groups_28S)){
  
  res =  inter_pdists_28S %>%
    dplyr::filter(var2 == groups_28S[i] & value <= THRESH1)
  
  if(length(row.names(res) >= 1))
    print(res)
  
}

# Visualise this on plots - everything below the threshold line can be considered "conspecific", or within the same taxonomic group

# threshold = 0.01
threshold_plot_1 = inter_pdists_28S %>%
  dplyr::filter(var1 != "OUTGROUP") %>%
  ggplot(., aes(x = var2, y = value)) +
  geom_point(size = 1) +
  ylab("P-distance (K2P)") + 
  xlab("Phylogenetic clade") +
  geom_hline(yintercept = THRESH1, linetype="dashed", color = "darkred") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(margin = margin(r = 20)) ) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.01)) +
  geom_text_repel( aes( label=ifelse(value <= THRESH1, as.character(var1), "") ), cex = 1, col = "darkred", max.overlaps = 100 ) +
  ggtitle("A")

# threshold = 0.015
threshold_plot_2 = inter_pdists_28S %>%
  dplyr::filter(var1 != "OUTGROUP") %>%
  ggplot(., aes(x = var2, y = value)) +
  geom_point(size = 1) +
  ylab("P-distance (K2P)") + 
  xlab("Phylogenetic clade") +
  geom_hline(yintercept = THRESH2, linetype="dashed", color = "darkred") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(margin = margin(r = 20)) ) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.01)) +
  geom_text_repel( aes( label=ifelse(value <= THRESH2, as.character(var1), "") ), cex = 1, col = "darkred", max.overlaps = 100 ) +
  ggtitle("B")

# threshold = 0.02
threshold_plot_3 = inter_pdists_28S %>%
  dplyr::filter(var1 != "OUTGROUP") %>%
  ggplot(., aes(x = var2, y = value)) +
  geom_point(size = 1) +
  ylab("P-distance (K2P)") + 
  xlab("Phylogenetic clade") +
  geom_hline(yintercept = THRESH3, linetype="dashed", color = "darkred") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(margin = margin(r = 20)) ) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.01)) +
  geom_text_repel( aes( label=ifelse(value <= THRESH3, as.character(var1), "") ), cex = 1, col = "darkred", max.overlaps = 100 ) +
  ggtitle("C")

# threshold = 0.025
threshold_plot_4 = inter_pdists_28S %>%
  dplyr::filter(var1 != "OUTGROUP") %>%
  ggplot(., aes(x = var2, y = value)) +
  geom_point(size = 1) +
  ylab("P-distance (K2P)") + 
  xlab("Phylogenetic clade") +
  geom_hline(yintercept = THRESH4, linetype="dashed", color = "darkred") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(margin = margin(r = 20)) ) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.01)) +
  geom_text_repel( aes( label=ifelse(value <= THRESH4, as.character(var1), "") ), cex = 1, col = "darkred", max.overlaps = 100 ) +
  ggtitle("D")

# plot them all together
gridExtra::grid.arrange(threshold_plot_1, threshold_plot_2, threshold_plot_3, threshold_plot_4)


####################
# COI
####################

##########################
# INTERSPECIFIC DISTANCES
##########################

inter_pdists_COI = read.csv("p_distances/COI_interspecific_dists.csv", row.names = 1)
inter_pdists_COI = reshape2::melt(as.matrix(inter_pdists_COI), na.rm = TRUE)
inter_pdists_COI = janitor::clean_names(inter_pdists_COI)
inter_pdists_COI_ingroup = dplyr::filter(inter_pdists_COI, var1 != "OUTGROUPS", var2 != "OUTGROUPS")
summary(inter_pdists_COI_ingroup$value)

# Create a heatmap with conditional colouring

ggplot(inter_pdists_COI %>% dplyr::filter(var1 != "GS.ETRICH.2", var2 != "GS.ETRICH.2"), 
       aes(x = var1, y=var2, fill=value, label=round(value*100, 1)))+
  geom_tile(fill = ifelse(inter_pdists_COI$value <= 0.03, "red", NA), color = "black", alpha = 1) +
  geom_tile(fill = ifelse(inter_pdists_COI$value <= 0.05, "red", NA), color = "black", alpha = 0.5) +
  geom_tile(fill = ifelse(inter_pdists_COI$value <= 0.11, "red", NA), color = "black", alpha = 0.2) +
  geom_tile(fill = ifelse(inter_pdists_COI$value > 0.11 & inter_pdists_COI$value < 0.15, "orange", NA), color = "black", alpha = 0.6) +
  geom_tile(fill = ifelse(inter_pdists_COI$value >= 0.15 & inter_pdists_COI$value < 0.18, "white", NA), color = "black", alpha = 0.2) +
  geom_tile(fill = ifelse(inter_pdists_COI$value >= 0.18, "blue", NA), color = "black", alpha = 0.4) +
  geom_text(size=2,color="black") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

##########################
# INTRASPECIFIC DISTANCES
##########################

intra_pdists_COI = read.csv("p_distances/COI_intraspecific_dists.csv")
intra_pdists_COI$p_dist = as.numeric(intra_pdists_COI$p_dist)
intra_pdists_COI = na.omit(intra_pdists_COI)

##############################################################################
# Extract only the golden and non-golden-shouldered clades of interest
##############################################################################
# GOLDEN-SHOULDERED RELATIVE TO OTHER GROUPS

GS_only_COI_intra = intra_pdists_COI %>% dplyr::filter(group == "agayanus" | group == "ETEF" | group == "GS.EBI" | group == "GS.ECUR.1" | group == "GS.ECUR.2" | group == "GS.ECUR.3"
                                                        | group == "GS.ECUR.4" | group == "GS.EPLAN" | group == "GS.ERIG.1" | group == "GS.ERIG.2" | group == "GS.ETRICH.1" | group == "GS.SPOR.NAT"
                                                        | group == "GS.SPOR.PYR" | group == "ROM" | group == "NGS.HHIR.2" | group == "NGS.HHIR.3")

GS_only_COI_inter = inter_pdists_COI %>% 
  dplyr::filter(var1 == "agayanus" | var1 == "ETEF" | var1 == "GS.EBI" | var1 == "GS.ECUR.1" | var1 == "GS.ECUR.2" | var1 == "GS.ECUR.3"
                | var1 == "GS.ECUR.4" | var1 == "GS.EPLAN" | var1 == "GS.ERIG.1" | var1 == "GS.ERIG.2" | var1 == "GS.ETRICH.1" | var1 == "GS.SPOR.NAT"
                | var1 == "GS.SPOR.PYR" | var1 == "ROM" | var1 == "NGS.HHIR.2" | var1 == "NGS.HHIR.3" | var1 == "T.bambusae")

colors_coi = ifelse(GS_only_COI_inter$var2 == "agayanus" | GS_only_COI_inter$var2 =="HHIR.2" | GS_only_COI_inter$var2 =="HHIR.3",
                        "agayanus/hhir", "other")

##############################################################################
# BARCODE GAP PLOT 
##############################################################################

coi_barcode_gap = ggplot() +
  geom_histogram(data = GS_only_COI_intra, aes(x = p_dist), bins = 40, fill = "orange", alpha = 0.95) +
  geom_histogram(data = GS_only_COI_inter, aes(x = value, fill = colors_coi), bins = 35,  alpha = 0.75 ) +
  scale_fill_manual(values = c("black", "lightblue")) +
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.02)) +
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 14) ) +
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 14) ) +
  ylab("Frequency") + 
  xlab("P-distance (K2P)") +
  ggtitle("COI") ; coi_barcode_gap

##############################################################################
# THRESHOLD EXPERIMENTATION - FINDING GENETIC BOUNDARIES
##############################################################################

# Select a range of threshold values to start with
THRESHA = 0.03
THRESHB = 0.04
THRESHC = 0.05
THRESHD = 0.06

groups_COI = levels(inter_pdists_COI$var2)

# have a look at which groups (via pairwise comps) are conspecific when a particular p-dist threshold is set
for(i in 1:length(groups_COI)){
  
  res =  inter_pdists_COI %>%
    dplyr::filter(var2 == groups_COI[i] & value <= THRESHA)
  
  if(length(row.names(res) >= 1))
    print(res)
  
}

# Visualise this on plots - everything below the threshold line can be considered "conspecific", or within the same taxonomic group

# threshold = 0.03

threshold_plotA = inter_pdists_COI %>%
  dplyr::filter(var1 != "OUTGROUP") %>%
  ggplot(., aes(x = var2, y = value)) +
  geom_point(size = 1) +
  ylab("P-distance (K2P)") + 
  xlab("Phylogenetic clade") +
  geom_hline(yintercept = THRESHA, linetype="dashed", color = "royalblue") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(margin = margin(r = 20)) ) +
  scale_y_continuous(breaks = seq(0, 0.32, by = 0.02)) +
  ggrepel::geom_text_repel( aes( label=ifelse(value <= THRESHA, as.character(var1), "") ), cex = 1, col = "royalblue", max.overlaps = 100 ) +
  ggtitle("A")

# threshold = 0.04

threshold_plotB = inter_pdists_COI %>%
  dplyr::filter(var1 != "OUTGROUP") %>%
  ggplot(., aes(x = var2, y = value)) +
  geom_point(size = 1) +
  ylab("P-distance (K2P)") + 
  xlab("Phylogenetic clade") +
  geom_hline(yintercept = THRESHB, linetype="dashed", color = "royalblue") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(margin = margin(r = 20)) ) +
  scale_y_continuous(breaks = seq(0, 0.32, by = 0.02)) +
  ggrepel::geom_text_repel( aes( label=ifelse(value <= THRESHB, as.character(var1), "") ), cex = 1, col = "royalblue", max.overlaps = 100 ) +
  ggtitle("B")

# threshold = 0.05

threshold_plotC = inter_pdists_COI %>%
  dplyr::filter(var1 != "OUTGROUP") %>%
  ggplot(., aes(x = var2, y = value)) +
  geom_point(size = 1) +
  ylab("P-distance (K2P)") + 
  xlab("Phylogenetic clade") +
  geom_hline(yintercept = THRESHC, linetype="dashed", color = "royalblue") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(margin = margin(r = 20)) ) +
  scale_y_continuous(breaks = seq(0, 0.32, by = 0.02)) +
  ggrepel::geom_text_repel( aes( label=ifelse(value <= THRESHC, as.character(var1), "") ), cex = 1, col = "royalblue", max.overlaps = 100 ) +
  ggtitle("C")

# threshold = 0.06

threshold_plotD = inter_pdists_COI %>%
  dplyr::filter(var1 != "OUTGROUP") %>%
  ggplot(., aes(x = var2, y = value)) +
  geom_point(size = 1) +
  ylab("P-distance (K2P)") + 
  xlab("Phylogenetic clade") +
  geom_hline(yintercept = THRESHD, linetype="dashed", color = "royalblue") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(margin = margin(r = 20)) ) +
  scale_y_continuous(breaks = seq(0, 0.32, by = 0.02)) +
  ggrepel::geom_text_repel( aes( label=ifelse(value <= THRESHD, as.character(var1), "") ), cex = 1, col = "royalblue", max.overlaps = 100 ) +
  ggtitle("D")

# plot them all together
gridExtra::grid.arrange(threshold_plotA, threshold_plotB, threshold_plotC, threshold_plotD)

##############################################################################################
##############################################################################################
##############################################################################################
# PART 4: AMOVAS FOR THE GOLDEN-SHOULDERED TETRAMESA ON ERAGROSTIS CURVULA
##############################################################################################
##############################################################################################
##############################################################################################

#########################
# 28S SEQUENCE DATA
#########################

seqs_28S = ape::read.FASTA("fasta_files/28S_only_aligned.fasta")

#########################
# COI SEQUENCE DATA
#########################

seqs_COI = ape::read.FASTA("fasta_files/COI_only_aligned.fas")

# Extract only the golden-shouldered Tetramesa from E. curvula
GS_eragrostis_curvula = dplyr::filter(source.mods, morphotype == "GS" & grepl("Eragrostis curvula", host))
# create an empty fasta file
file.create("COI_all_GS_eragrostis_curvula.fasta")

# loop through the subset, and generate a fasta files with these sequences
for(i in 1:nrow(GS_eragrostis_curvula)){
  
  for(j in 1:length(seqs_COI)){
    if(GS_eragrostis_curvula$sequence_id[i] == labels(seqs_COI[j]))
      ape::write.FASTA( seqs_COI[j], "COI_all_GS_eragrostis_curvula.fasta", append = T )
  }
  
}

# get the sample info for each sequence
GS_eragrostis_curvula_location_and_host = dplyr::select(GS_eragrostis_curvula, sequence_id, province, host)
# remove some samples
GS_eragrostis_curvula_location_and_host = GS_eragrostis_curvula_location_and_host[-c(5, 11, 22), ]
# write this summary to an excel file
write.csv(GS_eragrostis_curvula_location_and_host, "GS_curvula_hosts_localities.csv", row.names = F)
# find how many provinces are represented
curvula_province = factor(GS_eragrostis_curvula_location_and_host$province)

# read the FASTA file back in 
ecurvula_seqs = ape::read.FASTA("COI_all_GS_eragrostis_curvula.fasta")
# generate a distance matrix
d = ape::dist.dna(ecurvula_seqs)
class(d)

# AMOVA with PEGAS library
amova_result = pegas::amova(d~curvula_province, nperm = 10000)
amova_result$tab
pegas::write.pegas.amova(amova_result, file = "pegas.results.csv")

##############################################################################################
##############################################################################################
##############################################################################################
# PART 5: MANTEL TESTS: GENETIC DISTANCE VS GEOGRAPHIC DISTANCE
##############################################################################################
##############################################################################################
##############################################################################################

library(geosphere)
library(ade4)
library(ggplot2)
library(ggmatplot)

# read in the COI p-distance matrix obtained from MEGA for the golden-shouldered curvula wasps
COI_GS_curvula_p_dists = read.csv("mantel_test/coi_gs_curvula_distance_data.csv", row.names = 1)

# read in the GPS coordinates for each sequence
COI_GS_curvula_GPS = read.csv("mantel_test/gps_coords.csv", row.names = 1)

# convert the GPS coordinates into a pairwise distance matrix
dist_mat = geosphere::distm(COI_GS_curvula_GPS, fun = distGeo)

# run the Mantel test using the vegan package
mantel_test = vegan::mantel(
  as.dist(COI_GS_curvula_p_dists),
  as.dist(dist_mat),
  method = "spearman",
  permutations = 9999,
  na.rm = TRUE
)

mantel_test

# Mantel test statistic (this is the R2 value)
mantel_test$statistic

# Mantel correlation p-value
mantel_test$signif

# plot distance against genetic k2P distance
mantel_curvula = ggmatplot::ggmatplot( as.matrix( dist_mat/1000 ), as.matrix(COI_GS_curvula_p_dists), color = "black", shape = 16 ) +
  #geom_abline(intercept = 1.367e-02, slope = 2.536e-05, color = "red") +
  theme_classic() +
  theme(legend.position="none") +
  xlab("Distance (km)") + 
  ylab("Genetic K2P distance") +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 16, face="bold") ) +
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 16, face="bold") ) +
  theme(text = element_text(size=16) ) +
  ylim(0, 0.08) +
  scale_x_continuous(breaks = seq(0, 1200, by =200)) ; mantel_curvula

##############################################################################################
##############################################################################################
##############################################################################################
# PART 6: COPHYLOGENETIC ANALYSES
##############################################################################################
##############################################################################################
##############################################################################################

####################################################################
# A general cophylogenetic plot between wasps and their grass hosts,
# using grass chloroplast genes and wasp COI sequence data
####################################################################

# read in the full grass ML phylogeny - chloroplast genes
TreeGrass.full <- ape::read.nexus("cophylogenies/grasses_chloroplast.nex")
# read in the full COI wasp ML phylogeny
TreeTetramesa.full <- ape::read.nexus("cophylogenies/COI_tetramesa_all.nex")

# read in the tip label linkers
dat = as.matrix(read.csv("cophylogenies/COI_tiplabs.csv", header = F))

# generate a cophylogenetic plot
tetra_grass_cophylo.full = phytools::cophylo(TreeTetramesa.full, TreeGrass.full, assoc = dat)

# plot the figure with some tweaks
cophylo_plot = plot(tetra_grass_cophylo.full, link.lwd=2.5, link.lty="solid",
                    link.col=make.transparent("red", 0.3), fsize = c(0.1, 0.1),
                    pts=FALSE)

# create an svg image file
svg("coplot_chloroplast.svg")
plot(tetra_grass_cophylo.full, link.lwd=2.5, link.lty="solid",
     link.col=make.transparent("red", 0.3), fsize = c(0.1, 0.1),
     pts=FALSE)
dev.off()

####################################################################
# Run a formal cophylogenetic analysis to look for statistical 
# relationships
####################################################################

library(paco)

# read in a tree of representative sequences for the grasses and wasps
host_grass = ape::read.tree("cophylogenies/grass_representatives.nwk")
tet_parasite = ape::read.tree("cophylogenies/COI_tetramesa_representatives.nwk")
plot(tet_parasite)
# read in a CSV file with the herbivore-host associations
links = read.csv("cophylogenies/co_matrix.csv", row.names=1)

# set up the analysis
host_grass_co = ape::cophenetic.phylo(host_grass)
tet_parasite_co = ape::cophenetic.phylo(tet_parasite)
D = paco::prepare_paco_data(H = host_grass_co, P = tet_parasite_co, HP = links)
D = paco::add_pcoord(D, correction = "cailliez")

# run the cophylogenetic analysis
D = paco::PACo(D, nperm = 1000, seed = 12, method = "r2", symmetric = FALSE)

# extract the residuals
res = paco::residuals_paco(D$proc) ; res
res_df = as.data.frame(res)

# goodness of fit
D$gof

# Make a list out of the interaction matrix
association <- data.frame(tetra=rownames(links)[which(links==1, arr.ind=TRUE)[,"row"]], grass=colnames(links)[which(links==1, arr.ind=TRUE)[,"col"]])
weight <- (res^-2)/1500
# plot an interaction matrix, where lines are thickened based on their residual value
ape::cophyloplot(y = ladderize(tet_parasite, right = F), x = ladderize(host_grass, right = F), 
                 assoc = association, show.tip.label=T, use.edge.length=FALSE,
                 lwd=weight, col='steelblue', length.line=0, gap=-27, space=40, type = "phylogram")

# Visualise residuals 
ggplot(res_df, aes(x=res))+
  geom_density(fill='grey80')+
  theme_bw()+
  geom_vline(data=res_df, aes(xintercept=res), col='darkorange1')+
  theme(panel.grid=element_blank())+
  xlab('Procrustes residuals')+
  ylab('Frequency')+
  scale_x_continuous(limits=c(0, 0.075))+
  scale_y_continuous(limits=c(0, 45)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 14, face="bold") ) +
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 14, face="bold") ) +
  theme(text = element_text(size=16) ) 

# extract the residuals of all the Eragrostis wasps
eragrostis_resids = res_df %>% dplyr::filter(grepl("Eragrostis", rownames(res_df))) ; eragrostis_resids

# extract the residuals of the GS and NGS Eragrostis wasps
eragrostis_resids_GS = eragrostis_resids %>% dplyr::filter(grepl('_GS', rownames(eragrostis_resids))) ; eragrostis_resids_GS
eragrostis_resids_NGS = eragrostis_resids %>% dplyr::filter(grepl('_NGS', rownames(eragrostis_resids))) ; eragrostis_resids_NGS

# run some t-tests to see whether there are significant differences between the GS and NGS wasps on Eragrostis
t.test(eragrostis_resids_GS, eragrostis_resids_NGS, paired = FALSE, alternative = "less")
t.test(eragrostis_resids_NGS, eragrostis_resids_GS, paired = FALSE, alternative = "greater")

# extract the residuals for Sporobolus and Arundo wasps
spor_resids = res_df %>% filter(grepl('Sporobolus', rownames(res_df)))
arundo_resids = res_df %>% filter(grepl('Arundo', rownames(res_df)))

# create a matrix of all the residuals of interest
resid_combo = rbind(data.frame(eragrostis_resids_GS, level = "Eragrostis GS"), 
                    data.frame(eragrostis_resids_NGS, level = "Eragrostis NGS"),
                    data.frame(spor_resids, level = "S. pyramidalis"),
                    data.frame(arundo_resids, level = "A. donax") )

# generate a boxplot for the residuals of each broad group
ggplot( data = resid_combo, aes(x = level, y = res, fill = level) ) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_brewer(palette='Paired')+
  scale_x_discrete(labels=c("A. donax", 'Eragrostis GS','Eragrostis NGS', "S. pyramidalis" )) +
  xlab("") +
  ylab("Procrustes residuals") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 14, face="bold") ) +
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 14, face="bold") ) +
  theme(text = element_text(size=16) ) 

##############################################################################################
##############################################################################################
##############################################################################################
# PART 7: SPEDE-SAMPLER
##############################################################################################
##############################################################################################
##############################################################################################

#################
#################
# COI RESULTS
################
#################

data_sheet_COI = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/COI_info_sheet.csv")
length( unique(data_sheet_COI$morphospecies) )

count(data_sheet_COI, "morphospecies")

# CLUSTERS

clusters = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/clusters_amalgamated.csv", header = T, check.names = F)
clusters_melt = reshape2::melt(clusters)
colnames(clusters_melt) = c("file_name", "data_perc", "clusters")

# ENTITIES

entities = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/entities_amalgamated.csv", header = T, check.names = F)
entities_melt = reshape2::melt(entities)
colnames(entities_melt) = c("file_name", "data_perc", "entities")

clusters_melt$entities = entities_melt$entities 

clusts_ents = reshape::melt(clusters_melt)

# BOXPLOT

clust_ent_boxplot = ggplot(data = clusts_ents, aes(x = data_perc, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  stat_summary(fun = mean, aes(group = variable), geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  scale_fill_manual(values=c("#1EAD4B", "#AC72D9"), name = "", labels = c("Clusters", "Entities")) +
  xlab("Data Percentage") +
  ylab("Number of clusters or entities") +
  ggtitle("A") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  geom_hline(yintercept= 14, linetype="dashed") ; clust_ent_boxplot
#guides(fill=guide_legend(title=""))

# PLOT JUST CLUSTERS, WITH AN ACCUMULATION CURVE

clust_accum = ggplot(data = clusters_melt, aes(x = as.numeric( data_perc ), y = clusters)) + 
  geom_point() +
  geom_smooth(span = 3, col = "blue", fill = "lightblue", alpha = 0.5) +
  #geom_boxplot(fill = "#1EAD4B") +
  #stat_summary(fun = mean, geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  xlab("Data Percentage") +
  ylab("Number of clusters") +
  ggtitle("A1") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  ylim(0,25) +
  scale_x_continuous(breaks = as.numeric(clusters_melt$data_perc), labels = clusters_melt$data_perc) +
  geom_hline(yintercept= 14, linetype="dashed") ; clust_accum

# PLOT JUST ENTITIES, WITH AN ACCUMULATION CURVE

ent_accum = ggplot(data = clusters_melt, aes(x = as.numeric( data_perc ), y = entities)) + 
  geom_point() +
  geom_smooth(span = 3, col = "blue", fill = "lightblue", alpha = 0.5) +
  #geom_boxplot(fill = "#1EAD4B") +
  #stat_summary(fun = mean, geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  xlab("Data Percentage") +
  ylab("Number of entities") +
  ggtitle("A2") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  ylim(0,25) +
  geom_hline(yintercept= 14, linetype="dashed") ; ent_accum

# SPLITTING RATIOS

splitting = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/oversplitting_ratio_amalgamated.csv", header = T, check.names = F)
splitting_melt = reshape2::melt(splitting)
colnames(splitting_melt) = c("file_name", "data_perc", "splitting_ratio")

splitting_excl_sing = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/oversplitting_excl_singles_amalgamated.csv", header = T, check.names = F)
splitting_excl_sing_melt = reshape2::melt(splitting_excl_sing)
colnames(splitting_excl_sing_melt) = c("file_name", "data_perc", "splitting_ratio_excl")

splitting_melt$splitting_excl = splitting_excl_sing_melt$splitting_ratio_excl

splitting_combo = reshape::melt( splitting_melt )

# BOXPLOT

splitting_boxplot = ggplot(data = splitting_combo, aes(x = data_perc, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  stat_summary(fun = mean, aes(group = variable), geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  scale_fill_manual(values=c("#2D86F0", "#F02D5A"), name = "", labels = c("+ singletons", "- singletons")) +
  xlab("Data Percentage") +
  ylab("Splitting ratio") +
  ggtitle("B") +
  theme_classic() +
  theme(legend.position = "none") +
  guides(fill=guide_legend(title="")) +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  geom_hline(yintercept= 1, linetype="dashed") +
  theme(plot.title = element_text(size=16, face = "bold")) ; splitting_boxplot

# PERCENTAGE SINGLETONS

singletons = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/percentage_single_sample_GMYC_species_amalgamated.csv", header = T, check.names = F)
singletons_melt = reshape2::melt(singletons)
colnames(singletons_melt) = c("file_name", "data_perc", "perc_singletons")

singletons_boxplot = ggplot(data = singletons_melt, aes(x = data_perc, y = perc_singletons)) +
  geom_boxplot(fill = "#B5D5FA") +
  stat_summary(fun = mean, geom="point", shape=17, size=4, color="black", position=position_dodge(0.77)) +
  xlab("Data Percentage") +
  ylab("Percentage singletons") +
  ggtitle("C") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(plot.title = element_text(size=16, face = "bold")) ; singletons_boxplot

# SPLITTING PER SPECIES

oversplits_per_group = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/mean_oversplits_per_group_full_100.csv", header = T, check.names = F)

splitting_bar = ggplot(data = oversplits_per_group, aes(x = predef_unique, y = Freq)) +
  geom_bar(stat = "summary", fill = "lightgrey", col = "black") +
  xlab("Morphospecies") +
  ylab("Mean splitting ratio per morphospecies \n group, on the full dataset (100%)") +
  theme_classic() +
  ggtitle("D") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  geom_hline(yintercept= 1, linetype="dashed") ; splitting_bar


# PERCENTAGE MATCHES

percentage_matches = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/percentage_match_amalgamated.csv", check.names = F)
percentage_matches_m = reshape2::melt(percentage_matches)

percentage_matches_excl_s = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_COI/percentage_match_excl_singles_amalgamated.csv", check.names = F)
percentage_matches_excl_s_melt = reshape2::melt(percentage_matches_excl_s)

percentage_matches_m$ex_sing = percentage_matches_excl_s_melt$value
colnames(percentage_matches_m) = c("file_name", "data_perc", "inc_singletons", "ex_singletons")
percentage_matches_combo = reshape2::melt(percentage_matches_m)

perc_match_boxplot = ggplot(data = percentage_matches_combo, aes(x = data_perc, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  stat_summary(fun = mean, aes(group = variable), geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  scale_fill_manual(values=c("#2D86F0", "#F02D5A"), name = "", labels = c("+ singletons", "- singletons")) +
  xlab("Data Percentage") +
  ylab("Percentage match") +
  ggtitle("E") +
  theme_classic() +
  theme(legend.position = "none") +
  guides(fill=guide_legend(title="")) +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  ylim(0,100) +
  theme(plot.title = element_text(size=16, face = "bold")) ; perc_match_boxplot 

# COMBINE THE PLOTS

combo_plots_COI = gridExtra::grid.arrange(clust_accum,
                                          ent_accum,
                                          splitting_boxplot, 
                                          singletons_boxplot, 
                                          splitting_bar, 
                                          perc_match_boxplot)

# SAVE THE OUTPUT

ggsave(plot = combo_plots_COI, width = 25, height = 35, dpi = 350, filename = "combo_plots_COI.svg", units = "cm")

#################
#################
# 28S RESULTS
################
#################

data_sheet_28S = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/28S_info_sheet.csv")
length( unique(data_sheet_28S$morphospecies) )

count(data_sheet_28S, "morphospecies")

# CLUSTERS

clusters = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/clusters_amalgamated.csv", header = T, check.names = F)
clusters_melt = reshape2::melt(clusters)
colnames(clusters_melt) = c("file_name", "data_perc", "clusters")

# ENTITIES

entities = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/entities_amalgamated.csv", header = T, check.names = F)
entities_melt = reshape2::melt(entities)
colnames(entities_melt) = c("file_name", "data_perc", "entities")

clusters_melt$entities = entities_melt$entities 

clusts_ents = reshape::melt(clusters_melt)

# BOXPLOT

clust_ent_boxplot = ggplot(data = clusts_ents, aes(x = data_perc, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  stat_summary(fun = mean, aes(group = variable), geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  scale_fill_manual(values=c("#1EAD4B", "#AC72D9"), name = "", labels = c("Clusters", "Entities")) +
  xlab("Data Percentage") +
  ylab("Number of clusters or entities") +
  ggtitle("A") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  geom_hline(yintercept= 5, linetype="dashed") +
  geom_hline(yintercept= 11, linetype="dashed") ; clust_ent_boxplot
#guides(fill=guide_legend(title=""))

# PLOT JUST CLUSTERS, WITH AN ACCUMULATION CURVE

clust_accum = ggplot(data = clusters_melt, aes(x = as.numeric( data_perc ), y = clusters)) + 
  geom_point() +
  geom_smooth(span = 3, col = "blue", fill = "lightblue", alpha = 0.5) +
  #geom_boxplot(fill = "#1EAD4B") +
  #stat_summary(fun = mean, geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  xlab("Data Percentage") +
  ylab("Number of clusters") +
  ggtitle("A1") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  ylim(0,20) +
  geom_hline(yintercept= 14, linetype="dashed") +
  scale_x_continuous(breaks = as.numeric(clusters_melt$data_perc), labels = clusters_melt$data_perc) ; clust_accum

# PLOT JUST ENTITIES, WITH AN ACCUMULATION CURVE

ent_accum = ggplot(data = clusters_melt, aes(x = as.numeric( data_perc ), y = entities)) + 
  geom_point() +
  geom_smooth(span = 3, col = "blue", fill = "lightblue", alpha = 0.5) +
  #geom_boxplot(fill = "#1EAD4B") +
  #stat_summary(fun = mean, geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  xlab("Data Percentage") +
  ylab("Number of entities") +
  ggtitle("A2") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  ylim(0,20) +
  geom_hline(yintercept= 14, linetype="dashed") ; ent_accum

# SPLITTING RATIOS

splitting = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/oversplitting_ratio_amalgamated.csv", header = T, check.names = F)
splitting_melt = reshape2::melt(splitting)
colnames(splitting_melt) = c("file_name", "data_perc", "splitting_ratio")

splitting_excl_sing = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/oversplitting_excl_singles_amalgamated.csv", header = T, check.names = F)
splitting_excl_sing_melt = reshape2::melt(splitting_excl_sing)
colnames(splitting_excl_sing_melt) = c("file_name", "data_perc", "splitting_ratio_excl")

splitting_melt$splitting_excl = splitting_excl_sing_melt$splitting_ratio_excl

splitting_combo = reshape::melt( splitting_melt )

# BOXPLOT

splitting_boxplot = ggplot(data = splitting_combo, aes(x = data_perc, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  stat_summary(fun = mean, aes(group = variable), geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  scale_fill_manual(values=c("#2D86F0", "#F02D5A"), name = "", labels = c("+ singletons", "- singletons")) +
  xlab("Data Percentage") +
  ylab("Splitting ratio") +
  ggtitle("B") +
  theme_classic() +
  theme(legend.position = "none") +
  guides(fill=guide_legend(title="")) +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  geom_hline(yintercept= 1, linetype="dashed") +
  theme(plot.title = element_text(size=16, face = "bold")) ; splitting_boxplot

# PERCENTAGE SINGLETONS

singletons = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/percentage_single_sample_GMYC_species_amalgamated.csv", header = T, check.names = F)
singletons_melt = reshape2::melt(singletons)
colnames(singletons_melt) = c("file_name", "data_perc", "perc_singletons")

singletons_boxplot = ggplot(data = singletons_melt, aes(x = data_perc, y = perc_singletons)) +
  geom_boxplot(fill = "#B5D5FA") +
  stat_summary(fun = mean, geom="point", shape=17, size=4, color="black", position=position_dodge(0.77)) +
  xlab("Data Percentage") +
  ylab("Percentage singletons") +
  ggtitle("C") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(plot.title = element_text(size=16, face = "bold")) ; singletons_boxplot

# SPLITTING PER SPECIES

oversplits_per_group = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/mean_oversplits_per_group_full_100.csv", header = T, check.names = F)

splitting_bar = ggplot(data = oversplits_per_group, aes(x = predef_unique, y = Freq)) +
  geom_bar(stat = "summary", fill = "lightgrey", col = "black") +
  xlab("Morphospecies") +
  ylab("Mean splitting ratio per morphospecies \n group, on the full dataset (100%)") +
  theme_classic() +
  ggtitle("D") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 4, by = 0.5)) +
  geom_hline(yintercept= 1, linetype="dashed") ;splitting_bar


# PERCENTAGE MATCHES

percentage_matches = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/percentage_match_amalgamated.csv", check.names = F)
percentage_matches_m = reshape2::melt(percentage_matches)

percentage_matches_excl_s = read.csv("spede_sampler/spede_tetramesa_plotting/RESULTS_28S/percentage_match_excl_singles_amalgamated.csv", check.names = F)
percentage_matches_excl_s_melt = reshape2::melt(percentage_matches_excl_s)

percentage_matches_m$ex_sing = percentage_matches_excl_s_melt$value
colnames(percentage_matches_m) = c("file_name", "data_perc", "inc_singletons", "ex_singletons")
percentage_matches_combo = reshape2::melt(percentage_matches_m)

perc_match_boxplot = ggplot(data = percentage_matches_combo, aes(x = data_perc, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  stat_summary(fun = mean, aes(group = variable), geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  scale_fill_manual(values=c("#2D86F0", "#F02D5A"), name = "", labels = c("+ singletons", "- singletons")) +
  xlab("Data Percentage") +
  ylab("Percentage match") +
  ggtitle("E") +
  theme_classic() +
  theme(legend.position = "none") +
  guides(fill=guide_legend(title="")) +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  ylim(0,100) +
  theme(plot.title = element_text(size=16, face = "bold")) ; perc_match_boxplot 

# COMBINE THE PLOTS

combo_plots_28S = gridExtra::grid.arrange(clust_accum,
                                          ent_accum,
                                          splitting_boxplot, 
                                          singletons_boxplot, 
                                          splitting_bar, 
                                          perc_match_boxplot)

# SAVE THE OUTPUT

ggsave(plot = combo_plots_28S, width = 25, height = 35, dpi = 350, filename = "combo_plots_28S.svg", units = "cm")


##############################################################################################
##############################################################################################
##############################################################################################
# PART 8: MAP PLOTTING
##############################################################################################
##############################################################################################
##############################################################################################
# (1)
##############################################################################################
# Plot the global species description effort of the Tetramesa
# Thanks to Guy Sutton for the code that has been adapted below, and for the input file
##############################################################################################

theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill = NA),
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
                  legend.position = "none"))

# Load data
raw_data <- readxl::read_xlsx("maps/tetramesa_inventory_global.xlsx")

# Split the countries column into many columns (1 country per column)
countries_data <- splitstackshape::cSplit(raw_data, 
                         'countries_recorded', 
                         sep = ",", 
                         fixed = FALSE)

# Now reshape into long format
countries_long <- countries_data %>%
  pivot_longer(
    # cols = which columns do we want to pivot/move
    cols = starts_with("countries_recorded_"),
    # names_to = new column name that the names of cols above will be
    # moved to. This effectively creates your categorical
    # factor levels
    names_to = "country",
    # values_to = new column where the row values of cols will be stored
    values_to = "country_present")


# Count number of known Tetramesa species per country
tetra_dist_sum <- countries_long %>%
  drop_na(country_present) %>%
  summarise(no_rows = n()) 
tetra_dist_sum

tetra_dist_sum <- countries_long %>%
  drop_na(country_present) %>%
  dplyr::group_by(country_present) %>%
  dplyr::count()
tetra_dist_sum

# Load a world map 
# Load required libraries
library(sp)
library(sf)
library(raster)
library(rasterVis)
library(maptools)
library(dismo)
library(raster)
library(rgeos)
library(rJava)
library(rgdal)
library(geosphere)
library(scales)
library(maptools)
library(mapdata)
library(spThin)
library(ENMeval)
library(plyr)
library(grid)
library(gridSVG)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

# load up the world map
world <- ne_countries(scale = "medium", returnclass = "sf")

tetra_dist = tetra_dist_sum
tetra_dist$name = tetra_dist$country_present

tetra_dist = tetra_dist %>%
  dplyr::mutate(name = recode(name, USA = 'United States',
                                         PRC = "China",
                                         USSR = "Russia"))


world_data <- full_join(world, tetra_dist, by = "name")

world_data <- world_data %>%
  mutate(no_rows = as.numeric(n))

# Set the theme for the plot 
theme_opts<-list(theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_rect(fill = 'white', colour = NA),
                       plot.background = element_rect(),
                       axis.line = element_blank(),
                       axis.text.x = element_text(colour = "black"),
                       axis.text.y = element_text(colour = "black"),
                       axis.ticks = element_line(colour = "black"),
                       axis.title.x = element_text(colour = "black"),
                       axis.title.y = element_text(colour = "black"),
                       plot.title = element_text(colour = "black"),
                       panel.border = element_rect(fill = NA),
                       legend.key=element_blank()))

# Set 0 to NA
world_data$no_rows[world_data$no_rows == 0] <- NA

############################################################################
# PLOT THE WORLD MAP WITH THE NUMBER OF SPECIES DESCRIPTIONS
############################################################################

p <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = world_data, aes(fill = n)) +
  scale_fill_gradientn(colours = colorspace::heat_hcl(n=7, alpha = 0.5), na.value = "white") +
  labs(fill = "No. of species per country") + 
  coord_sf(xlim = c(-180, 180), 
           ylim = c(-90, 90), 
           crs = 4326, 
           expand = FALSE) +
  theme_opts +
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 1))

p + theme(legend.position = "bottom")

# colour the sampling effort
p + scale_fill_gradient(low="lightblue", high="red", na.value = "white") + theme(legend.position = "bottom")


############################################################################
############################################################################
# PLOT DISTRIBUTION MAPS OF GRASS SPECIES IN THE NATIVE AND INVADED RANGES
############################################################################
############################################################################

# (2)

# Load required packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  tidyverse,
  dismo,
  raster,
  here,
  corrplot,
  Hmisc,
  patchwork,
  ecospat,
  kuenm,
  gridSVG,
  gridExtra,
  grid,
  ENMeval,
  spThin,
  viridis,
  viridisLite,
  mapdata,
  maptools,
  scales,
  geosphere,
  rgdal,
  ggtext,
  rJava,
  rgeos,
  sp,
  sf,
  ggspatial,
  ecospat,
  rnaturalearth,
  rnaturalearthdata,
  #megaSDM,
  InformationValue,
  caret,
  spocc,
  scrubr
)

library(ggthemes)

# Change ggplot theme
theme_set(
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black",
                                  fill = NA),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = unit(c(2, 0, 0, 0),
                                                "mm")),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0),
                                                "mm")),
      legend.position = "none"
    )
)

# Set the theme for the maps
theme_opts <- list(
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    plot.background = element_rect(),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    plot.title = element_text(colour = "black"),
    panel.border = element_rect(fill = NA),
    legend.key = element_blank()
  )
)

################################################################################
# PLOT THJE DISTRIBUTION OF ERAGROSTIS CURVULA, FOR EXAMPLE, IN SOUTH AFRICA
################################################################################

# Get map of South Africa to project our model over
south_africa <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf") %>%
  dplyr::filter(name == "South Africa")


# List the countries we want records from 
gbifopts_rsa <- list(country = "ZA") 

# Extract records
df_rsa <- occ(query = "Eragrostis curvula", 
          from = c("gbif"),
          gbifopts = gbifopts_rsa,
          # Limit the maximum number of records
          limit = 10000)

# Process data 
df_comb_rsa <- occ2df(df_rsa)

################################################################################
# Plot base map
################################################################################

rsa_ecurvula = ggplot() +
  geom_sf(data = south_africa) +
  coord_sf(
    xlim = c(15, 33.5),
    ylim = c(-35,-21.5),
    crs = 4326,
    expand = FALSE
  ) +
  geom_point(data = df_comb_rsa, aes(x = longitude, y = latitude), alpha = 0.5, col = "darkgreen") + 
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.25, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Apply the theme for the map we defined above. 
  theme_opts +
  theme(legend.position = "right") +
  labs(
    title = "Native-range distribution",
    x = "Longitude",
    y = "Latitude"
  )

################################################################################
# PLOT THJE DISTRIBUTION OF ERAGROSTIS CURVULA IN AUSTRALIA, THE INVADED RANGE
################################################################################

# Get map of Australia to project our model over
australia <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf") %>%
  dplyr::filter(name == "Australia")


# List the countries we want records from 
gbifopts_au <- list(country = "AU") 

# Extract records
df_au <- occ(query = "Eragrostis curvula", 
          from = c("gbif"),
          gbifopts = gbifopts_au,
          # Limit the maximum number of records
          limit = 10000)

# Process data 
df_comb_au <- occ2df(df_au)

################################################################################
# Plot base map
################################################################################

# Plot base map
aus_ecurvula = ggplot() +
  geom_sf(data = australia) +
  coord_sf(
    xlim = c(109, 155),
    ylim = c(-45,-8),
    crs = 4326,
    expand = FALSE
  ) +
  geom_point(data = df_comb_au, aes(x = longitude, y = latitude), alpha = 0.5, col = "red") + 
  # Add scale bar to bottom-right of map
  annotation_scale(location = "bl", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.3, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Apply the theme for the map we defined above. 
  theme_opts +
  theme(legend.position = "right") +
  labs(
    title = "Invaded-range distribution",
    x = "Longitude",
    y = "Latitude"
  )

##################################################################################
# COMBINE MAPS
##################################################################################

map_combo = gridExtra::grid.arrange(rsa_ecurvula, aus_ecurvula, nrow = 1)

# Save image
ggsave("./curvula_invaded_range.png", plot = map_combo,
       dpi = 600,
       height = 10,
       width = 10)
