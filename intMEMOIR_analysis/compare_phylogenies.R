## ---------------------------
##
## Script for simulating phylogenies under inferred parameters and comparing them to observed lineage trees
##
## Author: Julia Pilarski
## Email: julia.pilarski@bsse.ethz.ch
##
## Date Created: 2025-06-02
##
## ---------------------------

setwd("~/Projects/P1_AgeDependentTrees/ADB-analysis/intMEMOIR_analysis")

library(dplyr)
library(ape)
library(treestats)
library(scTreeSim) ## tree simulator cloned from https://github.com/sccevo/scTreeSim.git
library(ggplot2)
library(patchwork)
theme_set(theme_classic(base_size = 14))
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")

## ---------------------------

# Function for extracting internal branches
get_internal_branches <- function(tree) {
  return(tree$edge.length[tree$edge[,2] > Ntip(tree)])
}

## ---------------------------

## read ground-truth phylogenies
# extracted from DREAM challenge publication (https://doi.org/10.1016/j.cels.2021.05.008) 
# and pre-processed following https://github.com/seidels/tidetree-material/tree/main/Fig3/1_preprocessing
ids = list.files(path = 'trees', pattern = ".newick") %>% sub("_data.newick", "", .)
trees = sapply(ids, function(i) read.tree(file = paste0("trees/", i, "_data.newick")), simplify = F, USE.NAMES = T) 

## load information on cell phylogenies
# from https://github.com/seidels/tidetree-material/tree/main/Fig3/dat
groundTruth = read.csv('groundTruthNewick.csv', sep = "\t")
dreamChallenge = read.csv('dreamChallengeDat.tsv', sep = "\t")
tree_data = full_join(groundTruth, dreamChallenge, by = c("file.name" = "fileName")) %>% 
  mutate(tree_id = sub("_data", "", file.name), rho = round(nCells / n.cells, digits = 2)) %>%
  select(tree_id, rho, n_full = n.cells, n_sampled = nCells, tree_full = newick, tree_sampled = ground) %>%
  filter(tree_id %in% ids)
all(names(trees) == tree_data$tree_id)

## ---------------------------

## simulate phylogenies with median phylodynamic parameters inferred under BD or ADB process
set.seed(123)
origin = 54 # assume initial cell is born at the start of the experiment
rho_seq = tree_data$rho
estimates = read.csv('estimates.csv')

# ADB
k_adb = estimates %>% filter(process == "ADB") %>% pull(shape)
l_adb = estimates %>% filter(process == "ADB") %>% pull(lifetime)
d_adb = estimates %>% filter(process == "ADB") %>% pull(deathprob)

trees_sim_adb = lapply(rho_seq, function(rho) {
  repeat { 
    tree = sim_adb_origin_samp(origin_time = origin, a = l_adb/k_adb, b = k_adb, d = d_adb, rho = rho) 
    if (!is.null(tree)) break
  }
  return(tree@phylo) }) 

# BD
l_bd = estimates %>% filter(process == "BD") %>% pull(lifetime)
d_bd = estimates %>% filter(process == "BD") %>% pull(deathprob) 

trees_sim_bd = lapply(rho_seq, function(rho) {
  repeat { 
    tree = sim_adb_origin_samp(origin_time = origin, a = l_bd, b = 1, d = d_bd, rho = rho) 
    if (!is.null(tree)) break
  }
  return(tree@phylo) }) 

## ---------------------------

## collect branch lengths
branch_lengths = rbind(lapply(trees, function(tree) data.frame(branch_length = get_internal_branches(tree))) %>% bind_rows(.id = 'tree') %>% mutate(process = "Empirical"),
                       lapply(trees_sim_adb, function(tree) data.frame(branch_length = get_internal_branches(tree))) %>% bind_rows(.id = 'tree') %>% mutate(process = "ADB"),
                       lapply(trees_sim_bd, function(tree) data.frame(branch_length = get_internal_branches(tree))) %>% bind_rows(.id = 'tree') %>% mutate(process = "BD")) %>%
  mutate(process = factor(process, levels = c("Empirical", "BD", "ADB")))


## plot distributions
g1 = ggplot(branch_lengths, aes(x = branch_length, fill = process)) + 
  geom_density(alpha = 0.6, color = NA) +
  scale_fill_manual(values = palette) +
  xlim(c(0,50)) +
  labs(title = "B", x = "Internal branch length", y = "Density", fill = NULL) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        legend.position = "inside", legend.position.inside = c(0.8,0.8))


## compare balance of simulated and observed phylogenies 
# B1 index https://www.rdocumentation.org/packages/treestats/versions/1.70.5/topics/b1
b1 = rbind(lapply(trees, function(tree) data.frame(b1_index = b1(tree, normalization = "tips"))) %>% bind_rows(.id = 'tree') %>% mutate(process = "Empirical"),
           lapply(trees_sim_adb, function(tree) data.frame(b1_index = b1(tree, normalization = "tips"))) %>% bind_rows(.id = 'tree') %>% mutate(process = "ADB"),
           lapply(trees_sim_bd, function(tree) data.frame(b1_index = b1(tree, normalization = "tips"))) %>% bind_rows(.id = 'tree') %>% mutate(process = "BD")) %>%
  mutate(process = factor(process, levels = c("Empirical", "BD", "ADB")))
g2 = ggplot(b1, aes(x = process, y = b1_index, fill = process)) + 
  geom_boxplot(outlier.size = 0.5, show.legend = F) +
  scale_y_continuous(limits = c(0,1), breaks = pretty_breaks()) +
  scale_fill_manual(values = palette) +
  labs(title = "C", x = NULL, y = "B1 balance index") +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot")

g1 + g2
ggsave("comparison_phylogenies.pdf", width = 8, height = 4)
