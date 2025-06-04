## ---------------------------
##
## Script for estimating the mean lifetime and death probability from complete intMEMOIR trees
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
library(treeman)

## ---------------------------

## get ground truth trees from https://github.com/seidels/tidetree-material/tree/main/Fig3/dat
# table belongs to the original intMEMOIR publication (https://doi.org/10.1126/science.abb3099)
# downloaded from https://data.caltech.edu/records/1444 (Fig 3/Fig3DandH/)
data = read.csv('groundTruthNewick.csv', sep = "\t") %>% filter(newick != "N/A")

# read in full trees
trees = lapply(data$newick, function(t) {
  tree = read.tree(text = t) 
  return(tree) })
names(trees) = sub("_data", "", data$file.name)

# count cell division & cell death events
death_count = 0
div_count = 0
for (tree in trees) {
  div_count = div_count + tree$Nnode # internal nodes represent cell divisions
  node_depth = sapply(strsplit(tree$tip.label, '_'), '[[', 1) %>% as.numeric() # first part of tip label records movie frame at death/ sampling
  node_depth = 216 - node_depth # compute backwards in time
  death_count = death_count + sum(node_depth > 0) # count tips which died before present
} # death_count 216, div_count 2039

# estimate death probability
d_emp = death_count / (death_count + div_count) # 0.09578714

# process trees in the same way as phylogenies, cf. https://github.com/seidels/tidetree-material/tree/main/Fig3/1_preprocessing 
# transform frames to hours (one frame stands for 15 min)
scale_by = 15/60
trees = lapply(trees, function(tree) {
  tree$edge.length = tree$edge.length * scale_by
  return(tree)
})

# each internal node represents a 15 min (0.25 h) interval 
# add 7.5 min to incoming and outgoing branches, effectively placing the true node height in the center of this time interval
change_branch_lengths <- function(childNode, tree){
  currentBrachLength = childNode['spn']
  if(childNode['tip']) {
    newBranchLength = currentBrachLength + 0.125
  } else {
    newBranchLength = currentBrachLength + 0.25
  }
  tree = setNdSpn(tree, id = node['id'], val = newBranchLength)
  return(tree)
}

for (i in c(1:length(trees))) {
  tree = as(trees[[i]], 'TreeMan')
  for (nodeID in c(tree['nds'], tree['tips'])) {
    node = tree[[nodeID]]
    if (node['root']) {
      next()
    } else {
      tree = change_branch_lengths(node, tree)
    }
  }
  trees[[i]] = as(tree, 'phylo')
}

# get all cell lifetimes and estimate mean
lifetimes = lapply(trees, function(tree) {
  ix = which(sapply(strsplit(tree$tip.label, '_'), '[[', 1) == "216") # exclude external branches (pruned cell lifetimes due to sampling)
  lifetimes = tree$edge.length[!(tree$edge[,2] %in% ix)] 
  return(lifetimes)
}) %>% unlist()
summary(lifetimes)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.125  10.750  12.000  12.168  13.500  31.250 
l_emp = mean(lifetimes) # 12.16842

# store estimates 
write.csv(data.frame(process = "empirical", shape = NA, lifetime = l_emp, deathprob = d_emp), 
          "estimates.csv", row.names = F, quote = F)
