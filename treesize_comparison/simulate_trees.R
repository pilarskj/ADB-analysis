## ---------------------------
##
## Script for simulating phylogenetic trees with different size under the Age-Dependent Branching model
##
## Author: Julia Pilarski
## Email: julia.pilarski@bsse.ethz.ch
##
## Date Created: 2025-06-25
##
## ---------------------------

setwd("~/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison")

## load packages
library(parallel)
library(scTreeSim)
library(ape)

## settings
set.seed(1)
nsim = 10 # number of trees per size
ntips = c(10, 100, 1000, 5000) # tree sizes 

## set model parameters
shape = 5
lifetime = 10
deathprob = 0.1
rho = 0.1
scale = lifetime / shape

## simulate trees 
tree_data = data.frame(tree = character(), ntips = numeric(), origin = numeric(), treeHeight = numeric(), treeLength = numeric())

for (n in ntips) {
  trees = mclapply(c(1:nsim), function(i) {
    repeat { 
      tree = sim_adb_ntaxa_samp(ntaxa = n, a = scale, b = shape, d = deathprob, rho = rho)
      if (!is.null(tree)) break
    }
    tree = tree@phylo
    tree$tip.label = c(1:Ntip(tree))
    return(tree) }, mc.cores = 8, mc.silent = F)
  
  # save trees
  lapply(c(1:nsim), function(i) write.tree(trees[[i]], paste0('trees/tree_n', n, '_', i, '.newick')))
  
  # add tree information
  tree_data = rbind(tree_data, data.frame(tree = c(1:nsim), ntips = n, 
                                          origin = sapply(trees, function(t) t$origin),
                                          treeHeight = sapply(trees, function(t) max(node.depth.edgelength(t))),
                                          treeLength = sapply(trees, function(t) sum(t$edge.length))))
  print(paste0('n = ', n, ' completed!'))
}

# store ground-truth  information
write.csv(tree_data, file = 'tree_data.csv', quote = F, row.names = F)


