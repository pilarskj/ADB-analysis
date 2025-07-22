## ---------------------------
##
## Script for simulating phylogenetic trees under the Age-Dependent Branching model for validation
##
## Author: Julia Pilarski
## Email: julia.pilarski@bsse.ethz.ch
##
## Date Created: 2025-05-30
##
## ---------------------------

setwd("~/Projects/P1_AgeDependentTrees/ADB-analysis/validation")

## install packages
library(devtools)
install("~/Projects/scTreeSim") ## tree simulator cloned from https://github.com/sccevo/scTreeSim.git
library(scTreeSim)
library(ape)

## settings
set.seed(2)
nsim = 100 # number of trees

## draw model parameters
shape = sample(1:100, size = nsim, replace = TRUE)
lifetime = rlnorm(n = nsim, meanlog = 2, sdlog = 1)
deathprob = rexp(n = nsim, rate = 10)
rho = rbeta(n = nsim, shape1 = 2, shape = 5)
scale = lifetime / shape
ntips = 100

## simulate trees
trees = lapply(c(1:nsim), function(i) {
  repeat { 
    tree = sim_adb_ntaxa_samp(ntaxa = ntips, a = scale[i], b = shape[i], d = deathprob[i], rho = rho[i])
    if (!is.null(tree)) break
  }
  tree = tree@phylo
  tree$tip.label = c(1:Ntip(tree))
  return(tree) }) 
lapply(c(1:nsim), function(i) write.tree(trees[[i]], paste0('trees/tree_', i, '.newick')))

## store ground-truth information
tree_data = data.frame(
  tree = c(1:nsim), 
  shape = shape,
  lifetime = lifetime,
  deathprob = deathprob,
  rho = rho,
  origin = sapply(trees, function(t) t$origin),
  treeHeight = sapply(trees, function(t) max(node.depth.edgelength(t))),
  treeLength = sapply(trees, function(t) sum(t$edge.length)))
write.csv(tree_data, file = 'tree_data.csv', quote = F, row.names = F)
