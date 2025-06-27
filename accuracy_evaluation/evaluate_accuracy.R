## ---------------------------
##
## Script for assessing log-likelihood curcves under ADB
##
## Author: Julia Pilarski
## Email: julia.pilarski@bsse.ethz.ch
##
## Date Created: 2025-06-26
##
## ---------------------------

setwd("~/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation")

library(TreeSim)
library(scTreeSim)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
theme_set(theme_classic(base_size = 14) + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))) 

## ---------------------------

## simulate BD process with TreeSim
shape = 1
scale = 5 # equal to lifetime
d = 0.1 # death probability
rho = 0.1 # sampling probability
lambda = (1 - d) / scale
mu = d / scale
set.seed(1)
tree = sim.bd.taxa.age(n = 100, numbsim = 1, lambda = lambda, mu = mu, frac = rho, age = 50)[[1]] 
tree$tip.label = c(1:Ntip(tree))
# write.tree(tree, 'tree_bd.newick') # this tree is also saved in https://github.com/pilarskj/ADB/test_data

## ---------------------------

## compare log-likelihood curves between ADB and BDMM-Prime (calculated in BEAST2, see https://github.com/pilarskj/ADB/blob/main/test/test/adbp/BDTest.java)

# load log-likelihood values (originally from ADB/test_data)
df = read.csv('loglikTreeBD.csv') 
true = data.frame(parameter = c('lifetime', 'deathprob', 'rho'), value = c(5, 0.1, 0.1))

methods = c('adb_approx' = 'ADB (approx)', 'adb_exact' = 'ADB (exact)', 'bdmm' = 'BDMM-Prime')
params = c("lifetime" = "'Mean lifetime'~italic(l)",
           "deathprob" = "'Death probability'~italic(d)",
           "rho" = "'Sampling probability'~rho")
df = df %>% mutate(method = factor(method, levels = names(methods)),
                   parameter = factor(parameter, levels = names(params)))
true = true %>% mutate(parameter = factor(parameter, levels = names(params)))

# plot log-likelihood match
ggplot(df, aes(x = value, y = logL, color = method, linetype = method)) + 
  geom_line(linewidth = 1) +
  scale_linetype_manual(labels = methods, values = c("solid", "dashed", "dotted")) +
  scale_colour_manual(labels = methods, values = c("lightblue", "black", "red")) +
  scale_x_continuous(breaks = pretty_breaks()) +
  geom_vline(data = true, aes(xintercept = value), linetype = 'dashed', color = 'darkgrey') +
  facet_wrap(~parameter, scales = 'free', labeller = as_labeller(params, default = label_parsed), strip.position = "bottom") +
  labs(x = NULL, y = "Log Lik", col = "", linetype = "") +
  theme(strip.placement = "outside", legend.position = "bottom", strip.background = element_rect(colour = NA))
ggsave('loglik_BDMM_ADB.pdf', width = 8, height = 3.5)
