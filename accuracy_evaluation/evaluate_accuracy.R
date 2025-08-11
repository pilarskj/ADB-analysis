## ---------------------------
##
## Script for assessing log-likelihood curves under ADB
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
library(tibble)
library(ggplot2)
library(scales)
library(patchwork)
theme_set(theme_classic(base_size = 14) + 
            theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                  strip.placement = "outside", strip.background = element_rect(colour = NA))) 

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
write.tree(tree, 'trees/treeBD.newick') # this tree is also saved in https://github.com/pilarskj/ADB/test_data

## compare log-likelihood curves between ADB and BDMM-Prime (calculated in BEAST2, see https://github.com/pilarskj/ADB/blob/main/test/test/adbp/BDTest.java)
methods = c('adb_approx' = 'ADB (approx)', 'adb_exact' = 'ADB (exact)', 'bdmm' = 'BDMM-Prime')
params = c("lifetime" = "'Mean lifetime'~italic(l)",
           "deathprob" = "'Death probability'~italic(d)",
           "rho" = "'Sampling probability'~rho")
true = data.frame(parameter = factor(c('lifetime', 'deathprob', 'rho'), levels = names(params)), 
                  value = c(5, 0.1, 0.1))

# load log-likelihood values (originally from ADB/test_data)
df = read.csv('loglik/loglikTreeBD.csv') %>% 
  mutate(method = factor(method, levels = names(methods)),
         parameter = factor(parameter, levels = names(params)))

# plot log-likelihood match
ggplot(df, aes(x = value, y = logL, color = method, linetype = method)) + 
  geom_line(linewidth = 1) +
  scale_linetype_manual(labels = methods, values = c("solid", "dashed", "dotted")) +
  scale_colour_manual(labels = methods, values = c("lightblue", "black", "red")) +
  scale_x_continuous(breaks = pretty_breaks()) +
  geom_vline(data = true, aes(xintercept = value), linetype = 'dashed', color = 'darkgrey') +
  facet_wrap(~parameter, scales = 'free', labeller = as_labeller(params, default = label_parsed), strip.position = "bottom") +
  labs(x = NULL, y = "Log Lik", col = "", linetype = "") +
  theme(legend.position = "bottom")
ggsave('loglik_BDMM_ADB.pdf', width = 8, height = 3.5)

## ---------------------------

## now investigate accuracy at low sampling 
# re-simulate BD process with low sampling probability
rho = 0.001 
set.seed(1)
tree = sim.bd.taxa.age(n = 100, numbsim = 1, lambda = lambda, mu = mu, frac = rho, age = 100)[[1]] 
tree$tip.label = c(1:Ntip(tree))
write.tree(tree, 'trees/treeBD_low_sampling.newick') 

# compare log-likelihood curves 
params = c("lifetime" = "'Mean lifetime'~italic(l)",
           "deathprob" = "'Death probability'~italic(d)")
true = data.frame(parameter = factor(c('lifetime', 'deathprob'), levels = names(params)), value = c(5, 0.1))
df = read.csv('loglik/loglikTreeBD_low_sampling.csv') %>%
  mutate(method = factor(method, levels = names(methods)),
         parameter = factor(parameter, levels = names(params)))

# plot ADB (exact and approx) vs. BD
g1 = ggplot(df %>% filter(!stepSizeP0 %in% c("2^16", "2^18")), aes(x = value, y = logL, color = method, linetype = method)) + 
  geom_line(linewidth = 1) +
  scale_linetype_manual(labels = methods, values = c("solid", "dashed", "dotted")) +
  scale_colour_manual(labels = methods, values = c("lightblue", "black", "red")) +
  scale_x_continuous(breaks = pretty_breaks()) +
  geom_vline(data = true, aes(xintercept = value), linetype = 'dashed', color = 'darkgrey') +
  facet_wrap(~parameter, scales = 'free', labeller = as_labeller(params, default = label_parsed)) +
  labs(x = NULL, y = "Log Lik", col = "", linetype = "") +
  theme(strip.text = element_blank(), axis.text.x = element_blank())

# plot ADB with increasing FFT step size for P0 vs. BD
g2 = ggplot(mapping = aes(x = value, y = logL)) + 
  geom_line(data = df %>% filter(method == "adb_approx"), aes(color = stepSizeP0), linewidth = 0.7) +
  scale_colour_manual(labels = c(expression(2^14), expression(2^16), expression(2^18)), values = c("#9ECAE1", "#2171B5", "#08306B")) +
  geom_line(data = df %>% filter(method == "bdmm"), aes(linetype = "BDMM-Prime"), color = "red", linewidth = 1) +#, linetype = "dotted") +
  scale_linetype_manual(values = c("BDMM-Prime" = "dotted")) +
  geom_vline(data = true, aes(xintercept = value), linetype = 'dashed', color = 'darkgrey') +
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_wrap(~parameter, scales = 'free', labeller = as_labeller(params, default = label_parsed), strip.position = "bottom") +
  labs(x = NULL, y = "Log Lik", col = "Step size\n(ADB approx)", linetype = "") + 
  theme(legend.title = element_text(size = 12))

g1 / g2 + plot_layout(axes = "collect")
ggsave('loglik_BDMM_ADB_low_sampling.pdf', width = 8, height = 6)

## ---------------------------

## investigate approximation error on likelihood curves for death probability
# simulate trees with different parameters
bs = c(1, 5, 100)
ds = c(0.01, 0.1, 0.25)
tree_data = expand.grid(shape = bs, deathprob = ds)
set.seed(1)
trees = lapply(c(1:nrow(tree_data)), function(i) {
  repeat { 
    tree = sim_adb_ntaxa_samp(ntaxa = 100, a = 10 / tree_data[i, "shape"], b = tree_data[i, "shape"], d = tree_data[i, "deathprob"], rho = 0.1)
    if (!is.null(tree)) break
  }
  tree = tree@phylo
  tree$tip.label = c(1:Ntip(tree))
  return(tree) }) 

# store trees with metadata
tree_data$newick = sapply(trees, write.tree)
tree_data$origin = sapply(trees, function(t) t$origin)
tree_data = rownames_to_column(tree_data, var = "tree")
write.table(tree_data, "trees/trees_grid_death.tsv", sep = "\t", row.names = F, quote = F)

# load likelihood values (calculated in BEAST2, see https://github.com/pilarskj/ADB/blob/main/test/test/adbp/GammaBranchingModelTest.java)
# tree_data = read.table("trees/trees_grid_death.tsv", sep = "\t", header = T) %>% mutate(tree = as.character(tree))
df = read.csv("loglik/loglik_grid_death.csv") %>%
  pivot_longer(cols = c("approx_logL", "exact_logL"), names_to = "method", values_to = "logL") %>%
  mutate(tree = as.character(tree), method = sub("_logL", "", method)) %>%
  left_join(tree_data %>% select(tree, shape, true_deathprob = deathprob))

ggplot(df, aes(x = deathprob, y = logL, linetype = method, color = method)) + 
  geom_line() +
  geom_vline(aes(xintercept = true_deathprob), linetype = 'dashed', color = 'darkgrey') +
  scale_linetype_manual(values = c("approx" = "solid", "exact" = "dashed")) + 
  scale_color_manual(values = c("approx" = "#9ECAE1", "exact" = "black")) + 
  facet_grid(rows = vars(shape), cols = vars(true_deathprob), 
             labeller = label_bquote(cols = italic(d):.(true_deathprob),
                                     rows = italic(k):.(shape))) +
  labs(x = expression(paste("Death probability ", italic(d))), y = "Log Lik", linetype = NULL, color = NULL) +
  theme(legend.position = "bottom")
ggsave('loglik_grid_death.pdf', width = 8, height = 8)

## ---------------------------

## check that for high sampling, lower step size for FFT is sufficient
# simulate trees with different parameters
bs = c(1, 5, 20, 50, 100)
rhos = c(0.1, 0.25, 0.5, 0.75, 1)
tree_data = expand.grid(shape = bs, rho = rhos)
set.seed(1)
trees = lapply(c(1:nrow(tree_data)), function(i) {
  repeat { 
    tree = sim_adb_ntaxa_samp(ntaxa = 100, a = 5 / tree_data[i, "shape"], b = tree_data[i, "shape"], d = 0.1, rho = tree_data[i, "rho"])
    if (!is.null(tree)) break
  }
  tree = tree@phylo
  tree$tip.label = c(1:Ntip(tree))
  return(tree) }) 

# store trees with metadata
tree_data$newick = sapply(trees, write.tree)
tree_data$origin = sapply(trees, function(t) t$origin)
tree_data = rownames_to_column(tree_data, var = "tree")
write.table(tree_data, "trees/trees_grid_sampling.tsv", sep="\t", row.names = F, quote = F)

# load likelihood values (calculated in BEAST2
# tree_data = read.table("trees/trees_grid_sampling.tsv", sep = "\t", header = T) %>% mutate(tree = as.character(tree))
df = read.csv("loglik/loglik_grid_sampling.csv") %>%
  pivot_longer(cols = contains("logL"), names_to = "step_size", values_to = "logL") %>%
  mutate(tree = as.character(tree), step_size = sub("_logL", "", step_size)) %>%
  left_join(tree_data %>% select(tree, shape, rho))

s1 = ggplot(df %>% filter(param == "lifetime"), aes(x = value, y = logL, linetype = step_size, color = step_size)) + # %>% filter(lifetime >= 4 & lifetime <= 6)
  geom_line() +
  geom_vline(xintercept = 5, linetype = 'dashed', color = 'darkgrey') +
  facet_grid(rows = vars(shape), cols = vars(rho), 
             labeller = label_bquote(cols = rho:.(rho),
                                     rows = italic(k):.(shape)), scales = "free_y") +
  scale_linetype_manual(labels = c(expression(2^10), expression(2^12), expression(2^14)), values = c("solid", "dashed", "dotted")) +
  scale_colour_manual(labels = c(expression(2^10), expression(2^12), expression(2^14)), values = c("#9ECAE1", "#2171B5", "#08306B")) +
  labs(x = expression(paste("Mean lifetime ", italic(l))), y = "Log Lik", linetype = NULL, color = NULL)

s2 = ggplot(df %>% filter(param == "deathprob"), aes(x = value, y = logL, linetype = step_size, color = step_size)) + # %>% filter(lifetime >= 4 & lifetime <= 6)
  geom_line() +
  geom_vline(xintercept = 0.1, linetype = 'dashed', color = 'darkgrey') +
  facet_grid(rows = vars(shape), cols = vars(rho), 
             labeller = label_bquote(cols = rho:.(rho),
                                     rows = italic(k):.(shape)), scales = "free_y") +
  scale_linetype_manual(labels = c(expression(2^10), expression(2^12), expression(2^14)), values = c("solid", "dashed", "dotted")) +
  scale_colour_manual(labels = c(expression(2^10), expression(2^12), expression(2^14)), values = c("#9ECAE1", "#2171B5", "#08306B")) +
  labs(x = expression(paste("Death probability ", italic(d))), y = "Log Lik", linetype = NULL, color = NULL) 

s1 / s2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave('loglik_grid_high_sampling.pdf', width = 10, height = 18)

# compare to BD
test = rbind(read.csv("loglik/loglikBDMM_high_sampling.csv") %>% mutate(tree = as.character(tree)),
             df %>% filter(tree %in% c("1", "6", "11")) %>% select(tree, param, value, logL, method = step_size)) %>%
  left_join(tree_data %>% select(tree, rho))

test_l = test %>% filter(param == "lifetime")
s1_bd = ggplot(mapping = aes(x = value, y = logL)) + 
  geom_line(data = test_l %>% filter(method != "bdmm"), aes(color = method), linewidth = 0.7) +
  scale_colour_manual(labels = c(expression(2^10), expression(2^12), expression(2^14)), values = c("#9ECAE1", "#2171B5", "#08306B")) +
  geom_line(data = test_l %>% filter(method == "bdmm"), aes(linetype = "BDMM-Prime"), color = "red", linewidth = 1) +
  scale_linetype_manual(values = c("BDMM-Prime" = "dotted")) +
  geom_vline(xintercept = 5, linetype = 'dashed', color = 'darkgrey') +
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_grid(rows = vars(rho), labeller = label_bquote(rows = rho:.(rho)), scales = "free_y") +
  labs(x = expression(paste("Mean lifetime ", italic(l))), y = "Log Lik", col = "Step size\n(ADB approx)", linetype = "") + 
  theme(legend.title = element_text(size = 12), strip.text = element_blank()) + 
  guides(colour = guide_legend(order = 1))

test_d = test %>% filter(param == "deathprob")
s2_bd = ggplot(mapping = aes(x = value, y = logL)) + 
  geom_line(data = test_d %>% filter(method != "bdmm"), aes(color = method), linewidth = 0.7) +
  scale_colour_manual(labels = c(expression(2^10), expression(2^12), expression(2^14)), values = c("#9ECAE1", "#2171B5", "#08306B")) +
  geom_line(data = test_d %>% filter(method == "bdmm"), aes(linetype = "BDMM-Prime"), color = "red", linewidth = 1) +
  scale_linetype_manual(values = c("BDMM-Prime" = "dotted")) +
  geom_vline(xintercept = 0.1, linetype = 'dashed', color = 'darkgrey') +
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_grid(rows = vars(rho), labeller = label_bquote(rows = rho:.(rho)), scales = "free_y") +
  labs(x = expression(paste("Death probability ", italic(d))), y = "Log Lik", col = "Step size\n(ADB approx)", linetype = "") + 
  theme(legend.title = element_text(size = 12)) + 
  guides(colour = guide_legend(order = 1))

s1_bd + s2_bd + plot_layout(guides = "collect", axes = "collect") 
ggsave('loglik_BDMM_ADB_high_sampling.pdf', width = 8.2, height = 8)
