## ---------------------------
##
## Script for comparing phylodynamic inference under the Age-Dependent Branching model (for different tree sizes)
##
## Author: Julia Pilarski
## Email: julia.pilarski@bsse.ethz.ch
##
## Date Created: 2025-06-26
##
## ---------------------------

setwd("~/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison")

library(tracerer)
library(HDInterval)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scales)
library(geomViolinDiscrete)
theme_set(theme_classic(base_size = 14) + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))) 
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")

## ---------------------------

## load tree information
ntips_ls = c(10, 100, 1000, 5000)
tree_data = read.csv("tree_data.csv") %>%
  mutate(ID = paste0('n', ntips, '_', tree), ntips = factor(ntips, levels = ntips_ls), tree = factor(tree, levels = c(1:10))) %>%
  select(ID, ntips, tree)

## ---------------------------

## compare runtime for exact and approximated log-likelihood calculations (per iteration)
## calculated in BEAST2 (https://github.com/pilarskj/ADB/blob/main/test/test/adbp/GammaBranchingModelTest.java)
df = read.csv("runtime.csv")
df = df %>% 
  pivot_longer(cols = c("exact_time", "approx_time"), names_to = "method", values_to = "time") %>%
  mutate(ntips = factor(ntips, levels = ntips_ls), method = as.factor(sub("_time", "", method))) 

g1 = ggplot(df, aes(x = ntips, y = time, color = method)) +
  geom_point(size = 1) +
  stat_summary(aes(group = method), fun = "median", geom = "line") +
  scale_color_manual(values = c("#009E73","#D55E00")) +
  geom_rect(xmin = 0.5, xmax = 4.5, ymin = 0, ymax = 300, fill = NA, color = "grey", linetype = "dashed", linewidth = 0.5) + # create inset
  labs(x = "Tree size (number of tips)", y = "Time for log-likelihood\ncomputation (ms)", color = NULL) 


g2 = ggplot(df, aes(x = ntips, y = time, color = method)) +
  geom_point(size = 1) +
  stat_summary(aes(group = method), fun = "median", geom = "line") +
  scale_color_manual(values = c("#009E73","#D55E00")) +
  coord_cartesian(ylim = c(0,300)) +
  labs(x = "Tree size (number of tips)", y = "Time for log-likelihood\ncomputation (ms)", color = NULL) 
  
g1 + g2 + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
ggsave("runtime_likelihood.pdf", width = 8, height = 4)

# optionally, evaluate also for MCMC chains the runtime until convergence (posterior ESS > 200) - g3

## ---------------------------

## load inferred parameters
log_fs = lapply(ntips_ls, function(n) paste0('inference/inference_n', n, '_', c(1:10), '.log')) %>% unlist()
parameters = c('shape', 'lifetime', 'deathprob', 'rho')

# get logs
data = lapply(log_fs, function(f) {
  log = remove_burn_ins(parse_beast_tracelog_file(f), burn_in_fraction = 0.1)
  log$ID = str_extract(sapply(str_split(f, '/'), tail, 1), "n[0-9]+_[0-9]+")
  return(log)
}) %>% 
  bind_rows() %>%
  left_join(tree_data)

# TODO
# plot estimates
ggplot(data, aes(x = ntips, y = deathprob)) +
  geom_violin(aes(fill = ntips),  position = 'identity', color = NA, alpha = 0.5, show.legend = F) + #group = interaction(ntips, tree)
  #geom_boxplot(width = 0.05, outliers = FALSE, fill = NA, linewidth = 0.3, show.legend = F) +
  geom_hline(yintercept = 0.1, linetype = 'dashed') +
  #coord_cartesian(ylim = c(0,50)) +
  scale_fill_manual(values = palette) 
  # geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.5, show.legend = F) +
  # geom_text(data = cov, mapping = aes(x = -Inf, y = -Inf, label = paste(cov, "%")), hjust = -6.5, vjust = -1, inherit.aes = F) + 
  # geom_line(aes(x = true, y = true), inherit.aes = F, alpha = 0.25) + 
  # scale_color_manual(values = c("FALSE" = "#F8766D", "TRUE" = "#009E73")) + 
  # scale_y_continuous(breaks = pretty_breaks()) +
  # labs(x = "True value", y = "Posterior median\nwith 95% HPD interval") +
  # facet_wrap(~parameter, scales = "free", labeller = as_labeller(
  #   c("shape" = "'Shape'~italic(k)", "lifetime" = "'Mean lifetime'~italic(l)", 
  #     "deathprob" = "'Death probability'~italic(d)", "rho" = "'Sampling probability'~italic(rho)"), 
  #   label_parsed))


