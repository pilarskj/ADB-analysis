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
df = read.csv("tree_likelihood.csv")
df_runtime = df %>% 
  select(-exact_logL, -approx_logL) %>%
  pivot_longer(cols = c("exact_time", "approx_time"), names_to = "method", values_to = "time") %>%
  mutate(ntips = factor(ntips, levels = ntips_ls), method = as.factor(sub("_time", "", method))) 

g1 = ggplot(df_runtime, aes(x = ntips, y = time, color = method)) +
  geom_point(size = 1) +
  stat_summary(aes(group = method), fun = "mean", geom = "line") +
  scale_color_manual(values = c("#009E73","#D55E00")) +
  geom_rect(xmin = 0.5, xmax = 4.5, ymin = 0, ymax = 300, fill = NA, color = "grey", linetype = "dashed", linewidth = 0.5) + # create inset
  labs(x = "Tree size", y = "Time for log-likelihood\ncomputation (ms)", color = NULL) 

g2 = ggplot(df_runtime, aes(x = ntips, y = time, color = method)) +
  geom_point(size = 1) +
  stat_summary(aes(group = method), fun = "mean", geom = "line") +
  scale_color_manual(values = c("#009E73","#D55E00")) +
  coord_cartesian(ylim = c(0,300)) +
  labs(x = "Tree size", y = "Time for log-likelihood\ncomputation (ms)", color = NULL) 
  
g1 + g2 + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
ggsave("likelihood_runtime.pdf", width = 8, height = 4)

## ---------------------------

## calculate error due to approximation
df_error = df %>% 
  mutate(abs_error = abs(exact_logL - approx_logL)) %>%
  mutate(rel_error = abs_error / ntips) %>%
  mutate(ntips = factor(ntips, levels = ntips_ls))

g1 = ggplot(df_error, aes(x = ntips, y = abs_error)) +
  geom_point(size = 1) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line") +
  labs(x = "Tree size", y = "Absolute error") 

g2 = ggplot(df_error, aes(x = ntips, y = rel_error)) +
  geom_point(size = 1) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line") +
  scale_y_continuous(limits = c(0, 0.05)) +
  labs(x = "Tree size", y = "Relative error") 

g1 + g2 + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
ggsave("likelihood_error.pdf", width = 8, height = 3.5)

## ---------------------------

## load inference logs
log_fs = lapply(ntips_ls, function(n) paste0('inference/inference_n', n, '_', c(1:10), '.log')) %>% unlist()
data_ls = lapply(log_fs, function(f) { remove_burn_ins(parse_beast_tracelog_file(f), burn_in_fraction = 0.1) })
names(data_ls) = lapply(log_fs, function(f) str_extract(sapply(str_split(f, '/'), tail, 1), "n[0-9]+_[0-9]+"))

## assess ESS
ess = lapply(data_ls, function(log) calc_esses(log, sample_interval = 1000)) %>% bind_rows(.id = "ID")
all(ess %>% select(-ID, -shape) >= 200)

## get runtime of MCMC chains
# helper function for converting elapsed time (D-HH:MM:SS) to seconds
elapsed_to_seconds <- function(elapsed) {
  if (str_detect(elapsed, "-")) {
    splits = str_split_1(elapsed, "-")
    days = as.numeric(splits[1])
    elapsed = splits[2]
  } else {
    days = 0
  }
  time_parts = str_split_1(elapsed, ":")
  seconds = days * 24*3600 + as.numeric(time_parts[1]) * 3600 + as.numeric(time_parts[2]) * 60 + as.numeric(time_parts[3])
  return(seconds)
}

# calculate proportion of samples needed for convergence
cfactor = sapply(data_ls, function(log) {
  for (i in c(2:nrow(log))) {
    current_ess = calc_ess(log$posterior[1:i], sample_interval = 1000)
    if (current_ess >= 200) {
      proportion = i / nrow(log)
      return(proportion)
    }
  }
})

# time recording from cluster
df_time = read.csv("job_duration.csv") %>%
  rowwise() %>%
  mutate(chainID = sapply(str_split(JobID, '_'), '[[', 2),
         runtime = elapsed_to_seconds(Elapsed)) %>%
  group_by(chainID) %>%
  summarise(runtime = sum(runtime)) %>%
  ungroup() %>%
  right_join(tree_data %>% rownames_to_column("chainID"), by = "chainID") %>%
  left_join(ess %>% select(ID, ESS = posterior), by = "ID") %>%
  mutate(cfactor = cfactor[ID]) %>% # simplefactor: 200 / ESS, alternatively, use subsampling
  mutate(ctime = runtime * cfactor) # time to convergence


ggplot(df_time, aes(x = ntips, y = ctime)) +
  geom_point(size = 1, alpha = 0.8) + #  position = position_jitter(0.1)
  scale_y_continuous(limits = c(1, 25)*3600, breaks = seq(5, 25, by = 5)*3600, labels = function(x) sprintf("%.0f h", x/3600)) + 
  stat_summary(aes(group = 1), fun = "mean", geom = "line") +
  labs(x = "Tree size", y = "Time to convergence (h)") 
ggsave("inference_runtime.pdf", width = 4.5, height = 3.5)


## evaluate parameter estimates
parameters = c('shape', 'lifetime', 'deathprob', 'rho')
data = data_ls %>% bind_rows(.id = "ID") %>%
  left_join(tree_data) %>% 
  pivot_longer(cols = parameters, names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = parameters))

# plot with true parameters as reference
g1 = ggplot(data_long %>% filter(parameter == "shape"), aes(x = interaction(tree, ntips), y = value, fill = ntips)) +
  geom_violin_discrete(show.legend = F) +
  geom_hline(yintercept = 5, linetype = 'dashed') +
  scale_fill_manual(values = palette) +
  coord_cartesian(ylim = c(0,50)) +
  labs(x = "Tree size", y = expression(paste("Shape ", italic(k)))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

g2 = ggplot(data_long %>% filter(parameter == "lifetime"), aes(x = ntips, y = value, fill = ntips, color = ntips, group = interaction(ntips, tree))) +
  geom_violin(position = "dodge", linewidth = 0.1, show.legend = F, scale = "width") + 
  geom_hline(yintercept = 10, linetype = 'dashed', alpha = 0.8) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) + 
  labs(x = "Tree size", y = expression(paste("Mean lifetime ", italic(l)))) 

g3 = ggplot(data_long %>% filter(parameter == "deathprob"), aes(x = ntips, y = value, fill = ntips, color = ntips, group = interaction(ntips, tree))) +
  geom_violin(position = "dodge", linewidth = 0.1, show.legend = F, scale = "width") + 
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.8) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) + 
  labs(x = "Tree size", y = expression(paste("Death probability ", italic(d)))) 

g4 = ggplot(data_long %>% filter(parameter == "rho"), aes(x = ntips, y = value, fill = ntips, color = ntips, group = interaction(ntips, tree))) +
  geom_violin(position = "dodge", linewidth = 0.1, show.legend = F, scale = "width") + 
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.8) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) + 
  labs(x = "Tree size", y = expression(paste("Sampling probability ", rho))) 

g1 + g2 + g3 + g4 + plot_layout(ncol = 2, nrow = 2, axes = "collect_x")
ggsave("inference_results.pdf", width = 10, height = 8)
