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

r1 = ggplot(df_runtime, aes(x = ntips, y = time, color = method)) +
  geom_point(size = 1) +
  stat_summary(aes(group = method), fun = "mean", geom = "line") +
  scale_color_manual(values = c("#009E73","#D55E00")) +
  geom_rect(xmin = 0.5, xmax = 4.5, ymin = 0, ymax = 300, fill = NA, color = "grey", linetype = "dashed", linewidth = 0.5) + # create inset
  labs(title = "B", x = "Tree size", y = "Time for log-likelihood\ncomputation (ms)", color = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.25,0.8),
        plot.title = element_text(face= "bold"), plot.title.position = "plot")

r2 = ggplot(df_runtime, aes(x = ntips, y = time, color = method)) +
  geom_point(size = 1) +
  stat_summary(aes(group = method), fun = "mean", geom = "line") +
  scale_color_manual(values = c("#009E73","#D55E00")) +
  coord_cartesian(ylim = c(0,300)) +
  labs(x = "Tree size", y = "Time for log-likelihood\ncomputation (ms)", color = NULL) +
  theme(legend.position = "none")


## ---------------------------

## calculate error due to approximation
df_error = df %>% 
  mutate(abs_error = abs(exact_logL - approx_logL)) %>%
  mutate(rel_error = abs_error / ntips) %>%
  mutate(ntips = factor(ntips, levels = ntips_ls))

e1 = ggplot(df_error, aes(x = ntips, y = abs_error)) +
  geom_point(size = 1) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line") +
  labs(x = "Tree size", y = "Absolute error") 

e2 = ggplot(df_error, aes(x = ntips, y = rel_error)) +
  geom_point(size = 1) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line") +
  scale_y_continuous(limits = c(0, 0.05)) +
  labs(x = "Tree size", y = "Relative error") 

e1 + e2 + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
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

# add to runtime plot
r3 = ggplot(df_time, aes(x = ntips, y = ctime)) +
  geom_point(size = 1, alpha = 0.8) + #  position = position_jitter(0.1)
  scale_y_continuous(limits = c(1, 25)*3600, breaks = seq(5, 25, by = 5)*3600, labels = function(x) sprintf("%.0f h", x/3600)) + 
  stat_summary(aes(group = 1), fun = "mean", geom = "line") +
  labs(title = "C", x = "Tree size", y = "Time to convergence (h)") +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot")

r1 + r2 + r3 + plot_layout(axes = "collect_y")
ggsave("runtime.pdf", width = 10, height = 4)


## evaluate parameter estimates
parameters = c('shape', 'lifetime', 'deathprob', 'rho')
data = data_ls %>% bind_rows(.id = "ID") %>%
  left_join(tree_data) %>% 
  pivot_longer(cols = parameters, names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = parameters))

# plot with true parameters as reference
g1 = ggplot(data %>% filter(parameter == "shape"), aes(x = interaction(tree, ntips), y = value, fill = ntips)) +
  geom_violin_discrete(show.legend = F) +
  geom_hline(yintercept = 5, linetype = 'dashed') +
  scale_fill_manual(values = palette) +
  coord_cartesian(ylim = c(0,50)) +
  labs(x = "Tree size", y = expression(paste("Shape ", italic(k)))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

g2 = ggplot(data %>% filter(parameter == "lifetime"), aes(x = ntips, y = value, fill = ntips, color = ntips, group = interaction(ntips, tree))) +
  geom_violin(position = "dodge", linewidth = 0.1, show.legend = F, scale = "width") + 
  geom_hline(yintercept = 10, linetype = 'dashed', alpha = 0.8) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) + 
  labs(x = "Tree size", y = expression(paste("Mean lifetime ", italic(l)))) 

g3 = ggplot(data %>% filter(parameter == "deathprob"), aes(x = ntips, y = value, fill = ntips, color = ntips, group = interaction(ntips, tree))) +
  geom_violin(position = "dodge", linewidth = 0.1, show.legend = F, scale = "width") + 
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.8) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) + 
  labs(x = "Tree size", y = expression(paste("Death probability ", italic(d)))) 

g4 = ggplot(data %>% filter(parameter == "rho"), aes(x = ntips, y = value, fill = ntips, color = ntips, group = interaction(ntips, tree))) +
  geom_violin(position = "dodge", linewidth = 0.1, show.legend = F, scale = "width") + 
  geom_hline(yintercept = 0.1, linetype = 'dashed', alpha = 0.8) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) + 
  labs(x = "Tree size", y = expression(paste("Sampling probability ", rho))) 

g1 + ggtitle("A") + g2 + g3 + g4 + plot_layout(ncol = 2, nrow = 2, axes = "collect_x") &
  theme(plot.title = element_text(face= "bold"), plot.title.position = "plot")
ggsave("inference_results.pdf", width = 10, height = 6)


## check coverage
# load get_sumstats_from_log from ADB-analysis/validation/evaluate_inference.R
truth = data.frame(parameter = factor(parameters, levels = parameters), true = c(5, 10, 0.1, 0.1))
test = lapply(data_ls, function(log) { get_sumstats_from_log(log, parameters) }) %>% 
  bind_rows(.id = "ID") %>%
  left_join(truth, by = "parameter") %>%
  left_join(tree_data, by = "ID") %>%
  mutate(correct = ifelse(true >= lower & true <= upper, T, F),
         parameter = factor(parameter, levels = parameters))

ggplot(test, aes(x = ntips, y = median, group = ID, color = correct)) +
  geom_point(size = 0.7, show.legend = F, position = position_dodge(0.9)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.5, show.legend = F, position = position_dodge(0.9)) +
  geom_hline(data = truth, aes(yintercept = true), alpha = 0.25) + 
  scale_color_manual(values = c("FALSE" = "#F8766D", "TRUE" = "#009E73")) + 
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(x = "Tree size", y = "Posterior median\nwith 95% HPD interval") +
  facet_wrap(~parameter, scales = "free", labeller = as_labeller(
    c("shape" = "'Shape'~italic(k)", "lifetime" = "'Mean lifetime'~italic(l)", 
      "deathprob" = "'Death probability'~italic(d)", "rho" = "'Sampling probability'~italic(rho)"), 
    label_parsed)) +
  theme(axis.text.x = element_blank())
ggsave("inference_coverage.pdf", width = 10, height = 6)

## ---------------------------

## demonstrate on likelihood curves the systematic bias in inferring death probabilities from large trees
test_death = test %>% 
  filter(ntips == 5000, parameter == "deathprob") %>%
  select(tree, lower, median, upper, correct) %>%
  filter(tree %in% c(2:4))

df_lik = read.csv("tree_likelihood_death.csv") %>%
  pivot_longer(cols = c("approx_logL", "exact_logL"), names_to = "method", values_to = "logL") %>%
  mutate(method = sub("_logL", "", method)) %>%
  filter(tree %in% c(2:4)) %>%
  distinct()

# view maxima
df_lik %>% group_by(tree, method) %>% filter(logL == max(logL))

# plot
ribbon_label = "posterior median with 95% HPD interval"
ggplot() + 
  geom_line(data = df_lik, aes(x = deathprob, y = logL, linetype = method)) + # %>% filter(deathprob <= 0.25)
  geom_vline(xintercept = 0.1, alpha = 0.3) +
  geom_vline(data = test_death, aes(xintercept = median, color = correct), linetype = "dotted", alpha = 0.5) +
  geom_rect(data = test_death, aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf, fill = correct), alpha = 0.1) + 
  scale_linetype_manual(values = c("approx" = "solid", "exact" = "dashed")) + 
  scale_color_manual(values = c("FALSE" = "#F8766D", "TRUE" = "#009E73"), labels = ribbon_label) + 
  scale_fill_manual(values = c("FALSE" = "#F8766D", "TRUE" = "#009E73"), labels = ribbon_label) + 
  facet_wrap(~tree, nrow = 1) +
  labs(x = expression(paste("Death probability ", italic(d))), y = "Log Lik", linetype = NULL, color = NULL, fill = NULL) +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "bottom") +
  guides(linetype = guide_legend(order = 1))
ggsave("approximation_error.pdf", width = 8, height = 3.5)

