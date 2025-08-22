## ---------------------------
##
## Script for evaluating phylodynamic inference under BD and ADB from intMEMOIR data 
##
## Author: Julia Pilarski
## Email: julia.pilarski@bsse.ethz.ch
##
## Date Created: 2025-06-02
##
## ---------------------------

setwd("~/Projects/P1_AgeDependentTrees/ADB-analysis/intMEMOIR_analysis")

library(tracerer)
library(HDInterval)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(geomViolinDiscrete)
library(ggdist)
library(scales)
library(patchwork)
theme_set(theme_classic(base_size = 14))
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")

## ---------------------------

## set empirical estimates from ground-truth cell population trees (see estimate_parameters.R)
estimates = read.csv('estimates.csv')
d_emp = estimates %>% filter(process == "empirical") %>% pull(deathprob)
l_emp = estimates %>% filter(process == "empirical") %>% pull(lifetime)


## load inference logs and assure convergence
bdsky_log_fs = paste0("inference_bdsky_10clocks_20scars_fixedTrees_condRoot/chain", c(1:5), ".log")
bdsky_log = lapply(bdsky_log_fs, function(f) remove_burn_ins(parse_beast_tracelog_file(f), burn_in_fraction = 0.1)) %>% 
  bind_rows(.id = "chain") %>%
  select(-starts_with("treeLikelihood"))
all(calc_esses(bdsky_log %>% select(-chain), sample_interval = 1000) > 200, na.rm = T)

adb_log_fs = paste0("inference_adb_10clocks_20scars_fixedTrees_condRoot/chain", c(1:5), ".log") 
adb_log = lapply(adb_log_fs, function(f) remove_burn_ins(parse_beast_tracelog_file(f), burn_in_fraction = 0.1)) %>% 
  bind_rows(.id = "chain") %>%
  select(-starts_with("treeLikelihood")) 
all(calc_esses(adb_log %>% select(-chain), sample_interval = 1000) > 200, na.rm = T)


## evaluate phylodynamic parameters
# store median estimates
estimates = rbind(estimates, 
                  data.frame(process = "BD", shape = NA, lifetime = median(bdsky_log$lifetime), deathprob = median(bdsky_log$deathprob)),
                  data.frame(process = "ADB", shape = median(adb_log$shape), lifetime = median(adb_log$lifetime), deathprob = median(adb_log$deathprob)))
write.csv(estimates, "estimates.csv", row.names = F, quote = F)

# extract distributions                           
set.seed(1) # for priors 
data = rbind(adb_log %>% select(lifetime, deathprob, shape) %>% mutate(model = "ADB"), 
             bdsky_log %>% select(lifetime, deathprob) %>% mutate(shape = 1, model = "BD"),
             # add prior distributions as reference for plotting
             data.frame(lifetime = rlnorm(n = 5e4, meanlog = 2.5, sdlog = 1), 
                        deathprob = rbeta(n = 5e4, shape1 = 2, shape2 = 5), 
                        shape = round(rlnorm(n = 5e4, meanlog = 2, sdlog = 1), 0), model = "Prior")) %>%
  mutate(model = factor(model, levels = c("Prior", "BD", "ADB")))

# view summary 
for (param in c("lifetime", "deathprob", "shape")) {
  print(param)
  summary = data %>% 
    group_by(model) %>%
    summarise(lower = HDInterval::hdi(.data[[param]], credMass = 0.95)[1], 
              median = median(.data[[param]]),
              upper = HDInterval::hdi(.data[[param]], credMass = 0.95)[2])
  print(summary)
}

# plot estimates
g1 = ggplot(data, aes(x = model, y = shape, fill = model)) + 
  geomViolinDiscrete::geom_violin_discrete(scale = "width", color = NA, show.legend = F) +
  stat_pointinterval(data = data %>% filter(model == "ADB"), mapping = aes(x = model, y = shape, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, show.legend = F) +
  scale_fill_manual(values = palette) + 
  scale_y_continuous(limits = c(0,50), breaks = seq(10,50,10)) +
  labs(x = NULL, y = expression(paste("Shape ", italic(k))), fill = NULL) +
  theme(axis.text.x = element_blank())

g2 = ggplot(data, aes(x = model, y = lifetime, fill = model)) +
  geom_violin(color = NA) + 
  stat_pointinterval(data = data %>% filter(model != "Prior"), mapping = aes(x = model, y = lifetime, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, show.legend = F) +
  geom_hline(yintercept = l_emp, linetype = "dashed") + # empirical estimate as reference
  scale_fill_manual(values = palette) + 
  ylim(c(1,20)) +
  labs(x = NULL, y = expression(paste("Mean lifetime ", italic(l))), fill = NULL) +
  theme(axis.text.x = element_blank())

g3 = ggplot(data, aes(x = model, y = deathprob, fill = model)) +
  geom_violin(color = NA) + 
  stat_pointinterval(data = data %>% filter(model != "Prior"), mapping = aes(x = model, y = deathprob, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, show.legend = F) +
  geom_hline(yintercept = d_emp, linetype = "dashed") + # empirical estimate as reference
  scale_fill_manual(values = palette) + 
  ylim(c(0,0.5)) +
  labs(x = NULL, y = expression(paste("Death probability ", italic(d))), fill = NULL) +
  theme(axis.text.x = element_blank())

g1 + ggtitle("A") + g2 + g3 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", plot.title = element_text(face= "bold"), plot.title.position = "plot")
ggsave("comparison_phylodynamic_parameters.pdf", width = 8, height = 4)


## evaluate editing and clock model parameters
# inferred editing rates per site
data_clock = rbind(adb_log %>% select(starts_with("clockRate")) %>% mutate(model = "ADB"),
                   bdsky_log %>% select(starts_with("clockRate")) %>% mutate(model = "BD")) %>%
  pivot_longer(cols = starts_with("clockRate"), names_to = "site", values_to = "clockRate") %>%
  mutate(site = factor(sub("clockRate_site_", "", site), levels = as.character(c(1:10)))) 

data_clock = rbind(data_clock, 
                   data.frame(clockRate = rlnorm(n = 1e5, meanlog = -5, sdlog = 1), model = "Prior", site = "")) %>%
  mutate(model = factor(model, levels = c("Prior", "BD", "ADB"))) %>%
  mutate(site = factor(site, levels = c("", 1:10)))

g4 = ggplot(data_clock, aes(x = site, y = clockRate, fill = model)) +
  geom_violin(color = NA, position = position_dodge(0.9)) +
  stat_pointinterval(data = data_clock %>% filter(model != "Prior"), mapping = aes(x = site, y = clockRate, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, position = position_dodge(0.9)) +
  scale_fill_manual(values = palette) + 
  ylim(c(0, 0.025)) +
  labs(x = expression(paste("Site ", italic(i))), y = expression(paste("Editing rate ", italic(r)[italic(i)])), fill = NULL) 

# inferred edit-outcome rate multipliers per site
data_scars = rbind(adb_log %>% select(starts_with("scarringRate")) %>% mutate(model = "ADB"), 
                   bdsky_log %>% select(starts_with("scarringRate")) %>% mutate(model = "BD")) %>%
  pivot_longer(cols = starts_with("scarringRate_site_"), names_to = "var", values_to = "scarringRate") %>%
  mutate(site = sub("scarringRate_site_([0-9]+)\\..*", "\\1", var),
         scar = sub(".*\\.([0-9]+)", "\\1", var)) %>%
  select(-var)

data_scars = rbind(data_scars, 
                   data.frame(scarringRate = rlnorm(n = 1e5, meanlog = 0, sdlog = 1), model = "Prior", site = "", scar = c(rep("1", 5e4), rep("2", 5e4)))) %>%
  mutate(model = factor(model, levels = c("Prior", "BD", "ADB"))) %>%
  mutate(site = factor(site, levels = c("", 1:10))) %>%
  mutate(mode = factor(scar, levels = c("1", "2")))

g5 = ggplot(data_scars, aes(x = site, y = scarringRate, fill = model)) +
  geom_violin(color = NA, position = position_dodge(0.9)) +
  stat_pointinterval(data = data_scars %>% filter(model != "Prior"), mapping = aes(x = site, y = scarringRate, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, position = position_dodge(0.9)) +
  scale_fill_manual(values = palette) + 
  facet_wrap(~scar, ncol = 2, labeller = labeller(scar = c("1" = "Inversion", "2" = "Deletion"))) +
  ylim(c(0, 2)) +
  labs(x = expression(paste("Site ", italic(i))), y = expression(paste("Scar multiplier ", italic(s)[italic(i)])), fill = NULL) 

g4 / g5 + plot_layout(guides = "collect")
ggsave("comparison_editing_parameters.pdf", width = 8, height = 7)

## ---------------------------

## assess inference of sampling probabilities under ADB
rho_log_fs = paste0("inference_adb_fixedTrees_condRoot_withRho/chain", c(1:5), ".log") 
rho_log = lapply(rho_log_fs, function(f) remove_burn_ins(parse_beast_tracelog_file(f), burn_in_fraction = 0.1)) %>% 
  bind_rows(.id = "chain") 
all(calc_esses(rho_log %>% select(-chain), sample_interval = 1000) > 200, na.rm = T)

# load ground truth
truth = read.table("tree_data.tsv", header = T)

# load get_sumstats_from_log from ADB-analysis/validation/evaluate_inference.R
parameters = colnames(rho_log)[grepl("rho", colnames(rho_log))]
df = get_sumstats_from_log(rho_log, parameters) %>%
  mutate(tree_id = sub("rho_", "", sub("_data", "", parameter))) %>%
  left_join(truth) %>%
  select(tree_id, true = rho, lower, median, upper) 

# assess coverage: consider that 1 can be approached but unlikely to be sampled
df = df %>% 
  mutate(correct = ifelse(true >= lower & true <= upper, T, F)) %>%
  mutate(correct_trunc = case_when(true < 1 & true >= lower & true <= upper ~ T,
                                   true == 1 & upper >= 0.999 ~ T,
                                   .default = F))

r1 = ggplot(df, aes(x = true, y = median, color = correct_trunc)) + 
  geom_point(size = 0.7, alpha = 0.8, show.legend = F) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.2, show.legend = F) +
  geom_abline(alpha = 0.1) +
  scale_x_continuous(limits = c(0,1), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(0,1), breaks = pretty_breaks()) +
  scale_color_manual(values = c("FALSE" = "#F8766D", "TRUE" = "#009E73")) +
  labs(title = "A", x = expression(paste("True sampling probability ", rho)), y = "Posterior median\nwith 95% HPD interval") +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot")

# plot distribution (true vs. inferred)
df_long = df %>% 
  pivot_longer(cols = c("true", "median"), names_to = "parameter", values_to = "value") %>%
  select(tree_id, parameter, value) %>%
  mutate(parameter = factor(parameter, levels = c("true", "median")))
r2 = ggplot(df_long, aes(x = parameter, y = value, fill = parameter)) +
  geom_boxplot(outlier.size = 0.7, show.legend = F) +
  scale_y_continuous(limits = c(0,1), breaks = pretty_breaks()) +
  scale_fill_manual(values = c("true" = palette[3], "median" = palette[5])) +
  labs(title = "B", x = NULL, y = expression(paste("Sampling probability ", rho))) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot")

# compare lifetime and death to previous estimates
df_cf = rbind(adb_log %>% select(shape, lifetime, deathprob) %>% mutate(model = "fixed"), 
              rho_log %>% select(shape, lifetime, deathprob) %>% mutate(model = "inferred"))

cf1 = ggplot(df_cf, aes(x = model, y = shape, fill = model)) + 
  geomViolinDiscrete::geom_violin_discrete(scale = "width", color = NA, show.legend = F) +
  stat_pointinterval(mapping = aes(x = model, y = shape, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, show.legend = F) +
  scale_fill_manual(values = palette[c(3,5)]) + 
  scale_y_continuous(limits = c(0,50), breaks = seq(10,50,10)) +
  labs(x = NULL, y = expression(paste("Shape ", italic(k))), fill = NULL) +
  theme(axis.text.x = element_blank())

cf2 = ggplot(df_cf, aes(x = model, y = lifetime, fill = model)) +
  geom_violin(color = NA) + 
  stat_pointinterval(mapping = aes(x = model, y = lifetime, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, show.legend = F) +
  geom_hline(yintercept = l_emp, linetype = "dashed") + # empirical estimate as reference
  scale_fill_manual(values = palette[c(3,5)]) + 
  ylim(c(1,20)) +
  labs(x = NULL, y = expression(paste("Mean lifetime ", italic(l))), fill = NULL) +
  theme(axis.text.x = element_blank())

cf3 = ggplot(df_cf, aes(x = model, y = deathprob, fill = model)) +
  geom_violin(color = NA) + 
  stat_pointinterval(mapping = aes(x = model, y = deathprob, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, show.legend = F) +
  geom_hline(yintercept = d_emp, linetype = "dashed") + # empirical estimate as reference
  scale_fill_manual(values = palette[c(3,5)]) + 
  ylim(c(0,0.5)) +
  labs(x = NULL, y = expression(paste("Death probability ", italic(d))), fill = NULL) +
  theme(axis.text.x = element_blank())

(r1 + r2 + plot_layout(axis_titles = "collect", widths = c(2,1))) /
(cf1 + ggtitle("C") + cf2 + cf3 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"), plot.title.position = "plot"))
ggsave("comparison_sampling.pdf", width = 8, height = 8)

## ---------------------------

## add-on: verify inference from simulations
sim_log_fs = paste0("inference_sim/chain", c(1:5), ".log") 
sim_log = lapply(sim_log_fs, function(f) remove_burn_ins(parse_beast_tracelog_file(f), burn_in_fraction = 0.1)) %>% 
  bind_rows(.id = "chain") 
all(calc_esses(sim_log %>% select(-chain), sample_interval = 1000) > 200, na.rm = T)

# get parameters
truth = read.csv("tree_sim_data.csv") %>% mutate(tree = as.character(tree))
parameters = paste0("rho_", c(1:106))
df = get_sumstats_from_log(sim_log, parameters) %>%
  mutate(tree = sub("rho_", "", parameter)) %>%
  left_join(truth) %>%
  select(tree, true = rho, lower, median, upper) 

# assess coverage: consider that 1 can be approached but unlikely to be sampled
df = df %>% 
  mutate(correct = ifelse(true >= lower & true <= upper, T, F)) %>%
  mutate(correct_trunc = case_when(true < 1 & true >= lower & true <= upper ~ T,
                                   true == 1 & upper >= 0.999 ~ T,
                                   .default = F))

ggplot(df, aes(x = true, y = median, color = correct_trunc)) + # r2 with correct_trunc
  geom_point(size = 0.7, alpha = 0.8, show.legend = F) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.2, show.legend = F) +
  geom_line(aes(x = true, y = true), inherit.aes = F, alpha = 0.1) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("FALSE" = "#F8766D", "TRUE" = "#009E73")) +
  labs(x = expression(paste("True sampling probability ", rho)), y = "Posterior median\nwith 95% HPD interval")
#r1 + r2 + plot_layout(axis_titles = "collect")

