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
             bdsky_log %>% select(lifetime, deathprob) %>% mutate(shape = NA, model = "BD"),
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
    summarise(lower = hdi(.data[[param]], credMass = 0.95)[1], 
              median = median(.data[[param]]),
              upper = hdi(.data[[param]], credMass = 0.95)[2])
  print(summary)
}

# plot estimates
g1 = ggplot(data %>% filter(model != "BD"), aes(x = model, y = shape, fill = model)) +
  geomViolinDiscrete::geom_violin_discrete(scale = "width", color = NA, show.legend = F) +
  stat_pointinterval(data = data %>% filter(model == "ADB"), mapping = aes(x = model, y = shape, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5, show.legend = F) +
  scale_fill_manual(values = palette[c(1,3)]) + 
  scale_y_continuous(limits = c(0,50), breaks = seq(10,50,10)) +
  labs(x = NULL, y = expression(paste("Shape ", italic(k))), fill = NULL) +
  theme(axis.text.x = element_blank())

g2 = ggplot(data, aes(x = model, y = lifetime, fill = model)) +
  geom_violin(color = NA) + 
  stat_pointinterval(data = data %>% filter(model != "Prior"), mapping = aes(x = model, y = lifetime, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5) +
  geom_hline(yintercept = l_emp, linetype = "dashed") + # empirical estimate as reference
  scale_fill_manual(values = palette) + 
  ylim(c(1,20)) +
  labs(x = NULL, y = expression(paste("Mean lifetime ", italic(l))), fill = NULL) +
  theme(axis.text.x = element_blank())

g3 = ggplot(data, aes(x = model, y = deathprob, fill = model)) +
  geom_violin(color = NA) + 
  stat_pointinterval(data = data %>% filter(model != "Prior"), mapping = aes(x = model, y = deathprob, group = model), 
                     point_interval = median_hdi, .width = 0.95, point_size = 0.2, linewidth = 0.2, alpha = 0.5) +
  geom_hline(yintercept = d_emp, linetype = "dashed") + # empirical estimate as reference
  scale_fill_manual(values = palette) + 
  ylim(c(0,0.5)) +
  labs(x = NULL, y = expression(paste("Death probability ", italic(d))), fill = NULL) +
  theme(axis.text.x = element_blank())


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
  facet_wrap(~scar, ncol = 1, labeller = labeller(scar = c("1" = "Inversion", "2" = "Deletion"))) +
  ylim(c(0, 2)) +
  labs(x = expression(paste("Site ", italic(i))), y = expression(paste("Scar multiplier ", italic(s)[italic(i)])), fill = NULL) 


## joint plot
plot1 = g1 + ggtitle("A") + g2 + g3 + plot_layout(guides = "collect", widths = c(0.7,1,1)) & 
  theme(legend.position = "none", plot.title = element_text(face= "bold"), plot.title.position = "plot")
plot2 = g4 + g5 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
plot1 / plot2
ggsave("comparison_inferred_parameters.pdf", width = 10, height = 9)

