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
library(ggplot2)
library(scales)
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
write.tree(tree, 'treeBD.newick') # this tree is also saved in https://github.com/pilarskj/ADB/test_data

## compare log-likelihood curves between ADB and BDMM-Prime (calculated in BEAST2, see https://github.com/pilarskj/ADB/blob/main/test/test/adbp/BDTest.java)
methods = c('adb_approx' = 'ADB (approx)', 'adb_exact' = 'ADB (exact)', 'bdmm' = 'BDMM-Prime')
params = c("lifetime" = "'Mean lifetime'~italic(l)",
           "deathprob" = "'Death probability'~italic(d)",
           "rho" = "'Sampling probability'~rho")
true = data.frame(parameter = factor(c('lifetime', 'deathprob', 'rho'), levels = names(params)), 
                  value = c(5, 0.1, 0.1))

# load log-likelihood values (originally from ADB/test_data)
df = read.csv('loglikTreeBD.csv') %>% 
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
write.tree(tree, 'treeBD_low_sampling.newick') 

# compare log-likelihood curves 
params = c("lifetime" = "'Mean lifetime'~italic(l)",
           "deathprob" = "'Death probability'~italic(d)")
true = data.frame(parameter = factor(c('lifetime', 'deathprob'), levels = names(params)), value = c(5, 0.1))
df = read.csv('loglikTreeBD_low_sampling.csv') %>%
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
