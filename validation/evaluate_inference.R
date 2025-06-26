## ---------------------------
##
## Script for evaluating phylodynamic inference under the Age-Dependent Branching model (validation)
##
## Author: Julia Pilarski
## Email: julia.pilarski@bsse.ethz.ch
##
## Date Created: 2025-06-02
##
## ---------------------------

setwd("~/Projects/P1_AgeDependentTrees/ADB-analysis/validation")

library(tracerer)
library(HDInterval)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
theme_set(theme_classic(base_size = 14))

## ---------------------------

# Function for calculating the median and 95% HPD interval for all inferred parameters from a log table
get_sumstats_from_log <- function(log, parameters){
  stats = sapply(parameters, function(x) {
    # get posterior estimates
    ix = match(x, colnames(log))
    posterior = log[, ix]
    # calculate median and HPD interval  
    hpd = hdi(posterior, credMass = 0.95)
    stats = c('lower' = hpd[[1]], 'median' = median(posterior), 'upper' = hpd[[2]])
  }) %>% 
    t() %>% as.data.frame() %>% rownames_to_column('parameter')
  return(stats)
}

# Function for processing a log file
get_estimates <- function(log_fs, parameters) {
  
  # get parameter estimates
  data = lapply(log_fs, function(f) {
    
    # parse beast log 
    log = remove_burn_ins(parse_beast_tracelog_file(f), burn_in_fraction = 0.1)
    
    # get summary statistics
    stats = get_sumstats_from_log(log, parameters)
    
    # add ID and chain length
    stats$ID = str_extract(sapply(str_split(f, '/'), tail, 1), "[0-9]+")
    stats$chainLength = log$Sample[nrow(log)]
    
    # get ESS per parameter
    ess = calc_esses(log, sample_interval = 1000)
    stats$ESS = as.numeric(ess[match(stats$parameter, names(ess))])
    
    return(stats)
  }) %>% bind_rows()
  
  return(data)
}

## ---------------------------

## load ground truth parameters
truth = read.csv("tree_data.csv")

## load inferred parameters
log_fs = paste0('inference/inference_', c(1:100), '.log')
parameters = c('shape', 'lifetime', 'deathprob', 'rho')
data = get_estimates(log_fs, parameters)

## check convergence
data %>% filter(ESS < 200) %>% pull(ID) %>% unique() # 3 chains did not converge: 16, 38, 42
data %>% filter(is.na(ESS)) # all have shape = 1
truth[data %>% filter(is.na(ESS)) %>% pull(ID) %>% as.numeric(), ] # indeed, the real shape is 1

## align estimates & calculate coverage
df = data %>% 
  left_join(truth %>%
              pivot_longer(cols = parameters, names_to = "parameter", values_to = "true") %>%
              select(ID = tree, parameter, true) %>%
              mutate(ID = as.character(ID)), by = c("ID", "parameter")) %>%
  mutate(correct = ifelse(true >= lower & true <= upper, T, F),
         parameter = factor(parameter, levels = parameters)) %>%
  filter(!ID %in% c("16", "38", "42")) # remove erroneous chains from downstream analysis
cov = df %>% group_by(parameter) %>% summarize(cov = round(sum(correct), 1)) 

## plot results
ggplot(df, aes(x = true, y = median, color = correct)) +
  geom_point(size = 0.7, show.legend = F) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.5, show.legend = F) +
  geom_text(data = cov, mapping = aes(x = -Inf, y = -Inf, label = paste(cov, "%")), hjust = -6.5, vjust = -1, inherit.aes = F) + 
  geom_line(aes(x = true, y = true), inherit.aes = F, alpha = 0.25) + 
  scale_color_manual(values = c("FALSE" = "#F8766D", "TRUE" = "#009E73")) + 
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(x = "True value", y = "Posterior median\nwith 95% HPD interval") +
  facet_wrap(~parameter, scales = "free", labeller = as_labeller(
    c("shape" = "'Shape'~italic(k)", "lifetime" = "'Mean lifetime'~italic(l)", 
      "deathprob" = "'Death probability'~italic(d)", "rho" = "'Sampling probability'~italic(rho)"), 
    label_parsed))
ggsave("validation_95HPD.pdf", width = 8, height = 6)







