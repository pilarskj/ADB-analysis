library(scTreeSim)
library(castor)
library(tidyverse)
library(ggplot2)
theme_set(theme_classic(base_size=14))
library(ggh4x)
## load scTreeSim 
library(scTreeSim)
## funcs for calculting dLTT plots
source("ltt_funcs.R")

set.seed(123)

# Simulate LTTs using scTreeSim
sim_LTT <- function(i, shape, rho, d, scale, t_or){
	sim <- sim_adb_origin_samp(t_or,scale,shape,d,rho,origin_type=0, min_tips=2)
	if (is.null(sim)){ # eventually remove, since we will condition on survival
		df <- data.frame(it=i, time=t_seq, lineages=NA) 
	} else {
		ltt <- count_lineages_through_time(sim@phylo, min_time =0.0, max_time = t_or, Ntimes = 100, ultrametric=TRUE)		
		# correct for root time
		corr_time <- t_or-max(ltt$time)
		df <- data.frame(it = i, time= t_or - (ltt$times+corr_time), lineages =ltt$lineages)
	}
    # add pars
    df$scale <- scale
    df$shape <- shape
    df$rho <- rho
    df$d <- d

	return(df)	
}


#---------------------------------#
# Compare simulated LTT to dLTT
#---------------------------------#
# set-up
t_or <- 10
t_seq <- seq(0, t_or, by=0.1)
nsim <- 75

# pop params:
mean_lifetime <- 1 # set time-scale
shape <- c(5,500)
rho <- c(0.08,0.3)
d <- c(0)
it <- 1:nsim
pars <- crossing(i=it, shape=shape, rho=rho, d=d)
pars$scale <- mean_lifetime/pars$shape
n <- 80
dtau <- 0.005
# calc dLTTs:
print(paste("Running dLTT for", n, "time pts with resolution", dtau)) 
df_ltt <- pars %>% pmap_dfr(~get_dLTT_df(a=..5, b=..2, rho=..3, d=..4, t_or=t_or, n=n,dtau=dtau,
    controls=list(m=2^14, max_it=100)))

print("Simulating LTT curves")
# simulate LTTs:
res <- pars %>% pmap_dfr(sim_LTT,t_or=t_or)

#-----get mean over simulated LTTs---#
# first, put on same time grid (lin approx)
res2 <-  res %>% 
	group_by(it, shape,scale, rho, d) %>%
	do({data.frame(
		time=t_seq,
		y = approx(.$time, .$lineages, xout=t_seq)$y,
		it=.$it[1])}) %>%
	ungroup()

# average over iterations
mean_obs_ltt <- res2 %>%
	group_by(time,shape,scale, rho, d) %>%
	summarise(mean=mean(y,na.rm=TRUE),
		.groups='drop')


#-----PLOT-----#
# with fixed d:
cols <- c("dLTT" = "blue","Simulated Average"= "black")
gg_final <- ggplot(res,
    aes(x=time, y=lineages, group=it))+
    geom_line(alpha=0.05)+
    geom_line(data=df_ltt, inherit.aes=FALSE, aes(x=t, y=y, colour="dLTT"), 
    linewidth=1)+
	geom_line(data=mean_obs_ltt, inherit.aes=FALSE, aes(x=time, y=mean, colour="Simulated Average"),
	linewidth=1, linetype="dotdash")+
    scale_x_reverse()+
    facet_grid(shape~rho,labeller=label_bquote(cols=Sampling~Probability~rho:.(rho),
        rows=Shape~italic(k):.(shape)), scales="free")+
    labs(x="Time since present",col="", y="Number of lineages")+
	scale_colour_manual(values=cols)+theme(legend.position="bottom")
