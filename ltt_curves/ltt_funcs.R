##========================================##
# Funcs to calc expected lineages thru time
##========================================##
library(Rcpp)
sourceCpp("extinction.cpp") # to get P0 calculation


get_dLTT_df <- function(
    a, #scale
    b, #shape
    rho,#sampling
    d,#death
    t_or,#time of origin
    n,#number of time pts to evaluate dLTT
    dtau, #resolution for solver
	controls#control to be passed to P0 solver
){
	# set-up time-grid for dLTT
	dt <- t_or/(n-1)
	t <- seq(0,t_or, by=dt)
	# calc dLTT curve
	y <- get_dLTT(a,b,rho,d,t_or,n,dtau,controls)
	return(data.frame(t=t,y=y,scale=a, shape=b, d=d, rho=rho))
}

get_dLTT <- function(
    a, #scale
    b, #shape
    rho,#sampling
    d,#death
    t_or,#time of origin
    n,#number of time pts to evaluate dLTT
	dtau, #resolution for solver
    controls,#control to be passed to P0 solver
    method="BE",#method for solving dLTT--default: backward euler
    final_size=F){#if onyl want expected final size
	#-------------------#
	# set-up   
    dt <- t_or / (n - 1)
    m <- controls$m
    max_it <- controls$max_it
    P_0 <- get_P0(a, b, rho, d, t_or, controls)
    t_0 <- seq(0, t_or, by = t_or / (m - 1))
  	# calc 
    if (final_size) {
      LTT <- solve_integral_LTT_backward_euler(0, P_0, t_0, t_or, a, b, d, dtau,complete)
    } else {
      t_seq <- seq(0, t_or - dt, by = dt)
      LTT <- sapply(t_seq, solve_integral_LTT_backward_euler, 
			P_0, t_0, t_or, a, b, d, dtau, FALSE)
      LTT <- c(LTT, NA)
    }
  	# condition on survival
    cond <- 1 - P_0[length(P_0)]
    return(LTT / cond)
}


