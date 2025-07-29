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


#Backwards Euler
#A function numerically approximating the LTT curve using backward euler method
#This version is significantly faster

solve_integral_LTT_backward_euler <- function(t, P_0, t_0, t_or, a, b, d, dtau, complete) {
  tau_seq <- seq(t, t_or, by = dtau)
  N <- length(tau_seq)
  
  f_tau <- dgamma(tau_seq - t, shape = b, scale = a)
  if (complete) {
    P_vals <- rep(0, N)
  } else {
    P_vals <- approx(t_0, P_0, xout = tau_seq, method = "linear", rule = 2)$y
  }
  
  # Initialization
  M <- rep(0, N)
  M00 <- (1 - pgamma(tau_seq - t, shape = b, scale = a)) * (1 - P_vals)
  M_old <- M00
  max_it <- 100
  tol <- 1e-8
  
  M[1] <- 1/(1-2*(1-d)*dtau*dgamma(0, shape = b, scale = a)) * M00[1]
  M[2] <- 1/(1-2*(1-d)*dtau*dgamma(0, shape = b, scale = a)) * (M00[2] + 2*(1-d)*dtau*dgamma(dtau, shape = b, scale = a)*M[1])
  for (j in 3:N){
    M[j] <-  1/(1-2*(1-d)*dtau*dgamma(0, shape = b, scale = a)) * (M00[j] + 2*(1-d)*dtau*sum((rev(f_tau[2:j]))*M[1:(j-1)]))
  }
  
  return(M[length(M)])
}


#Trapezoid
#A function numerically approximating the LTT curve using trapezoid method
#This version is not used in practice as it is significantly slower

convolution_trapezoid <- function(tau, dt, ft, gt, s, tau_seq){
  idx_tau <- which(abs(tau_seq - tau) < dt/2)
  idx_range <- 1:idx_tau
  
   f_vals <- ft[idx_range]
  g_vals <- gt[idx_range]
  f_vals <- rev(f_vals)
  
  # Apply trapezoidal rule
  integral <- sum((f_vals[-1] + f_vals[-length(f_vals)]) / 2 *
                    (g_vals[-1] + g_vals[-length(g_vals)]) / 2) * dt
  
  return(integral)
}

solve_integral_LTT_trapezoid <- function(t,P_0,t_0, t_or,a,b,d, dt, complete){
  tau_seq <- seq(t, t_or, by=dt)
  if(t==15){
    cat(tau_seq)
  }
  f_tau <- dgamma(tau_seq-t, shape=b, scale=a)
  if (complete){
    P_0 <- 0
  } else {
    P_0 <- approx(t_0, P_0, xout=tau_seq, method="linear")$y 
  }

  M00 <- (1-pgamma(tau_seq - t, shape=b, scale=a))*(1-P_0)
  M0 <- M00#
  it <- 0
  err <- 1
  max_it <- 100

  while (err > 1e-8 & it < max_it){
    g_t <- M0
    I <- sapply(tau_seq, function(tau) {
      convolution_trapezoid(tau, dt, f_tau, g_t, t, tau_seq)
    })
    M <- M00 + 2*(1 - d)*I
    err <- norm(M - M0, type = "2")
    M0 <- M
    it <- it + 1
    if (it == max_it) message("warning: max iter reached")
  }
  
  return(M0[length(M0)])
}

