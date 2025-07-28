#include <string>
#include <unordered_map>
#include <cmath>   
#include <algorithm> 
#include <vector>
#include "RcppArmadillo.h"
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

//-----------------------------------//
// FFT SOLVERS (NONLINEAR)
//----------------------------------//

// partial convolution using fft
// used in iterative alg
arma::vec convolve_fft(
  const arma::cx_vec& fx, 
  const arma::vec& y, 
  int n, 
  double eps) {
    
  int nz = 2 * y.n_elem;	
  arma::vec y_ext = y;
  y_ext.resize(nz);  

  arma::cx_vec fy = fft(y_ext);   
  arma::cx_vec fz = fx % fy;      
  arma::cx_vec z = ifft(fz); 
  return real(z.head(n)) * eps;   // take real part & scale 
}

arma::vec branch_initial_guess(
  const arma::cx_vec& ft,
  double d,
  int n,
  double eps) {

  arma::cx_vec rhs = ft / (1.-2.*d * ft);
  arma::cx_vec b0 = ifft(rhs);

  return (1.-d) * real(b0.head(n)) * eps;
}

// helper func
arma::vec force_nonincr(arma::vec X) {
    for (size_t i = 1; i < X.n_elem; ++i) {
        if (X(i) > X(i - 1)) {
            X(i) = X(i - 1); 
        }
    }
    return X;
}
// Extinction prob
arma::vec get_X(
  double rho, 
  double a, 
  double b, 
  double d,
  const arma::vec &t, 
  const arma::vec &pdf,
  const arma::vec &cdf,
  const arma::cx_vec &F_t,
  double dx, 
  int maxit,
  double tol) {
    
   
  arma::vec X00 = (1.0 - rho) * (1.0 - cdf) + //lives to present, unsampled
        d*cdf; //death event
  double err = 1;
  int it = 0;
  int n = X00.n_elem;
  bool cond = (1.0 - rho) > d/(1.0-d);
  arma::vec X0 = X00;
  
  // solve integral equation:
  // prob of branching but no descendants in sample
  while (err > tol && it < maxit) {
    arma::vec I = convolve_fft(F_t, X0 % X0, n, dx);
    arma::vec Xi = X00 + (1.-d)*I; 
    if (cond) {
      	Xi = force_nonincr(Xi);
	}
	err = arma::norm(X0-Xi, 2);
    X0 = Xi;
    it++;
  }

    if (it == maxit) {
      Rcpp::warning("max iterations reached with error: %f", err);
    }

  return X0;
}

// Get P_0 directly
//[[Rcpp::export]]
NumericVector get_P0(
  double a,//scale 
  int b, //shape
  double rho, //sampling
  double d,
  double t_or, // time of origin
  List controls) //list of controls for method
{
  // Get extinction probabilities (must use fftsolve)
  int	m = controls["m"];
  int max_it = controls["max_it"];
  arma::vec t_seq_fft = arma::linspace(0.0, t_or, m);
  double dx = t_or/(m-1.0);
  NumericVector t_nv = wrap(t_seq_fft);
  NumericVector cdf_nv = Rcpp::pgamma(t_nv, b, a);
  NumericVector pdf_nv = Rcpp::dgamma(t_nv, b, a);
  arma::vec pdf = as<arma::vec>(pdf_nv);
  arma::vec cdf = as<arma::vec>(cdf_nv);

  arma::vec f_t = arma::join_cols(pdf, arma::zeros(pdf.n_elem));
  arma::cx_vec F_t = fft(f_t);// FFT 
  arma::vec p0;
  p0 = get_X(rho,a,b,d,t_seq_fft,pdf,cdf,F_t,dx,max_it, 1e-12);        

  return Rcpp::NumericVector(p0.begin(), p0.end());
}

