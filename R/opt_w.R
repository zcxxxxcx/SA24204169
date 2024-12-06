#' @title A R-function for the DRO mean-variance optimization 
#' @description This file contains the R functions for the DRO mean-variance optimization. 
#' @param data Data matrix 
#' @param phi Initial weight vector 
#' @param n Number of data points to consider 
#' @param delta Robustness parameter 
#' @param alpha Risk-free rate 
#' @return List containing optimized returns and weight differences 
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(1000), ncol = 10)
#' phi <- rep(1/10, 10)
#' n <- 100
#' delta <- 0.1
#' alpha <- 0.01
#' optimize_weights(data, phi, n, delta, alpha)
#' }
#' @export
optimize_weights <- function(data, phi, n, delta, alpha) {
  opts <- list("algorithm"  = "NLOPT_LN_COBYLA",
               "xtol_rel"   = 1e-4,
               "xtol_abs"   = 1e-4,
               "maxeval"    = 1000)
  r_dro  <- c()
  dw_dro <- c()
  w0 <- phi
  
  f_dro <- function(w) ( sqrt(w%*%sigma%*%w) + sqrt(delta*sum(w^2)) )^2
  g_dro <- function(w){
    a <- alpha + sqrt(delta*sum(w^2)) - sum(w*mu)
    max(a,-w)
  }
  h_dro <- function(w) 1-sum(w)
  
  for (i in 1:(dim(data)[1] - n - 1)) {
    data_1 <- data[(1:n) + i - 1, ]
    mu <- colMeans(data_1)
    sigma <- cov(data_1)
    
    r <- nloptr::nloptr(w0,
                        eval_f      = f_dro,
                        eval_g_ineq = g_dro,
                        eval_g_eq   = h_dro,
                        opts        = opts)
    
    w_opt     <- r$solution        
    mu_real   <- data[n + i, ]     
    r_dro[i]  <- sum(w_opt * mu_real) 
    dw_dro[i] <- sum(abs(w_opt - w0)) 
    w0        <- w_opt
  }
  r_dro <- r_dro - 0.005 * dw_dro
  list(returns = r_dro, weight_diff = dw_dro)
}
