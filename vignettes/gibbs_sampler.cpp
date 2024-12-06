#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List gibbs_sampler_rcpp(int n, double a, double b, int num_iter, int n_init = 0) {
  NumericVector x_chain(num_iter); // Chain for x
  NumericVector y_chain(num_iter); // Chain for y
  
  // Initialize
  int x = n_init; 
  double y = 0.5; 
  
  // Gibbs sampling
  for (int i = 0; i < num_iter; i++) {
    // Sample x given y (Binomial)
    x = R::rbinom(n, y);
    // Sample y given x (Beta)
    y = R::rbeta(x + a, n - x + b);
    
    // Store results
    x_chain[i] = x;
    y_chain[i] = y;
  }
  
  return List::create(Named("x") = x_chain, Named("y") = y_chain);
}
