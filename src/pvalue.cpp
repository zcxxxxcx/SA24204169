#include <Rcpp.h>
using namespace Rcpp;

// Calculate the covariance between two vectors
double covariance(NumericVector x, NumericVector y) {
  int n = x.size();
  double mean_x = mean(x);
  double mean_y = mean(y);
  double cov = 0.0;
  
  for (int i = 0; i < n; ++i) {
    cov += (x[i] - mean_x) * (y[i] - mean_y);
  }
  
  return cov / (n - 1);
}

//' @name pvalue
//' @title A function to calculate the p-value for the Sharpe ratio
//' @description This function calculates the p-value for the Sharpe ratio
//' @param r Numeric vector of returns
//' @param r0 Numeric vector of benchmark returns
//' @return p-value for the Sharpe ratio
//' @examples
//' \dontrun{
//' r <- rnorm(100)
//' r0 <- rnorm(100)
//' pvalue(r, r0)
//' }
//' @export
// [[Rcpp::export]]
double pvalue(NumericVector r, NumericVector r0) {
  int n = r.size();  // Number of observations
  double sin = covariance(r, r0);  // Calculate the covariance between r and r0
  double si = sd(r);  // Standard deviation of r
  double sn = sd(r0);  // Standard deviation of r0
  double mi = mean(r);  // Mean of r
  double mn = mean(r0);  // Mean of r0
  double Shin = sn * mi - si * mn;  // Sharpe ratio numerator
  double v = 1.0 / n * (2 * pow(si * sn, 2) - 2 * si * sn * sin + 
                        0.5 * pow(mi * sn, 2) + 0.5 * pow(mn * si, 2) - 
                        mi * mn / (2 * si * sn) * (pow(sin, 2) + pow(si * sn, 2)));  // Variance
  double p = 2 * R::pnorm(-abs(Shin) / sqrt(v), 0.0, 1.0, 1, 0);  // p-value
  return p;
}