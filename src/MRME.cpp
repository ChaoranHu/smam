// [[Rcpp::interfaces(r, cpp)]]

#include "movres.h"
#include "thmam_likelihood.h"


// [[Rcpp::export]]
double test(double x) {
  double result = dcoga2dim(x, 1, 1, 3, 2);
  return result; 
}
