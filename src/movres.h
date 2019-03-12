#ifndef MOVRES_H
#define MOVRES_H

#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

NumericVector h11(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl);

NumericVector h10(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl);

NumericVector h00(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl);

NumericVector h01(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl);

#endif
