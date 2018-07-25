#ifndef THMAM_LIKELIHOOD_H
#define THMAM_LIKELIHOOD_H

#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
#include <R_ext/Applic.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;


NumericVector ths_h00(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h01(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h02(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h10(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h11(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h12(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h20(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h21(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h22(NumericMatrix x, NumericVector t, NumericVector theta,
		      NumericVector integrControl);

NumericVector ths_h00_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

NumericVector ths_h01_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

NumericVector ths_h02_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

NumericVector ths_h10_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

NumericVector ths_h11_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

NumericVector ths_h12_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

NumericVector ths_h20_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

NumericVector ths_h21_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

NumericVector ths_h22_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize);

double nllk_fwd_ths(NumericVector &theta, NumericMatrix &data,
		    NumericVector &integrControl);

double nllk_fwd_ths_parallel(NumericVector &theta, NumericMatrix &data,
	                     NumericVector &integrControl, int grainSize);



#endif
