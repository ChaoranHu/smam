// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppGSL, RcppParallel)]]

#include "thmam_likelihood.h"

/******************************************************************************
 ****************** forward variables and backward variables ******************
 **************************** normalized variables ****************************
 ******************************* serial version *******************************
 ******************************************************************************/

// [[Rcpp::export]]
NumericMatrix fwd_bwd_ths(NumericVector &theta, NumericMatrix &data,
	                  NumericVector &integrControl) {
  // theta lambda0, lambda1, lambda2, sigma, p
  // data diff of t and x
  int n = data.nrow(); int dim = data.ncol() - 1;
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double p = theta[4];
  // if (lambda1 < lambda2) return NA_REAL;
  double ps0 = 1. / lambda0 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps1 = p / lambda1 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps2 = (1 - p) / lambda2 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  NumericVector tt = data.column(0);
  NumericMatrix x = data(Range(0, n - 1), Range(1, dim));

  // result matrix: frist three col for forward
  //                last  three col for backward
  NumericMatrix result(n + 1, 6);
  result(0, 0) = ps0; result(0, 1) = ps1; result(0, 2) = ps2;
  result(n, 3) =   1; result(n, 4) =   1; result(n, 5) =   1;
  NumericVector dx(n);

  // calculate all h functions
  NumericVector
    hresult00 = ths_h00(x, tt, theta, integrControl),
    hresult01 = ths_h01(x, tt, theta, integrControl),
    hresult02 = ths_h02(x, tt, theta, integrControl),
    hresult10 = ths_h10(x, tt, theta, integrControl),
    hresult11 = ths_h11(x, tt, theta, integrControl),
    hresult12 = ths_h12(x, tt, theta, integrControl),
    hresult20 = ths_h20(x, tt, theta, integrControl),
    hresult21 = ths_h21(x, tt, theta, integrControl),
    hresult22 = ths_h22(x, tt, theta, integrControl);
  
  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true(all(crow == 0.))) {
      hresult00[i] = 0.;
      hresult01[i] = 0.;
      hresult02[i] = 0.;
      hresult10[i] = 0.;
      hresult11[i] = exp(-lambda1 * tt[i]);
      hresult12[i] = 0.;
      hresult20[i] = 0.;
      hresult21[i] = 0.;
      hresult22[i] = exp(-lambda2 * tt[i]);
    }
  }

  // forward algorithm
  for (int i = 0; i < n; i++) {
    double sumf0 = result(i, 0) * hresult00[i] + result(i, 1) * hresult10[i] + result(i, 2) * hresult20[i];
    double sumf1 = result(i, 0) * hresult01[i] + result(i, 1) * hresult11[i] + result(i, 2) * hresult21[i];
    double sumf2 = result(i, 0) * hresult02[i] + result(i, 1) * hresult12[i] + result(i, 2) * hresult22[i];
    dx[i] = sumf0 + sumf1 + sumf2;
    result(i + 1, 0) = sumf0 / dx[i];
    result(i + 1, 1) = sumf1 / dx[i];
    result(i + 1, 2) = sumf2 / dx[i];
  }

  //backward algorithm
  for (int i = 0; i < n; i++) {
    double sumb0 = result(n-i, 3) * hresult00[n-i-1] + result(n-i, 4) * hresult01[n-i-1] + result(n-i, 5) * hresult02[n-i-1];
    double sumb1 = result(n-i, 3) * hresult10[n-i-1] + result(n-i, 4) * hresult11[n-i-1] + result(n-i, 5) * hresult12[n-i-1];
    double sumb2 = result(n-i, 3) * hresult20[n-i-1] + result(n-i, 4) * hresult21[n-i-1] + result(n-i, 5) * hresult22[n-i-1];
    result(n-i-1, 3) = sumb0 / dx[n-i-1];
    result(n-i-1, 4) = sumb1 / dx[n-i-1];
    result(n-i-1, 5) = sumb2 / dx[n-i-1];
  }

  return result;
}





/******************************************************************************
 ******************************* Viterbi Algorithm ****************************
 ******************************* serial version *******************************
 ****************** note that the result is log-likelihood of path ************
 ******************************************************************************/

// [[Rcpp::export]]
NumericMatrix viterbi_ths(NumericVector &theta, NumericMatrix &data,
	                  NumericVector &integrControl) {
  // theta lambda0, lambda1, lambda2, sigma, p
  // data diff of t and x
  int n = data.nrow(); int dim = data.ncol() - 1;
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double p = theta[4];
  // if (lambda1 < lambda2) return NA_REAL;
  double ps0 = 1. / lambda0 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps1 = p / lambda1 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps2 = (1 - p) / lambda2 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  NumericVector tt = data.column(0);
  NumericMatrix x = data(Range(0, n - 1), Range(1, dim));

  // result matrix: three cols stand for Viterbi prob of state 0,1,2
  //                at current time points. For numerical reason,
  //                the log-prob is returned.
  NumericMatrix result(n + 1, 3);
  result(0, 0) = log(ps0); result(0, 1) = log(ps1); result(0, 2) = log(ps2);
  NumericVector cartV = result.row(0);
  NumericVector cartW(3);

  // calculate all h functions
  NumericVector
    hresult00 = ths_h00(x, tt, theta, integrControl),
    hresult01 = ths_h01(x, tt, theta, integrControl),
    hresult02 = ths_h02(x, tt, theta, integrControl),
    hresult10 = ths_h10(x, tt, theta, integrControl),
    hresult11 = ths_h11(x, tt, theta, integrControl),
    hresult12 = ths_h12(x, tt, theta, integrControl),
    hresult20 = ths_h20(x, tt, theta, integrControl),
    hresult21 = ths_h21(x, tt, theta, integrControl),
    hresult22 = ths_h22(x, tt, theta, integrControl);
  
  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true(all(crow == 0.))) {
      hresult00[i] = 0.;
      hresult01[i] = 0.;
      hresult02[i] = 0.;
      hresult10[i] = 0.;
      hresult11[i] = exp(-lambda1 * tt[i]);
      hresult12[i] = 0.;
      hresult20[i] = 0.;
      hresult21[i] = 0.;
      hresult22[i] = exp(-lambda2 * tt[i]);
    }
  }

  // calculate Viterbi path
  for (int i = 1; i <= n; i++) {
    cartW[0] = cartV[0] + log(hresult00[i - 1]);
    cartW[1] = cartV[1] + log(hresult10[i - 1]);
    cartW[2] = cartV[2] + log(hresult20[i - 1]);
    result(i, 0) = max(cartW);

    cartW[0] = cartV[0] + log(hresult01[i - 1]);
    cartW[1] = cartV[1] + log(hresult11[i - 1]);
    cartW[2] = cartV[2] + log(hresult21[i - 1]);
    result(i, 1) = max(cartW);

    cartW[0] = cartV[0] + log(hresult02[i - 1]);
    cartW[1] = cartV[1] + log(hresult12[i - 1]);
    cartW[2] = cartV[2] + log(hresult22[i - 1]);
    result(i, 2) = max(cartW);

    cartV = result.row(i);
  }

  return result;
}


/***** Partial Viterbi Path *****/

// [[Rcpp::export]]
NumericMatrix partial_viterbi_ths(NumericVector &theta, NumericMatrix &data,
				  NumericVector &integrControl,
				  int &startpoint, int &pathlength){
  // theta lambda0, lambda1, lambda2, sigma, p
  // data diff of t and x
  // startpoint the start time point, note that
  //            the first time point in data is t0
  // pathlength the length of partial viterbi path
  int n = data.nrow(); int dim = data.ncol() - 1;
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double p = theta[4];
  double ps0 = 1. / lambda0 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps1 = p / lambda1 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps2 = (1 - p) / lambda2 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  NumericVector tt = data.column(0);
  NumericMatrix x = data(Range(0, n - 1), Range(1, dim));

  // bf_result matrix: frist three col for forward
  //                   last  three col for backward
  NumericMatrix bf_result(n + 1, 6);
  bf_result(0, 0) = ps0; bf_result(0, 1) = ps1; bf_result(0, 2) = ps2;
  bf_result(n, 3) =   1; bf_result(n, 4) =   1; bf_result(n, 5) =   1;
  NumericVector dx(n);
  
  // result matrix: three cols stand for Viterbi prob of state 0,1,2
  //                at current time points. For numerical reason,
  //                the log-prob is returned.
  NumericMatrix result(pathlength, 3);
  NumericVector cartV(3);
  NumericVector cartW(3);
  
  

  // calculate all h functions
  NumericVector
    hresult00 = ths_h00(x, tt, theta, integrControl),
    hresult01 = ths_h01(x, tt, theta, integrControl),
    hresult02 = ths_h02(x, tt, theta, integrControl),
    hresult10 = ths_h10(x, tt, theta, integrControl),
    hresult11 = ths_h11(x, tt, theta, integrControl),
    hresult12 = ths_h12(x, tt, theta, integrControl),
    hresult20 = ths_h20(x, tt, theta, integrControl),
    hresult21 = ths_h21(x, tt, theta, integrControl),
    hresult22 = ths_h22(x, tt, theta, integrControl);
  
  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true(all(crow == 0.))) {
      hresult00[i] = 0.;
      hresult01[i] = 0.;
      hresult02[i] = 0.;
      hresult10[i] = 0.;
      hresult11[i] = exp(-lambda1 * tt[i]);
      hresult12[i] = 0.;
      hresult20[i] = 0.;
      hresult21[i] = 0.;
      hresult22[i] = exp(-lambda2 * tt[i]);
    }
  }

  // forward algorithm
  for (int i = 0; i < n; i++) {
    double sumf0 = bf_result(i, 0) * hresult00[i] + bf_result(i, 1) * hresult10[i] + bf_result(i, 2) * hresult20[i];
    double sumf1 = bf_result(i, 0) * hresult01[i] + bf_result(i, 1) * hresult11[i] + bf_result(i, 2) * hresult21[i];
    double sumf2 = bf_result(i, 0) * hresult02[i] + bf_result(i, 1) * hresult12[i] + bf_result(i, 2) * hresult22[i];
    dx[i] = sumf0 + sumf1 + sumf2;
    bf_result(i + 1, 0) = sumf0 / dx[i];
    bf_result(i + 1, 1) = sumf1 / dx[i];
    bf_result(i + 1, 2) = sumf2 / dx[i];
  }

  // backward algorithm
  for (int i = 0; i < n; i++) {
    double sumb0 = bf_result(n-i, 3) * hresult00[n-i-1] + bf_result(n-i, 4) * hresult01[n-i-1] + bf_result(n-i, 5) * hresult02[n-i-1];
    double sumb1 = bf_result(n-i, 3) * hresult10[n-i-1] + bf_result(n-i, 4) * hresult11[n-i-1] + bf_result(n-i, 5) * hresult12[n-i-1];
    double sumb2 = bf_result(n-i, 3) * hresult20[n-i-1] + bf_result(n-i, 4) * hresult21[n-i-1] + bf_result(n-i, 5) * hresult22[n-i-1];
    bf_result(n-i-1, 3) = sumb0 / dx[n-i-1];
    bf_result(n-i-1, 4) = sumb1 / dx[n-i-1];
    bf_result(n-i-1, 5) = sumb2 / dx[n-i-1];
  }

  // prepare for viterbi path
  result(0, 0) = log(bf_result(startpoint, 0));
  result(0, 1) = log(bf_result(startpoint, 1));
  result(0, 2) = log(bf_result(startpoint, 2));
  cartV = result.row(0);
  int ite_stop = startpoint + pathlength - 2;

  // viterbi algorithm
  for (int i = startpoint; i < ite_stop; i++) {
    cartW[0] = cartV[0] + log(hresult00[i]);
    cartW[1] = cartV[1] + log(hresult10[i]);
    cartW[2] = cartV[2] + log(hresult20[i]);
    result(i - startpoint + 1, 0) = max(cartW);

    cartW[0] = cartV[0] + log(hresult01[i]);
    cartW[1] = cartV[1] + log(hresult11[i]);
    cartW[2] = cartV[2] + log(hresult21[i]);
    result(i - startpoint + 1, 1) = max(cartW);

    cartW[0] = cartV[0] + log(hresult02[i]);
    cartW[1] = cartV[1] + log(hresult12[i]);
    cartW[2] = cartV[2] + log(hresult22[i]);
    result(i - startpoint + 1, 2) = max(cartW);

    cartV = result.row(i - startpoint + 1);
  }

  // last step of viterbi algorithm
  cartW[0] = cartV[0] + log(hresult00[ite_stop]) + log(bf_result(ite_stop + 1, 3));
  cartW[1] = cartV[1] + log(hresult10[ite_stop]) + log(bf_result(ite_stop + 1, 3));
  cartW[2] = cartV[2] + log(hresult20[ite_stop]) + log(bf_result(ite_stop + 1, 3));
  result(pathlength - 1, 0) = max(cartW);

  cartW[0] = cartV[0] + log(hresult01[ite_stop]) + log(bf_result(ite_stop + 1, 4));
  cartW[1] = cartV[1] + log(hresult11[ite_stop]) + log(bf_result(ite_stop + 1, 4));
  cartW[2] = cartV[2] + log(hresult21[ite_stop]) + log(bf_result(ite_stop + 1, 4));
  result(pathlength - 1, 1) = max(cartW);

  cartW[0] = cartV[0] + log(hresult02[ite_stop]) + log(bf_result(ite_stop + 1, 5));
  cartW[1] = cartV[1] + log(hresult12[ite_stop]) + log(bf_result(ite_stop + 1, 5));
  cartW[2] = cartV[2] + log(hresult22[ite_stop]) + log(bf_result(ite_stop + 1, 5));
  result(pathlength - 1, 2) = max(cartW);

  return result;
}
