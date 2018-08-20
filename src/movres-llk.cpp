// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

/******************************************************************
 Densities for time spent in 1-cycle or 0-cycle in interval (0, t]
******************************************************************/

// To derive distribution of increments of the BMT process X(t) we need the
// joint distributions of the two pairs  M(t),S(t)  and  R(t),S(t)  for
// a given initial state S(0), where M(t) and R(t), t > 0, are the total time
// in interval (0, t] spent in the 1-cycles and in the 0-cycles, respectively;
// consequently, R(t) = t − M (t).
//
// Following the notations in \citet{Yan:etal:2014}, define
// \begin{equation*}
//  P_1\big[\cdot\big] = \Pr \big[\cdot\,|S(0) = 1\big],
//  \qquad\mbox{and}\qquad
//  P_0\big[\cdot\big] = \Pr \big[\cdot\,|S(0) = 0\big].
// \end{equation*}
// Then, for $0 < w < t$, we introduce the following (defective) densities
// \begin{align*}
//   p_{11}(w,t)\dif w &=P_1\big[M(t)\in \dif w, S(t)=1\big],\\
//   p_{10}(w,t)\dif w &=P_1\big[M(t)\in \dif w, S(t)=0\big],\\
//   p_{01}(w,t)\dif w &=P_0\big[R(t)\in \dif w, S(t)=1\big],\\
//   p_{00}(w,t)\dif w &=P_0\big[R(t)\in \dif w, S(t)=0\big].
// \end{align*}


double p11(double w, double t, double lambda1, double lambda0) {
    if (w > t || w < 0) return(0.);
    double lm = lambda1 * w, lr = lambda0 * (t - w);
    double z = 2 * sqrt(lm * lr); 
    return(exp(- lm - lr) * sqrt( lambda1 * lambda0 * w / (t - w)) *
	   R::bessel_i(z, 1, 1));
}

// [[Rcpp::export]]
NumericVector vp11(NumericVector vw, double t, double lambda1, double lambda0) {
    int n = vw.size();
    NumericVector dens(n);
    for (int i = 0; i < n; i++) {
	dens[i] = p11(vw[i], t, lambda1, lambda0);
    }
    return(dens);
}


double p10(double w, double t, double lambda1, double lambda0) {
  if (w > t || w < 0) return(0.);
  double lm = lambda1 * w, lr = lambda0 * (t - w);
  double z = 2 * sqrt(lm * lr); 
  return(lambda1 * exp(-lm - lr) * R::bessel_i(z, 0, 1));
}


// [[Rcpp::export]]
NumericVector vp10(NumericVector vw, double t, double lambda1, double lambda0) {
    int n = vw.size();
    NumericVector dens(n);
    for (int i = 0; i < n; i++) {
	dens[i] = p10(vw[i], t, lambda1, lambda0);
    }
    return(dens);
}


// for convenience; just switch lambda0 and lambda1 in p11
double p00(double w, double t, double lambda1, double lambda0) {
  double dens = p11(w, t, lambda0, lambda1);
  return(dens);
}

// [[Rcpp::export]]
NumericVector vp00(NumericVector vw, double t, double lambda1, double lambda0) {
  return(vp11(vw, t, lambda0, lambda1));
}

// For convenience; just switch lambda0 and lambda1 in p10
double p01(double w, double t, double lambda1, double lambda0) {
  double dens = p10(w, t, lambda0, lambda1);
  return(dens);
}

// [[Rcpp::export]]
NumericVector vp01(NumericVector vw, double t, double lambda1, double lambda0) {
  return(vp10(vw, t, lambda0, lambda1));
}


/***************************************************************
 Interface of integration
***************************************************************/
// integr_fn is the type of function Rdqags expects
// typedef void integr_fn(double *x, int n, void *ex);
// Access point for numerical integration from the C code of R
// void Rdqags(integr_fn f, void *ex, double *a, double *b,
// 	    double *epsabs, double *epsrel,
// 	    double *result, double *abserr, int *neval, int *ier,
// 	    int *limit, int *lenw, int *last,
// 	    int *iwork, double *work);

// integrand for the h.. functions; has to be of type integr_fn.
// pointer *ex carryies a double vector of
// (t, signa, lambda0, lambda1, dim, x)

void f11(double *w, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda1 = ptr[2];
  double lambda0 = ptr[3];
  int    dim     = (int) ptr[4];
  double *x      = ptr + 5;
  for (int i = 0; i < n; i++) {
    double temp = p11(w[i], t, lambda1, lambda0);
    double sd = sigma * sqrt(w[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    w[i] = temp;
  }
}

// only different from f11 by replacing  p11 with p10
// there must be a better way to do this: put a function pointer to *ex
void f10(double *w, int n, void *ex){
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda1 = ptr[2];
  double lambda0 = ptr[3];
  int    dim     = (int) ptr[4];
  double *x      = ptr + 5;
  for (int i = 0; i < n; i++){
    double temp = p10(w[i], t, lambda1, lambda0);
    double sd = sigma * sqrt(w[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    w[i] = temp;
  }
}

void f00(double *w, int n, void *ex){
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda1 = ptr[2];
  double lambda0 = ptr[3];
  int    dim     = (int) ptr[4];
  double *x      = ptr + 5;
  for (int i = 0; i < n; i++){
    double temp = p00(w[i], t, lambda1, lambda0);
    double sd = sigma * sqrt(t - w[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    w[i] = temp;
  }
}

// only different from f00 by replacing  p00 with p01
// there must be a better way to do this: put a function pointer to *ex
void f01(double *w, int n, void *ex){
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda1 = ptr[2];
  double lambda0 = ptr[3];
  int    dim     = (int) ptr[4];
  double *x      = ptr + 5;
  for (int i = 0; i < n; i++){
    double temp = p01(w[i], t, lambda1, lambda0);
    double sd = sigma * sqrt(t - w[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    w[i] = temp;
  }
}

// h functions needs numerical integration
// [[Rcpp::export]]
NumericVector h11(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2];
  /* set up for Rdqags */
  double *ex = Calloc(5 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda1; ex[3] = lambda0;
  ex[4] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    double sd = sigma * sqrt(t[i]);
    double prod = exp(-lambda1 * t[i]);
    for (int j = 0; j < dim; j++) {
      ex[5 + j] = x(i, j);
      prod *= R::dnorm(x(i, j), 0.0, sd, 0);
    }
    b = t[i]; ex[0] = t[i];
    Rdqags(f11, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result + prod;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector h10(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2];
  /* set up for Rdqags */
  double *ex = Calloc(5 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork =   Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda1; ex[3] = lambda0;
  ex[4] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[5 + j] = x(i, j);
    b = t[i]; ex[0] = t[i];
    Rdqags(f10, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector h00(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2];
  /* set up for Rdqags */
  double *ex = Calloc(5 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork =   Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda1; ex[3] = lambda0;
  ex[4] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[5 + j] = x(i, j);
    b = t[i]; ex[0] = t[i];
    Rdqags(f00, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector h01(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2];
  /* set up for Rdqags */
  double *ex = Calloc(5 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork =   Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda1; ex[3] = lambda0;
  ex[4] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[5 + j] = x(i, j);
    b = t[i]; ex[0] = t[i];
    Rdqags(f01, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}

/******************************************************************************
 Composite likelihood 
******************************************************************************/
// [[Rcpp::export]]
double ncllk_m1_inc(NumericVector &theta, NumericMatrix &data,
		    NumericVector &integrControl, LogicalVector &logtr) {
  if (logtr[0]) theta = exp(theta);
  if (is_true( any(theta <= 0.) )) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  double lambda1 = theta[0], lambda0 = theta[1];
  double pm = 1. / lambda1 / (1. / lambda1 + 1. / lambda0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix x  = data(Range(0, n - 1), Range(1, dim));
  // can be made more efficient by only doing this for the moving transitions
  NumericVector m1dens =
    pm * (h11(x, tt, theta, integrControl) + h10(x, tt, theta, integrControl)) +
    pr * (h01(x, tt, theta, integrControl) + h00(x, tt, theta, integrControl));
  // process resting points
  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true( all(crow == 0.) )) m1dens[i] = pr * exp(-lambda0 * tt[i]);
  }
  return( - sum(log(m1dens)));
}


/****************************************************************************
 likelihood function using the forward algorithm for hidden Markov model
****************************************************************************/
// [[Rcpp::export]]
double nllk_inc(NumericVector &theta, NumericMatrix &data,
		NumericVector &integrControl, LogicalVector &logtr) {
  if (logtr[0]) theta = exp(theta);
  if (is_true( any(theta <= 0.) )) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  double lambda1 = theta[0], lambda0 = theta[1];
  double pm = 1. / lambda1 / (1. / lambda1 + 1. / lambda0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix x  = data(Range(0, n - 1), Range(1, dim));
  NumericVector
    hmm = h11(x, tt, theta, integrControl),
    hmr = h10(x, tt, theta, integrControl),
    hrr = h00(x, tt, theta, integrControl),
    hrm = h01(x, tt, theta, integrControl);
  double alpha0 = pr, alpha1 = pm;
  
  // forward algorithm in this loop
  double llk = 0.;
  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true( all(crow == 0.) )) { // resting
      hmm[i] = 0.; hrr[i] = 0; hrm[i] = 0.;
      hrr[i] = exp( -lambda0 * tt[i]); 
    }
    // ending with resting
    double sumfr = alpha0 * hrr[i] + alpha1 * hmr[i];
    // ending with moving
    double sumfm = alpha0 * hrm[i] + alpha1 * hmm[i];
 
    double dx = sumfr + sumfm;
    alpha0 = sumfr / dx;
    alpha1 = sumfm / dx;
    llk += log(dx);
  }
  return( - llk);
}




/*******************************************
  forward variable and backward
*******************************************/
// [[Rcpp::export]]
NumericMatrix fwd_bwd_mr(NumericVector &theta, NumericMatrix &data,
	                 NumericVector &integrControl) {
  // theta lambdaM (lambda1), lambdaR(lambda0), sigma
  if (is_true( any(theta <= 0.) )) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  double lambda1 = theta[0], lambda0 = theta[1];
  double pm = 1. / lambda1 / (1. / lambda1 + 1. / lambda0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix x  = data(Range(0, n - 1), Range(1, dim));


  NumericVector
    hmm = h11(x, tt, theta, integrControl),
    hmr = h10(x, tt, theta, integrControl),
    hrr = h00(x, tt, theta, integrControl),
    hrm = h01(x, tt, theta, integrControl);

  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true(all(crow == 0.))) {
      hmm[i] = 0.; hrr[i] = 0; hrm[i] = 0.;
      hrr[i] = exp( -lambda0 * tt[i]);
    }
  }

  // result matrix: frist two col for forward
  //                last  two col for backward
  //                moving first, resting second
  NumericMatrix result(n + 1, 4);
  result(0, 0) = pm; result(0, 1) = pr;
  result(n, 2) =  1; result(n, 3) =  1;
  NumericVector dx(n);



    
  // forward algorithm
  for (int i = 0; i < n; i++) {
    // ending with resting
    double sumfr = result(i, 1) * hrr[i] + result(i, 0) * hmr[i];
    // ending with moving
    double sumfm = result(i, 1) * hrm[i] + result(i, 0) * hmm[i];

    dx[i] = sumfr + sumfm;
    result(i + 1, 1) = sumfr / dx[i];
    result(i + 1, 0) = sumfm / dx[i];
  }

  //backward algorithm
  for (int i = 0; i < n; i++) {
    // starting from resting
    double sumbr = result(n-i, 2) * hrm[n-i-1] + result(n-i, 3) * hrr[n-i-1];
    // starting from moving
    double sumbm = result(n-i, 2) * hmm[n-i-1] + result(n-i, 3) * hmr[n-i-1];

    result(n-i-1, 2) = sumbm / dx[n-i-1];
    result(n-i-1, 3) = sumbr / dx[n-i-1];
  }

  return result;
}


/*******************************
  full viterbi path
*******************************/

// [[Rcpp::export]]
NumericMatrix viterbi_mr(NumericVector &theta, NumericMatrix &data,
	                 NumericVector &integrControl) {
  // theta lambdaM (lambda1), lambdaR(lambda0), sigma
  if (is_true( any(theta <= 0.) )) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  double lambda1 = theta[0], lambda0 = theta[1];
  double pm = 1. / lambda1 / (1. / lambda1 + 1. / lambda0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix x  = data(Range(0, n - 1), Range(1, dim));


  NumericVector
    hmm = h11(x, tt, theta, integrControl),
    hmr = h10(x, tt, theta, integrControl),
    hrr = h00(x, tt, theta, integrControl),
    hrm = h01(x, tt, theta, integrControl);

  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true(all(crow == 0.))) {
      hmm[i] = 0.; hrr[i] = 0; hrm[i] = 0.;
      hrr[i] = exp( -lambda0 * tt[i]);
    }
  }

  // result matrix: three cols stand for Viterbi prob of state moving, resting
  //                at current time points. For numerical reason, the log-prob is returned.
  NumericMatrix result(n + 1, 2);
  result(0, 0) = log(pm); result(0, 1) = log(pr);
  NumericVector cartV = result.row(0);
  NumericVector cartW(2);

  // calculate Viterbi path
  for (int i = 1; i <= n; i++) {
    cartW[0] = cartV[0] + log(hmm[i - 1]);
    cartW[1] = cartV[1] + log(hrm[i - 1]);
    result(i, 0) = max(cartW);

    cartW[0] = cartV[0] + log(hmr[i - 1]);
    cartW[1] = cartV[1] + log(hrr[i - 1]);
    result(i, 1) = max(cartW);

    cartV = result.row(i);
  }

  return result;
}


/********************************
  partial viterbi path
*********************************/

// [[Rcpp::export]]
NumericMatrix partial_viterbi_mr(NumericVector &theta, NumericMatrix &data,
				 NumericVector &integrControl,
				 int &startpoint, int &pathlength){
  // data diff of t and x
  // startpoint the start time point, note that
  //            the first time point in data is t0
  // pathlength the length of partial viterbi path
  // theta lambdaM (lambda1), lambdaR(lambda0), sigma
  if (is_true( any(theta <= 0.) )) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  double lambda1 = theta[0], lambda0 = theta[1];
  double pm = 1. / lambda1 / (1. / lambda1 + 1. / lambda0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix x  = data(Range(0, n - 1), Range(1, dim));

  // bf_result matrix: frist two col for forward
  //                   last  two col for backward
  //                   moving first, resting second
  NumericMatrix bf_result(n + 1, 4);
  bf_result(0, 0) = pm; bf_result(0, 1) = pr;
  bf_result(n, 2) =  1; bf_result(n, 3) =  1;
  NumericVector dx(n);
  
  // result matrix: two cols stand for Viterbi prob of state moving, resting
  //                at current time points. For numerical reason, the log-prob is returned.
  NumericMatrix result(pathlength, 2);
  NumericVector cartV(2);
  NumericVector cartW(2);

  // calculate all h functions
  NumericVector
    hmm = h11(x, tt, theta, integrControl),
    hmr = h10(x, tt, theta, integrControl),
    hrr = h00(x, tt, theta, integrControl),
    hrm = h01(x, tt, theta, integrControl);

  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true(all(crow == 0.))) {
      hmm[i] = 0.; hrr[i] = 0; hrm[i] = 0.;
      hrr[i] = exp( -lambda0 * tt[i]);
    }
  }

  // forward algorithm
  for (int i = 0; i < n; i++) {
    // ending with resting
    double sumfr = bf_result(i, 1) * hrr[i] + bf_result(i, 0) * hmr[i];
    // ending with moving
    double sumfm = bf_result(i, 1) * hrm[i] + bf_result(i, 0) * hmm[i];

    dx[i] = sumfr + sumfm;
    bf_result(i + 1, 1) = sumfr / dx[i];
    bf_result(i + 1, 0) = sumfm / dx[i];
  }
  

  // backward algorithm
  for (int i = 0; i < n; i++) {
    // starting from resting
    double sumbr = bf_result(n-i, 2) * hrm[n-i-1] + bf_result(n-i, 3) * hrr[n-i-1];
    // starting from moving
    double sumbm = bf_result(n-i, 2) * hmm[n-i-1] + bf_result(n-i, 3) * hmr[n-i-1];

    bf_result(n-i-1, 2) = sumbm / dx[n-i-1];
    bf_result(n-i-1, 3) = sumbr / dx[n-i-1];
  }

  // prepare for viterbi path
  result(0, 0) = log(bf_result(startpoint, 0));
  result(0, 1) = log(bf_result(startpoint, 1));
  cartV = result.row(0);
  int ite_stop = startpoint + pathlength - 2;

  // viterbi algorithm
  for (int i = startpoint; i < ite_stop; i++) {
    cartW[0] = cartV[0] + log(hmm[i]);
    cartW[1] = cartV[1] + log(hrm[i]);
    result(i - startpoint + 1, 0) = max(cartW);

    cartW[0] = cartV[0] + log(hmr[i]);
    cartW[1] = cartV[1] + log(hrr[i]);
    result(i - startpoint + 1, 1) = max(cartW);

    cartV = result.row(i - startpoint + 1);
  }

  // last step of viterbi algorithm
  cartW[0] = cartV[0] + log(hmm[ite_stop]) + log(bf_result(ite_stop + 1, 2));
  cartW[1] = cartV[1] + log(hrm[ite_stop]) + log(bf_result(ite_stop + 1, 2));
  result(pathlength - 1, 0) = max(cartW);

  cartW[0] = cartV[0] + log(hmr[ite_stop]) + log(bf_result(ite_stop + 1, 3));
  cartW[1] = cartV[1] + log(hrr[ite_stop]) + log(bf_result(ite_stop + 1, 3));
  result(pathlength - 1, 1) = max(cartW);

  return result;
}
