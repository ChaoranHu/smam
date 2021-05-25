// [[Rcpp::interfaces(r, cpp)]]

#include "movres.h"

/***************************************************************
 Moving-Moving process implementation
***************************************************************/



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

// integrand for the h..star functions; has to be of type integr_fn.
// pointer *ex carryies a double vector of
// (t, sigma1, sigma0, lambda1, lambda0, dim, x)

void f11mm(double *w, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma1   = ptr[1];
  double sigma0  = ptr[2];
  double lambda1 = ptr[3];
  double lambda0 = ptr[4];
  int    dim     = (int) ptr[5];
  double *x      = ptr + 6;
  for (int i = 0; i < n; i++) {
    double temp = p11(w[i], t, lambda1, lambda0);
    double sd = sqrt(pow(sigma1, 2) * w[i] + pow(sigma0, 2) * (t - w[i]));
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    w[i] = temp;
  }
}

// only different from f11 by replacing  p11 with p10
// there must be a better way to do this: put a function pointer to *ex
void f10mm(double *w, int n, void *ex){
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma1   = ptr[1];
  double sigma0  = ptr[2];
  double lambda1 = ptr[3];
  double lambda0 = ptr[4];
  int    dim     = (int) ptr[5];
  double *x      = ptr + 6;
  for (int i = 0; i < n; i++){
    double temp = p10(w[i], t, lambda1, lambda0);
    double sd = sqrt(pow(sigma1, 2) * w[i] + pow(sigma0, 2) * (t - w[i]));
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    w[i] = temp;
  }
}

void f00mm(double *w, int n, void *ex){
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma1   = ptr[1];
  double sigma0  = ptr[2];
  double lambda1 = ptr[3];
  double lambda0 = ptr[4];
  int    dim     = (int) ptr[5];
  double *x      = ptr + 6;
  for (int i = 0; i < n; i++){
    double temp = p00(w[i], t, lambda1, lambda0);
    double sd = sqrt(pow(sigma0, 2) * w[i] + pow(sigma1, 2) * (t - w[i]));
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    w[i] = temp;
  }
}

// only different from f00 by replacing  p00 with p01
// there must be a better way to do this: put a function pointer to *ex
void f01mm(double *w, int n, void *ex){
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma1   = ptr[1];
  double sigma0   = ptr[2];
  double lambda1 = ptr[3];
  double lambda0 = ptr[4];
  int    dim     = (int) ptr[5];
  double *x      = ptr + 6;
  for (int i = 0; i < n; i++){
    double temp = p01(w[i], t, lambda1, lambda0);
    double sd = sqrt(pow(sigma0, 2) * w[i] + pow(sigma1, 2) * (t - w[i]));
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    w[i] = temp;
  }
}

// h functions needs numerical integration
// [[Rcpp::export]]
NumericVector h11mm(NumericMatrix x, NumericVector t, NumericVector theta,
		    NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma1 = theta[2], sigma0 = theta[3];
  /* set up for Rdqags */
  double *ex = R_Calloc(6 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = R_Calloc(limit, int);
  double *work = R_Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma1; ex[2] = sigma0; ex[3] = lambda1; ex[4] = lambda0;
  ex[5] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    double sd = sqrt(pow(sigma1, 2) * t[i]);
    double prod = exp(-lambda1 * t[i]);
    for (int j = 0; j < dim; j++) {
      ex[6 + j] = x(i, j);
      prod *= R::dnorm(x(i, j), 0.0, sd, 0);
    }
    b = t[i]; ex[0] = t[i];
    Rdqags(f11mm, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result + prod;
  }
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector h10mm(NumericMatrix x, NumericVector t, NumericVector theta,
		    NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma1 = theta[2], sigma0 = theta[3];
  /* set up for Rdqags */
  double *ex = R_Calloc(6 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork =   R_Calloc(limit, int);
  double *work = R_Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma1; ex[2] = sigma0; ex[3] = lambda1; ex[4] = lambda0;
  ex[5] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[6 + j] = x(i, j);
    b = t[i]; ex[0] = t[i];
    Rdqags(f10mm, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector h00mm(NumericMatrix x, NumericVector t, NumericVector theta,
		    NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma1 = theta[2], sigma0 = theta[3];
  /* set up for Rdqags */
  double *ex = R_Calloc(6 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork =   R_Calloc(limit, int);
  double *work = R_Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma1; ex[2] = sigma0; ex[3] = lambda1; ex[4] = lambda0;
  ex[5] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    double sd = sqrt(pow(sigma0, 2) * t[i]);
    double prod = exp(-lambda0 * t[i]);
    for (int j = 0; j < dim; j++) {
      ex[6 + j] = x(i, j);
      prod *= R::dnorm(x(i, j), 0.0, sd, 0);
    }
    b = t[i]; ex[0] = t[i];
    Rdqags(f00mm, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result + prod;
  }
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector h01mm(NumericMatrix x, NumericVector t, NumericVector theta,
		    NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma1 = theta[2], sigma0 = theta[3];
  /* set up for Rdqags */
  double *ex = R_Calloc(6 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork =   R_Calloc(limit, int);
  double *work = R_Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma1; ex[2] = sigma0; ex[3] = lambda1; ex[4] = lambda0;
  ex[5] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[6 + j] = x(i, j);
    b = t[i]; ex[0] = t[i];
    Rdqags(f01mm, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}


/****************************************************************************
 likelihood function using the forward algorithm for hidden Markov model
****************************************************************************/
// [[Rcpp::export]]
double nllk_inc_mm(NumericVector &theta, NumericMatrix &data,
		   NumericVector &integrControl, LogicalVector &logtr) {
  if (logtr[0]) theta = exp(theta);
  if (is_true( any(theta <= 0.) )) return(NA_REAL);
  if (theta[2] <= theta[3]) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  double lambda1 = theta[0], lambda0 = theta[1];
  double pm = 1. / lambda1 / (1. / lambda1 + 1. / lambda0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix x  = data(Range(0, n - 1), Range(1, dim));
  NumericVector
    hmm = h11mm(x, tt, theta, integrControl),
    hmr = h10mm(x, tt, theta, integrControl),
    hrr = h00mm(x, tt, theta, integrControl),
    hrm = h01mm(x, tt, theta, integrControl);
  double alpha0 = pr, alpha1 = pm;
  
  // forward algorithm in this loop
  double llk = 0.;
  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    
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

