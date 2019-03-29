// [[Rcpp::interfaces(r, cpp)]]

#include "movres.h"
#include "thmam_likelihood.h"



// reduce h function to one dim situation
double h11_mrme(double x, double t, NumericVector theta,
		NumericVector integrControl) {
  NumericVector vec_x = {x};
  vec_x.attr("dim") = Dimension(1, 1);
  NumericMatrix mat_x = as<NumericMatrix>(vec_x);
  NumericVector vec_t = {t};
  double result = h11(mat_x, vec_t, theta, integrControl)[0];

  return result;
}

double h10_mrme(double x, double t, NumericVector theta,
		NumericVector integrControl) {
  NumericVector vec_x = {x};
  vec_x.attr("dim") = Dimension(1, 1);
  NumericMatrix mat_x = as<NumericMatrix>(vec_x);
  NumericVector vec_t = {t};
  double result = h10(mat_x, vec_t, theta, integrControl)[0];

  return result;
}

double h01_mrme(double x, double t, NumericVector theta,
		NumericVector integrControl) {
  NumericVector vec_x = {x};
  vec_x.attr("dim") = Dimension(1, 1);
  NumericMatrix mat_x = as<NumericMatrix>(vec_x);
  NumericVector vec_t = {t};
  double result = h01(mat_x, vec_t, theta, integrControl)[0];

  return result;
}

double h00_mrme(double x, double t, NumericVector theta,
		NumericVector integrControl) {
  NumericVector vec_x = {x};
  vec_x.attr("dim") = Dimension(1, 1);
  NumericMatrix mat_x = as<NumericMatrix>(vec_x);
  NumericVector vec_t = {t};
  double result = h00(mat_x, vec_t, theta, integrControl)[0];

  return result;
}

// evaluation of g function
// g_{ij}(z, t) = P_i(Z(t) \in dz, S(t) = j)/dz
// we hold one dim first.
// then, calculate high dim with product.

// pointer *ex carryies a double array of
// (t, sigma, lambda1, lambda0, sig_err, integrControl(3), z(1))

void g11_integrand_mrme(double *w, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  NumericVector theta = {ptr[2], ptr[3], ptr[1]}; //lam1, lam0, sigma
  double sig_err = ptr[4];
  NumericVector integrControl = {ptr[5], ptr[6], ptr[7]};
  double z       = ptr[8];

  for (int i = 0; i < n; i++) {
    double temp = h11_mrme(z-w[i], t, theta, integrControl);
    temp *= R::dnorm(w[i], 0.0, 2*pow(sig_err, 2), 0);
    w[i] = temp;
  }
}

void g10_integrand_mrme(double *w, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  NumericVector theta = {ptr[2], ptr[3], ptr[1]}; //lam1, lam0, sigma
  double sig_err = ptr[4];
  NumericVector integrControl = {ptr[5], ptr[6], ptr[7]};
  double z       = ptr[8];

  for (int i = 0; i < n; i++) {
    double temp = h10_mrme(z-w[i], t, theta, integrControl);
    temp *= R::dnorm(w[i], 0.0, 2*pow(sig_err, 2), 0);
    w[i] = temp;
  }
}

void g01_integrand_mrme(double *w, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  NumericVector theta = {ptr[2], ptr[3], ptr[1]}; //lam1, lam0, sigma
  double sig_err = ptr[4];
  NumericVector integrControl = {ptr[5], ptr[6], ptr[7]};
  double z       = ptr[8];

  for (int i = 0; i < n; i++) {
    double temp = h01_mrme(z-w[i], t, theta, integrControl);
    temp *= R::dnorm(w[i], 0.0, 2*pow(sig_err, 2), 0);
    w[i] = temp;
  }
}

void g00_integrand_mrme(double *w, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  NumericVector theta = {ptr[2], ptr[3], ptr[1]}; //lam1, lam0, sigma
  double sig_err = ptr[4];
  NumericVector integrControl = {ptr[5], ptr[6], ptr[7]};
  double z       = ptr[8];

  for (int i = 0; i < n; i++) {
    double temp = h00_mrme(z-w[i], t, theta, integrControl);
    temp *= R::dnorm(w[i], 0.0, 2*pow(sig_err, 2), 0);
    w[i] = temp;
  }
}

// g function: theta(lam1, lam0, sigma, sig_err)
// [[Rcpp::export]]
NumericVector g11_mrme(NumericMatrix z, NumericVector t,
		       NumericVector theta, NumericVector integrControl) {
  int dim = z.ncol(), n = z.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2], sig_err = theta[3];
  /* set up for Rdqagi */
  double *ex = Calloc(9, double);
  double bound = 0; // not use
  int inf   = 2; // integration interval (-inf, inf)
  // input
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
  ex[4] = sig_err; ex[5] = integrControl[0]; ex[6] = integrControl[1];
  ex[7] = integrControl[2];
  NumericVector value(n);
  double cart = 1;
  for (int i = 0; i < n; i++) {
    ex[0] = t[i];
    for (int j = 0; j < dim; j++) {
      ex[8] = z(i, j);
      Rdqagi(g11_integrand_mrme, ex, &bound, &inf, &epsabs, &epsrel,
	     &result, &abserr, &neval, &ier, &limit, &lenw, &last,
	     iwork, work);
      cart *= result;
    }
    value[i] = cart;
    cart = 1;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector g10_mrme(NumericMatrix z, NumericVector t,
		       NumericVector theta, NumericVector integrControl) {
  int dim = z.ncol(), n = z.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2], sig_err = theta[3];
  /* set up for Rdqagi */
  double *ex = Calloc(9, double);
  double bound = 0; // not use
  int inf   = 2; // integration interval (-inf, inf)
  // input
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
  ex[4] = sig_err; ex[5] = integrControl[0]; ex[6] = integrControl[1];
  ex[7] = integrControl[2];
  NumericVector value(n);
  double cart = 1;
  for (int i = 0; i < n; i++) {
    ex[0] = t[i];
    for (int j = 0; j < dim; j++) {
      ex[8] = z(i, j);
      Rdqagi(g10_integrand_mrme, ex, &bound, &inf, &epsabs, &epsrel,
	     &result, &abserr, &neval, &ier, &limit, &lenw, &last,
	     iwork, work);
      cart *= result;
    }
    value[i] = cart;
    cart = 1;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector g01_mrme(NumericMatrix z, NumericVector t,
		       NumericVector theta, NumericVector integrControl) {
  int dim = z.ncol(), n = z.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2], sig_err = theta[3];
  /* set up for Rdqagi */
  double *ex = Calloc(9, double);
  double bound = 0; // not use
  int inf   = 2; // integration interval (-inf, inf)
  // input
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
  ex[4] = sig_err; ex[5] = integrControl[0]; ex[6] = integrControl[1];
  ex[7] = integrControl[2];
  NumericVector value(n);
  double cart = 1;
  for (int i = 0; i < n; i++) {
    ex[0] = t[i];
    for (int j = 0; j < dim; j++) {
      ex[8] = z(i, j);
      Rdqagi(g01_integrand_mrme, ex, &bound, &inf, &epsabs, &epsrel,
	     &result, &abserr, &neval, &ier, &limit, &lenw, &last,
	     iwork, work);
      cart *= result;
    }
    value[i] = cart;
    cart = 1;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector g00_mrme(NumericMatrix z, NumericVector t,
		       NumericVector theta, NumericVector integrControl) {
  int dim = z.ncol(), n = z.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2], sig_err = theta[3];
  /* set up for Rdqagi */
  double *ex = Calloc(9, double);
  double bound = 0; // not use
  int inf   = 2; // integration interval (-inf, inf)
  // input
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
  ex[4] = sig_err; ex[5] = integrControl[0]; ex[6] = integrControl[1];
  ex[7] = integrControl[2];
  NumericVector value(n);
  double cart = 1;
  for (int i = 0; i < n; i++) {
    ex[0] = t[i];
    for (int j = 0; j < dim; j++) {
      ex[8] = z(i, j);
      Rdqagi(g00_integrand_mrme, ex, &bound, &inf, &epsabs, &epsrel,
	     &result, &abserr, &neval, &ier, &limit, &lenw, &last,
	     iwork, work);
      result += exp(-lambda0*ex[0])*R::dnorm(ex[8], 0.0, 2*pow(sig_err, 2), 0);
      cart *= result;
    }
    value[i] = cart;
    cart = 1;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// evaluation of t function
// t_{ij}(t) = P_i(S(t) = j)
// param: t: time diff; theta: (lam1, lam0, sigma, sig_err).
// [[Rcpp::export]]
NumericVector t11_mrme(NumericVector t, NumericVector theta) {
  double lam1 = theta[0], lam0 = theta[1];
  int n = t.length();
  NumericVector result(n);
  double cartA, cartB, cartlast;
  int j;

  for (int i = 0; i < n; i++) {
    j = 1;
    cartA = 1 - R::pgamma(t[i], 1, 1/lam1, 1, 0);
    cartB = 0; cartlast = 0;
    while(TRUE) {
      cartB = pcoga2dim_diff_shape(t[i], j, j, lam1, lam0);
      if (cartB == R_PosInf || R_IsNaN(cartB)) {
	warning("Inf or NaN happened, not converge!");
	break;
      }
      cartA += cartB;
      if (cartB == 0 && cartlast >= cartB && j > 1) break;
      cartlast = cartB;
      j++;
    }
    result[i] = cartA;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector t00_mrme(NumericVector t, NumericVector theta) {
  double lam1 = theta[0], lam0 = theta[1];
  int n = t.length();
  NumericVector result(n);
  double cartA, cartB, cartlast;
  int j;

  for (int i = 0; i < n; i++) {
    j = 1;
    cartA = 1 - R::pgamma(t[i], 1, 1/lam0, 1, 0);
    cartB = 0; cartlast = 0;
    while(TRUE) {
      cartB = pcoga2dim_diff_shape(t[i], j, j, lam0, lam1);
      if (cartB == R_PosInf || R_IsNaN(cartB)) {
	warning("Inf or NaN happened, not converge!");
	break;
      }
      cartA += cartB;
      if (cartB == 0 && cartlast >= cartB && j > 1) break;
      cartlast = cartB;
      j++;
    }
    result[i] = cartA;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector t10_mrme(NumericVector t, NumericVector theta) {
  double lam1 = theta[0], lam0 = theta[1];
  int n = t.length();
  NumericVector result(n);
  double cartA, cartB, cartlast;
  int j;

  for (int i = 0; i < n; i++) {
    j = 0;
    cartA = 0;
    cartB = 0; cartlast = 0;
    while(TRUE) {
      cartB = pcoga2dim_diff_shape(t[i], j, j+1, lam0, lam1);
      if (cartB == R_PosInf || R_IsNaN(cartB)) {
	warning("Inf or NaN happened, not converge!");
	break;
      }
      cartA += cartB;
      if (cartB == 0 && cartlast >= cartB && j > 1) break;
      cartlast = cartB;
      j++;
    }
    result[i] = cartA;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector t01_mrme(NumericVector t, NumericVector theta) {
  double lam1 = theta[0], lam0 = theta[1];
  int n = t.length();
  NumericVector result(n);
  double cartA, cartB, cartlast;
  int j;

  for (int i = 0; i < n; i++) {
    j = 0;
    cartA = 0;
    cartB = 0; cartlast = 0;
    while(TRUE) {
      cartB = pcoga2dim_diff_shape(t[i], j, j+1, lam1, lam0);
      if (cartB == R_PosInf || R_IsNaN(cartB)) {
	warning("Inf or NaN happened, not converge!");
	break;
      }
      cartA += cartB;
      if (cartB == 0 && cartlast >= cartB && j > 1) break;
      cartlast = cartB;
      j++;
    }
    result[i] = cartA;
  }
  return result;
}


// negative log-likelihood of MVME
// theta: c(lam1, lam0, sigma, sig_err)
// data: diff time locations

/*
// [[Rcpp::export]]
double nllk_mrme(NumericVector &theta, NumericMatrix &data,
		 NumericVector &integrControl)
*/
