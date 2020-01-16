// [[Rcpp::interfaces(r, cpp)]]

#include "movres.h"
#include "thmam_likelihood.h"


// evaluate int_R dnorm(z-x, 0, b)*dnorm(x, 0, d) dx
// ex carryies (z, b, d)
// b and d here are s.d.
void norm_integrand_mrme(double *w, int n, void *ex) {
  double *ptr = (double *) ex;
  double z = ptr[0];
  double b = ptr[1];
  double d = ptr[2];

  for (int i = 0; i < n; i++) {
    double temp = R::dnorm(z-w[i], 0 , b, 0) * R::dnorm(w[i], 0, d, 0);
    w[i] = temp;
  }
}

// [[Rcpp::export]]
double norm_mrme(double z, double b, double d,
		 NumericVector integrControl) {
  double *ex = Calloc(3, double);
  double bound = 0; // not use
  int inf = 2; // integration interval (-inf, inf)
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

  ex[0] = z; ex[1] = b; ex[2] = d;
  Rdqagi(norm_integrand_mrme, ex, &bound, &inf, &epsabs, &epsrel,
	 &result, &abserr, &neval, &ier, &limit, &lenw, &last,
	 iwork, work);
  Free(ex); Free(iwork); Free(work);
  return(result);
}


// evaluation of g function
// g_{ij}(z, t) = P_i(Z(t) \in dz, S(t) = j)/dz
// we hold one dim first.
// then, calculate high dim with product.

// pointer *ex carryies a double array of
// (t, sigma, lambda1, lambda0, sig_err, integrControl(3), dim, z(1), z(2), ...)

// evaluate g10
void g10_integrand_mrme(double *w, int n, void *ex) {
  double *ptr       = (double *) ex;
  double t          = ptr[0];
  double sigma      = ptr[1];
  double lambda1    = ptr[2];
  double lambda0    = ptr[3];
  double sig_err    = ptr[4];
  NumericVector integrControl = {ptr[5], ptr[6], ptr[7]};
  int dim           = (int) ptr[8];
  double *z         = ptr + 9;
  for (int i = 0; i < n; i++) {
    double temp = p10(w[i], t, lambda1, lambda0);
    double sd1 = sqrt(w[i] * pow(sigma, 2));
    double sd2 = sqrt(2    * pow(sig_err, 2));
    for (int j = 0; j < dim; j++){
      temp *= norm_mrme(z[j], sd1, sd2, integrControl);
    }
    w[i] = temp;
  }
}

// [[Rcpp::export]]
NumericVector g10_mrme(NumericMatrix z, NumericVector t,
		       NumericVector theta, NumericVector integrControl) {
  int dim = z.ncol(), n = z.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2], sig_err = theta[3];
  /* set up for Rdqags */
  double *ex = Calloc(9 + dim, double);
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
  ex[4] = sig_err; ex[5] = integrControl[0]; ex[6] = integrControl[1];
  ex[7] = integrControl[2]; ex[8] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[9 + j] = z(i, j);
    b = t[i]; ex[0] = t[i];
    Rdqags(g10_integrand_mrme, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// evaluate g01
void g01_integrand_mrme(double *w, int n, void *ex) {
  double *ptr       = (double *) ex;
  double t          = ptr[0];
  double sigma      = ptr[1];
  double lambda1    = ptr[2];
  double lambda0    = ptr[3];
  double sig_err    = ptr[4];
  NumericVector integrControl = {ptr[5], ptr[6], ptr[7]};
  int dim           = (int) ptr[8];
  double *z         = ptr + 9;
  for (int i = 0; i < n; i++) {
    double temp = p01(w[i], t, lambda1, lambda0);
    double sd1 = sqrt((t - w[i]) * pow(sigma, 2));
    double sd2 = sqrt(2    * pow(sig_err, 2));
    for (int j = 0; j < dim; j++){
      temp *= norm_mrme(z[j], sd1, sd2, integrControl);
    }
    w[i] = temp;
  }
}

// [[Rcpp::export]]
NumericVector g01_mrme(NumericMatrix z, NumericVector t,
		       NumericVector theta, NumericVector integrControl) {
  int dim = z.ncol(), n = z.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2], sig_err = theta[3];
  /* set up for Rdqags */
  double *ex = Calloc(9 + dim, double);
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
  ex[4] = sig_err; ex[5] = integrControl[0]; ex[6] = integrControl[1];
  ex[7] = integrControl[2]; ex[8] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[9 + j] = z(i, j);
    b = t[i]; ex[0] = t[i];
    Rdqags(g01_integrand_mrme, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// evaluate g00
void g00_integrand_mrme(double *w, int n, void *ex) {
  double *ptr       = (double *) ex;
  double t          = ptr[0];
  double sigma      = ptr[1];
  double lambda1    = ptr[2];
  double lambda0    = ptr[3];
  double sig_err    = ptr[4];
  NumericVector integrControl = {ptr[5], ptr[6], ptr[7]};
  int dim           = (int) ptr[8];
  double *z         = ptr + 9;
  for (int i = 0; i < n; i++) {
    double temp = p00(w[i], t, lambda1, lambda0);
    double sd1 = sqrt((t - w[i]) * pow(sigma, 2));
    double sd2 = sqrt(2    * pow(sig_err, 2));
    for (int j = 0; j < dim; j++){
      temp *= norm_mrme(z[j], sd1, sd2, integrControl);
    }
    w[i] = temp;
  }
}

// [[Rcpp::export]]
NumericVector g00_mrme(NumericMatrix z, NumericVector t,
		       NumericVector theta, NumericVector integrControl) {
  int dim = z.ncol(), n = z.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2], sig_err = theta[3];
  /* set up for Rdqags */
  double *ex = Calloc(9 + dim, double);
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
  ex[4] = sig_err; ex[5] = integrControl[0]; ex[6] = integrControl[1];
  ex[7] = integrControl[2]; ex[8] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    double atom = exp(-lambda0*t[i]);
    for (int j = 0; j < dim; j++) {
      ex[9 + j] = z(i, j);
      atom *= R::dnorm(z(i, j), 0, sqrt(2 * pow(sig_err, 2)), 0);
    }
    b = t[i]; ex[0] = t[i];
    Rdqags(g00_integrand_mrme, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    
    value[i] = result + atom;
  }
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// evaluate g11
void g11_integrand_mrme(double *w, int n, void *ex) {
  double *ptr       = (double *) ex;
  double t          = ptr[0];
  double sigma      = ptr[1];
  double lambda1    = ptr[2];
  double lambda0    = ptr[3];
  double sig_err    = ptr[4];
  NumericVector integrControl = {ptr[5], ptr[6], ptr[7]};
  int dim           = (int) ptr[8];
  double *z         = ptr + 9;
  for (int i = 0; i < n; i++) {
    double temp = p11(w[i], t, lambda1, lambda0);
    double sd1 = sqrt(w[i] * pow(sigma, 2));
    double sd2 = sqrt(2    * pow(sig_err, 2));
    for (int j = 0; j < dim; j++){
      temp *= norm_mrme(z[j], sd1, sd2, integrControl);
    }
    w[i] = temp;
  }
}

// [[Rcpp::export]]
NumericVector g11_mrme(NumericMatrix z, NumericVector t,
		       NumericVector theta, NumericVector integrControl) {
  int dim = z.ncol(), n = z.nrow();
  double lambda1 = theta[0], lambda0 = theta[1], sigma = theta[2], sig_err = theta[3];
  /* set up for Rdqags */
  double *ex = Calloc(9 + dim, double);
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
  ex[4] = sig_err; ex[5] = integrControl[0]; ex[6] = integrControl[1];
  ex[7] = integrControl[2]; ex[8] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    double atom = exp(-lambda1*t[i]);
    for (int j = 0; j < dim; j++) {
      ex[9 + j] = z(i, j);
      atom *= norm_mrme(z(i, j), sqrt(t[i] * pow(sigma, 2)), sqrt(2 * pow(sig_err, 2)), integrControl);
    }
    b = t[i]; ex[0] = t[i];
    Rdqags(g11_integrand_mrme, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    
    value[i] = result + atom;
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


// negative log-likelihood of MRME
// theta: c(lam1, lam0, sigma, sig_err)
// data: diff time locations
// nllk_mrme = nllk_chain1 + nllk_chain2
// chain1 starts from Z_0
// chain2 starts from Z_1

// [[Rcpp::export]]
double nllk_mrme(NumericVector &theta, NumericMatrix &data,
		 NumericVector &integrControl) {
  if (is_true(any(theta <= 0))) return(NA_REAL);
  if (theta[2] <= theta[3]) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  if (n < 2) {
    warning("Sample size is too small to process, should be at least 3. Return nllk as 0.");
    return(0);
  }
  double lam1 = theta[0], lam0 = theta[1];
  double pm = 1. / lam1 / (1. / lam1 + 1. / lam0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix z  = data(Range(0, n - 1), Range(1, dim));
  NumericVector
    gmm = g11_mrme(z, tt, theta, integrControl),
    grr = g00_mrme(z, tt, theta, integrControl),
    grm = g01_mrme(z, tt, theta, integrControl),
    gmr = g10_mrme(z, tt, theta, integrControl);
  NumericVector
    tmm = t11_mrme(tt, theta),
    trr = t00_mrme(tt, theta),
    trm = t01_mrme(tt, theta),
    tmr = t10_mrme(tt, theta);

  // forward algorithm for the first chain
  // start from Z_0
  double alpha0 = pr, alpha1 = pm;//alpha's are normalized forward variable
  double llk1 = 0;
  double cartrr = 0, cartrm = 0, cartmm = 0, cartmr = 0;
  double dx = 0, sumfr = 0, sumfm = 0;
  for (int i = 0; i < std::floor(static_cast <double> (n) / 2); i++) {
    cartrr = trm[2*i]*gmr[2*i+1] + trr[2*i]*grr[2*i+1];
    cartmr = tmm[2*i]*gmr[2*i+1] + tmr[2*i]*grr[2*i+1];
    cartrm = trm[2*i]*gmm[2*i+1] + trr[2*i]*grm[2*i+1];
    cartmm = tmr[2*i]*grm[2*i+1] + tmm[2*i]*gmm[2*i+1];
    sumfr = cartrr * alpha0 + cartmr * alpha1;
    sumfm = cartrm * alpha0 + cartmm * alpha1;
    dx = sumfr + sumfm;
    alpha0 = sumfr / dx;
    alpha1 = sumfm / dx;
    llk1 += log(dx);
  }

  // forward algorithm for the second chain
  // start from Z_1
  alpha0 = pr, alpha1 = pm;//alpha's are normalized forward variable
  double llk2 = 0;
  cartrr = 0, cartrm = 0, cartmm = 0, cartmr = 0;
  dx = 0, sumfr = 0, sumfm = 0;
  // we have to reset all carts first.
  for(int i = 0; i < std::floor(static_cast <double> (n-1) / 2); i++) {
    cartrr = trm[2*i+1]*gmr[2*i+2] + trr[2*i+1]*grr[2*i+2];
    cartmr = tmm[2*i+1]*gmr[2*i+2] + tmr[2*i+1]*grr[2*i+2];
    cartrm = trm[2*i+1]*gmm[2*i+2] + trr[2*i+1]*grm[2*i+2];
    cartmm = tmr[2*i+1]*grm[2*i+2] + tmm[2*i+1]*gmm[2*i+2];
    sumfr = cartrr * alpha0 + cartmr * alpha1;
    sumfm = cartrm * alpha0 + cartmm * alpha1;
    dx = sumfr + sumfm;
    alpha0 = sumfr / dx;
    alpha1 = sumfm / dx;
    llk2 += log(dx);
  }

  return(-llk1-llk2);
}


// negative naive composite log-likelihood
// [[Rcpp::export]]
double nllk_mrme_naive_cmp(NumericVector &theta, NumericMatrix &data,
			   NumericVector &integrControl) {
  if (is_true(any(theta <= 0))) return(NA_REAL);
  if (theta[2] <= theta[3]) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  double lam1 = theta[0], lam0 = theta[1];
  double pm = 1. / lam1 / (1. / lam1 + 1. / lam0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix z  = data(Range(0, n - 1), Range(1, dim));
  NumericVector
    gmm = g11_mrme(z, tt, theta, integrControl),
    grr = g00_mrme(z, tt, theta, integrControl),
    grm = g01_mrme(z, tt, theta, integrControl),
    gmr = g10_mrme(z, tt, theta, integrControl);

  double llk = 0;

  for(int i = 0; i < n; i++) {
    llk += log(pm * gmm[i] + pm * gmr[i] + pr * grm[i] + pr * grr[i]);
  }

  return(-llk);
}








/***********************************************************
      Following code is for testing purpose only
***********************************************************/

// [[Rcpp::export]]
double nllk_mrme_fixed_sig_err(NumericVector &theta, double sig_err,
			       NumericMatrix &data,
			       NumericVector &integrControl){
  // the theta here only contains lam1, lam0, sigma
  theta.push_back(sig_err);
  return(nllk_mrme(theta, data, integrControl));
}


// [[Rcpp::export]]
double nllk_mrme_one_chain(NumericVector &theta, NumericMatrix &data,
			   NumericVector &integrControl) {
  if (is_true(any(theta <= 0))) return(NA_REAL);
  if (theta[2] <= theta[3]) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  if (n < 2) {
    warning("Sample size is too small to process, should be at least 3. Return nllk as 0.");
    return(0);
  }
  double lam1 = theta[0], lam0 = theta[1];
  double pm = 1. / lam1 / (1. / lam1 + 1. / lam0), pr = 1. - pm;
  NumericVector tt = data.column(0);
  NumericMatrix z  = data(Range(0, n - 1), Range(1, dim));
  NumericVector
    gmm = g11_mrme(z, tt, theta, integrControl),
    grr = g00_mrme(z, tt, theta, integrControl),
    grm = g01_mrme(z, tt, theta, integrControl),
    gmr = g10_mrme(z, tt, theta, integrControl);
  NumericVector
    tmm = t11_mrme(tt, theta),
    trr = t00_mrme(tt, theta),
    trm = t01_mrme(tt, theta),
    tmr = t10_mrme(tt, theta);

  // forward algorithm for the first chain
  // start from Z_0
  double alpha0 = pr, alpha1 = pm;//alpha's are normalized forward variable
  double llk1 = 0;
  double cartrr = 0, cartrm = 0, cartmm = 0, cartmr = 0;
  double dx = 0, sumfr = 0, sumfm = 0;
  for (int i = 0; i < std::floor(static_cast <double> (n) / 2); i++) {
    cartrr = trm[2*i]*gmr[2*i+1] + trr[2*i]*grr[2*i+1];
    cartmr = tmm[2*i]*gmr[2*i+1] + tmr[2*i]*grr[2*i+1];
    cartrm = trm[2*i]*gmm[2*i+1] + trr[2*i]*grm[2*i+1];
    cartmm = tmr[2*i]*grm[2*i+1] + tmm[2*i]*gmm[2*i+1];
    sumfr = cartrr * alpha0 + cartmr * alpha1;
    sumfm = cartrm * alpha0 + cartmm * alpha1;
    dx = sumfr + sumfm;
    alpha0 = sumfr / dx;
    alpha1 = sumfm / dx;
    llk1 += log(dx);
  }

  return(-llk1);
}


// [[Rcpp::export]]
double nllk_mrme_one_chain_fixed_sig_err(NumericVector &theta,
					 double sig_err,
					 NumericMatrix &data,
					 NumericVector &integrControl){
  // the theta here only contains lam1, lam0, sigma
  theta.push_back(sig_err);
  return(nllk_mrme_one_chain(theta, data, integrControl));
}
