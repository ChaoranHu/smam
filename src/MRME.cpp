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
// nllk_mrme = nllk_chain1 + nllk_chain2
// chain1 starts from Z_0
// chain2 starts from Z_1

// [[Rcpp::export]]
double nllk_mrme(NumericVector &theta, NumericMatrix &data,
		 NumericVector &integrControl) {
  if (is_true(any(theta <= 0))) return(NA_REAL);
  if (theta[3] <= theta[4]) return(NA_REAL);
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
    tmm = t11_mrme(z, theta),
    trr = t00_mrme(z, theta),
    trm = t01_mrme(z, theta),
    tmr = t10_mrme(z, theta);

  // forward algorithm for the first chain
  // start from Z_0
  double alpha0 = pr, alpha1 = pm;//alpha's are normalized forward variable
  double llk1 = 0;
  double cartrr = 0, cartrm = 0, cartmm = 0, cartmr = 0;
  double dx = 0, sumfr = 0, sumfm = 0;
  for (int i = 0; i < floor(n/2); i++) {
    cartrr = trm[2*i]*gmr[2*i+1] + trr[2*i]*grr[2*i+1];
    cartmr = tmm[2*i]*gmr[2*i+1] + tmr[2*i]*grr[2*i+1];
    cartrm = trm[2*i]*gmm[2*i+1] + trr[2*i]*grm[2*i+1];
    cartmr = tmr[2*i]*grr[2*i+1] + tmm[2*i]*gmr[2*i+1];
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
  for(int i = 0; i < floor((n-1)/2); i++) {
    cartrr = trm[2*i+1]*gmr[2*i+2] + trr[2*i+1]*grr[2*i+2];
    cartmr = tmm[2*i+1]*gmr[2*i+2] + tmr[2*i+1]*grr[2*i+2];
    cartrm = trm[2*i+1]*gmm[2*i+2] + trr[2*i+1]*grm[2*i+2];
    cartmr = tmr[2*i+1]*grr[2*i+2] + tmm[2*i+1]*gmr[2*i+2];
    sumfr = cartrr * alpha0 + cartmr * alpha1;
    sumfm = cartrm * alpha0 + cartmm * alpha1;
    dx = sumfr + sumfm;
    alpha0 = sumfr / dx;
    alpha1 = sumfm / dx;
    llk2 += log(dx);
  }

  return(-llk1-llk2);
}


// the following code is for testing purpose only
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
  if (theta[3] <= theta[4]) return(NA_REAL);
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
    tmm = t11_mrme(z, theta),
    trr = t00_mrme(z, theta),
    trm = t01_mrme(z, theta),
    tmr = t10_mrme(z, theta);

  // forward algorithm for the first chain
  // start from Z_0
  double alpha0 = pr, alpha1 = pm;//alpha's are normalized forward variable
  double llk1 = 0;
  double cartrr = 0, cartrm = 0, cartmm = 0, cartmr = 0;
  double dx = 0, sumfr = 0, sumfm = 0;
  for (int i = 0; i < floor(n/2); i++) {
    cartrr = trm[2*i]*gmr[2*i+1] + trr[2*i]*grr[2*i+1];
    cartmr = tmm[2*i]*gmr[2*i+1] + tmr[2*i]*grr[2*i+1];
    cartrm = trm[2*i]*gmm[2*i+1] + trr[2*i]*grm[2*i+1];
    cartmr = tmr[2*i]*grr[2*i+1] + tmm[2*i]*gmr[2*i+1];
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
