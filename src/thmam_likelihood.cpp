// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppGSL, RcppParallel)]]

#include "thmam_likelihood.h"

/***************** basic formula for coga *****************/
double dcoga2dim(double x, double shape1, double shape2,
		       double rate1, double rate2) {
  // transfer rate to scale
  double beta1 = 1 / rate1;
  double beta2 = 1 / rate2;
  
  // gsl_set_error_handler_off();
  // double lgam = shape1 + shape2;
  // double parx = (1/beta1 - 1/beta2) * x;
  // double result = pow(x, lgam - 1);
  // result *= gsl_sf_hyperg_1F1(shape2, lgam, parx);
  // result /= pow(beta1, shape1) * pow(beta2, shape2);
  // result /= exp(R::lgammafn(lgam) + (x / beta1));
  // return result;

  gsl_set_error_handler_off();
  double lgam = shape1 + shape2;
  double parx = (1/beta1 - 1/beta2) * x;
  double result = gsl_sf_hyperg_1F1(shape2, lgam, parx);
  result *= R::dgamma(x, lgam, beta1, 0);
  result *= pow(beta1 / beta2, shape2);
  return result;
}

double pcoga2dim_diff_shape (double x,
			     double shape1, double shape2,
			     double rate1, double rate2) {
  gsl_set_error_handler_off();
  // double result = pow(rate1, shape1) * pow(rate2, shape2);
  // double lgam = shape1 + shape2 + 1;
  // double parx = x * (rate1 - rate2);
  // result *= pow(x, lgam - 1);
  // result /= exp(R::lgammafn(lgam) + (x * rate1));
  // result *= gsl_sf_hyperg_1F1(shape2, lgam, parx);
  // return result;

  double result = pow(rate2 / rate1, shape2);
  double lgam = shape1 + shape2 + 1;
  double parx = x * (rate1 - rate2);
  result /= rate1;
  result *= R::dgamma(x, lgam, 1 / rate1, 0);
  result *= gsl_sf_hyperg_1F1(shape2, lgam, parx);
  return result;
}


/**************************************************************************
 ****************************** calculate p  ******************************
 **************************************************************************/

/****************************** calculate p00 *****************************/
double sumT_p00(double s, double t,
	    double lambda1, double lambda2,
	    double p, int n) {
  double result = 0;
  double cart = pow(1 - p, n);
  for (int k = 0; k < n + 1; ++k) {
    result += dcoga2dim(t - s, k, n - k, lambda1, lambda2) * cart;
    cart *= p * (n - k) / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p00(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = R::pgamma(s, n, 1/lambda0, 1, 0) - R::pgamma(s, n + 1, 1/lambda0, 1, 0);
    cart *= sumT_p00(s, t, lambda1, lambda2, p, n);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0 && cartlast >= cart && n > 1) break;
    cartlast = cart;
    n++;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector ths_vp00(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p00(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}


/****************************** calculate p01 *****************************/

double sumT_p01(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  double pdsx = t - s;
  double cart = pow(1 - p, n);
  for (int k = 0; k < n + 1; ++k) {
    result += pcoga2dim_diff_shape(pdsx, k, n - k, lambda1, lambda2) * cart;
    cart *= (n - k) * p / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p01(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 0;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = p * R::dgamma(s, n + 1, 1/lambda0, 0);
    cart *= sumT_p01(s, t, lambda1, lambda2, p, n);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0 && cartlast >= cart && n > 1) break;
    cartlast = cart;
    n++;
  }

  return result;
}

// [[Rcpp::export]]
NumericVector ths_vp01(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p01(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/****************************** calculate p02 *****************************/

double ths_p02(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  return ths_p01(s, t, lambda0, lambda2, lambda1, 1 - p);
}

// [[Rcpp::export]]
NumericVector ths_vp02(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  return ths_vp01(vs, t, lambda0, lambda2, lambda1, 1 - p);
}


/****************************** calculate p10 *****************************/
double sumT_p10(double s, double t,
	    double lambda1, double lambda2,
	    double p, int n) {
  double result = 0;
  double cart = pow(1 - p, n);
  for (int k = 0; k < n + 1; ++k) {
    result += dcoga2dim(t - s, k + 1, n - k, lambda1, lambda2) * cart;
    cart *= p * (n - k) / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p10(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 0;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = R::pgamma(s, n, 1/lambda0, 1, 0) - R::pgamma(s, n + 1, 1/lambda0, 1, 0);
    cart *= sumT_p10(s, t, lambda1, lambda2, p, n);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0 && cartlast >= cart && n > 1) break;
    cartlast = cart;
    n++;
  }

  return result;
}

// [[Rcpp::export]]
NumericVector ths_vp10(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p10(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}


/****************************** calculate p11 *****************************/

double sumT_p11(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  double pdsx = t - s;
  double cart = pow(1 - p, n - 1);
  for (int k = 0; k < n; ++k) {
    result += pcoga2dim_diff_shape(pdsx, k + 1, n - k - 1, lambda1, lambda2) * cart;
    cart *= (n - k - 1) * p / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p11(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = p * R::dgamma(s, n, 1/lambda0, 0);
    cart *= sumT_p11(s, t, lambda1, lambda2, p, n);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0 && cartlast >= cart && n > 1) break;
    cartlast = cart;
    n++;
  }

  return result;
}

// [[Rcpp::export]]
NumericVector ths_vp11(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p11(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}


/****************************** calculate p12 *****************************/

double sumT_p12(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  double pdsx = t - s;
  double cart = pow(1 - p, n - 1);
  for (int k = 0; k < n; ++k) {
    result += pcoga2dim_diff_shape(pdsx, n - k - 1, k + 1, lambda2, lambda1) * cart;
    cart *= (n - k - 1) * p / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p12(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = (1 - p) * R::dgamma(s, n, 1/lambda0, 0);
    cart *= sumT_p12(s, t, lambda1, lambda2, p, n);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0 && cartlast >= cart && n > 1) break;
    cartlast = cart;
    n++;
  }

  return result;
}

// [[Rcpp::export]]
NumericVector ths_vp12(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p12(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/****************************** calculate p20 *****************************/

double ths_p20(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  return ths_p10(s, t, lambda0, lambda2, lambda1, 1 - p);
}

// [[Rcpp::export]]
NumericVector ths_vp20(NumericVector vs, double t,
			double lambda0, double lambda1, double lambda2,
			double p) {
  return ths_vp10(vs, t, lambda0, lambda2, lambda1, 1 - p);
}


/****************************** calculate p21 *****************************/

double ths_p21(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  return ths_p11(s, t, lambda0, lambda2, lambda1, 1 - p);
}

// [[Rcpp::export]]
NumericVector ths_vp21(NumericVector vs, double t,
			double lambda0, double lambda1, double lambda2,
			double p) {
  return ths_vp11(vs, t, lambda0, lambda2, lambda1, 1 - p);
}

/****************************** calculate p22 *****************************/

double ths_p22(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  return ths_p12(s, t, lambda0, lambda2, lambda1, 1 - p);
}

// [[Rcpp::export]]
NumericVector ths_vp22(NumericVector vs, double t,
			double lambda0, double lambda1, double lambda2,
			double p) {
  return ths_vp12(vs, t, lambda0, lambda2, lambda1, 1 - p);
}


/****************************************************************************
 ************************** prepare for integrate  **************************
 ****************************************************************************/

void ths_f00(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p00(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f01(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p01(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f02(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p02(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f10(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p10(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f11(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p11(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f12(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p12(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f20(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p20(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f21(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p21(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f22(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p22(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

/*****************************************************************************
 ****************************** integrate for h ******************************
 ***************************** serial version ********************************
 *****************************************************************************/

// [[Rcpp::export]]
NumericVector ths_h00(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    // for atom
    double sd = sigma * sqrt(t[i]);
    double prod = exp(-lambda0 * t[i]);
    for (int j = 0; j < dim; j++) {
      ex[7 + j] = x(i, j);
      prod *= R::dnorm(x(i, j), 0.0, sd, 0);
    }
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f00, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result + prod;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector ths_h01(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f01, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h02(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f02, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector ths_h10(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f10, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h11(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f11, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h12(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f12, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h20(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f20, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h21(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f21, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector ths_h22(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = R_Calloc(7 + dim, double);
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

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f22, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  return(value);
}


/*****************************************************************************
 ****************************** integrate for h ******************************
 ****************************** parallel version *****************************
 *****************************************************************************/

// parallel h00 **************************

struct THS_h00_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h00_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      // for atom
      double sd = sigma * sqrt(t[i]);
      double prod = exp(-lambda0 * t[i]);
      for (int j = 0; j < dim; j++) {
        ex[7 + j] = x(i, j);
        prod *= R::dnorm(x(i, j), 0.0, sd, 0);
      }
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f00, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result + prod;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h00_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h00_p ths_h00_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h00_p, grainSize);

  return value;
  
}

// parallel h01 **************************

struct THS_h01_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h01_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f01, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h01_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h01_p ths_h01_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h01_p, grainSize);

  return value;
  
}

// parallel h02 **************************

struct THS_h02_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h02_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f02, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h02_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h02_p ths_h02_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h02_p, grainSize);

  return value;
  
}

// parallel h10 **************************

struct THS_h10_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h10_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f10, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h10_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h10_p ths_h10_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h10_p, grainSize);

  return value;
  
}

// parallel h11 **************************

struct THS_h11_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h11_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f11, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h11_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h11_p ths_h11_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h11_p, grainSize);

  return value;
  
}

// parallel h12 **************************

struct THS_h12_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h12_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f12, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h12_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h12_p ths_h12_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h12_p, grainSize);

  return value;
  
}

// parallel h20 **************************

struct THS_h20_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h20_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f20, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h20_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h20_p ths_h20_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h20_p, grainSize);

  return value;
  
}

// parallel h21 **************************

struct THS_h21_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h21_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f21, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h21_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h21_p ths_h21_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h21_p, grainSize);

  return value;
  
}

// parallel h22 **************************

struct THS_h22_p : public Worker {

  // input
  const RMatrix<double> x;
  const RVector<double> t;
  const RVector<double> theta;
  const RVector<double> integrControl;

  // output
  RVector<double> value;

  // initialize
  THS_h22_p(const NumericMatrix x, const NumericVector t,
	    const NumericVector theta, const NumericVector integrControl,
	    NumericVector value)
    : x(x), t(t), theta(theta), integrControl(integrControl), value(value) {}

  // function call in specified range
  void operator()(std::size_t begin, std::size_t end) {
    // int dim = x.ncol(), n = x.nrow();
    int dim = x.ncol();
    double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
    double sigma = theta[3], p = theta[4];

    /* set up for Rdqags */
    double *ex = R_Calloc(7 + dim, double);
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

    ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
    ex[5] = p; ex[6] = (double) dim;

    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
      // do integrate
      b = t[i]; ex[0] = t[i];
      Rdqags(ths_f22, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	     &limit, &lenw, &last, iwork, work);
      value[i] = result;
    }
  // free memory
  R_Free(ex); R_Free(iwork); R_Free(work);
  }
};


// [[Rcpp::export]]
NumericVector ths_h22_paral(NumericMatrix x, NumericVector t, NumericVector theta,
		            NumericVector integrControl, int grainSize) {
  // allocate return vector
  NumericVector value(x.nrow());

  // creat the worker
  THS_h22_p ths_h22_p(x, t, theta, integrControl, value);

  // call parallelFor
  parallelFor(0, x.nrow(), ths_h22_p, grainSize);

  return value;
  
}


/******************************************************************************
 ********** negative log likelihood function via forward algorithm ************
 ***************************** serial version *********************************
 ******************************************************************************/

// convert vector to matrix
// [[Rcpp::export]]
NumericMatrix con_v_m(NumericVector x) {
  int nx = x.size();
  NumericMatrix y(1, nx);
  y(0, _ ) = x;
  return y;
}

//convert double to vector
// [[Rcpp::export]]
NumericVector con_n_v(double x) {
  NumericVector y(1, x);
  return y;
}

// [[Rcpp::export]]
double nllk_fwd_ths(NumericVector &theta, NumericMatrix &data,
	              NumericVector &integrControl) {
  // theta lambda0, lambda1, lambda2, sigma, p
  // data diff of t and x
  int n = data.nrow(); int dim = data.ncol() - 1;
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double p = theta[4];
  if (lambda1 < lambda2) return NA_REAL;
  double ps0 = 1. / lambda0 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps1 = p / lambda1 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps2 = (1 - p) / lambda2 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  NumericVector tt = data.column(0);
  NumericMatrix x = data(Range(0, n - 1), Range(1, dim));
  double hresult00, hresult01, hresult02;
  double hresult10, hresult11, hresult12;
  double hresult20, hresult21, hresult22;
  double alpha0 = ps0, alpha1 = ps1, alpha2 = ps2;

  // do forward algorithm
  NumericMatrix crowm(1, dim);
  NumericVector ttv(1);
  double llk = 0.;
  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true(all(crow == 0.))) {
      hresult00 = 0.;
      hresult01 = 0.;
      hresult02 = 0.;
      hresult10 = 0.;
      hresult11 = exp(-lambda1 * tt[i]);
      hresult12 = 0.;
      hresult20 = 0.;
      hresult21 = 0.;
      hresult22 = exp(-lambda2 * tt[i]);
    } else {
      crowm = con_v_m(crow);
      ttv   = con_n_v(tt[i]);
      hresult00 = ths_h00(crowm, ttv, theta, integrControl)[0],
      hresult01 = ths_h01(crowm, ttv, theta, integrControl)[0],
      hresult02 = ths_h02(crowm, ttv, theta, integrControl)[0],
      hresult10 = ths_h10(crowm, ttv, theta, integrControl)[0],
      hresult11 = ths_h11(crowm, ttv, theta, integrControl)[0],
      hresult12 = ths_h12(crowm, ttv, theta, integrControl)[0],
      hresult20 = ths_h20(crowm, ttv, theta, integrControl)[0],
      hresult21 = ths_h21(crowm, ttv, theta, integrControl)[0],
      hresult22 = ths_h22(crowm, ttv, theta, integrControl)[0];
    }
    double sumf0 = alpha0 * hresult00 + alpha1 * hresult10 + alpha2 * hresult20;
    double sumf1 = alpha0 * hresult01 + alpha1 * hresult11 + alpha2 * hresult21;
    double sumf2 = alpha0 * hresult02 + alpha1 * hresult12 + alpha2 * hresult22;

    double dx = sumf0 + sumf1 + sumf2;
    alpha0 = sumf0 / dx;
    alpha1 = sumf1 / dx;
    alpha2 = sumf2 / dx;
    llk += log(dx);
  }
  return(-llk);
}



/******************************************************************************
 ********** negative log likelihood function via forward algorithm ************
 *************************** parallel version *********************************
 ******************************************************************************/

// parallel version of likelihood
// [[Rcpp::export]]
double nllk_fwd_ths_parallel(NumericVector &theta, NumericMatrix &data,
	                     NumericVector &integrControl, int grainSize) {
  // theta lambda0, lambda1, lambda2, sigma, p
  // data diff of t and x
  int n = data.nrow(); int dim = data.ncol() - 1;
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double p = theta[4];
  if (lambda1 < lambda2) return NA_REAL;
  double ps0 = 1. / lambda0 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps1 = p / lambda1 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  double ps2 = (1 - p) / lambda2 / (1. / lambda0 + p / lambda1 + (1 - p) / lambda2);
  NumericVector tt = data.column(0);
  NumericMatrix x = data(Range(0, n - 1), Range(1, dim));
  NumericVector
    hresult00 = ths_h00_paral(x, tt, theta, integrControl, grainSize),
    hresult01 = ths_h01_paral(x, tt, theta, integrControl, grainSize),
    hresult02 = ths_h02_paral(x, tt, theta, integrControl, grainSize),
    hresult10 = ths_h10_paral(x, tt, theta, integrControl, grainSize),
    hresult11 = ths_h11_paral(x, tt, theta, integrControl, grainSize),
    hresult12 = ths_h12_paral(x, tt, theta, integrControl, grainSize),
    hresult20 = ths_h20_paral(x, tt, theta, integrControl, grainSize),
    hresult21 = ths_h21_paral(x, tt, theta, integrControl, grainSize),
    hresult22 = ths_h22_paral(x, tt, theta, integrControl, grainSize);
  double alpha0 = ps0, alpha1 = ps1, alpha2 = ps2;

  // do forward algorithm
  double llk = 0.;
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
    double sumf0 = alpha0 * hresult00[i] + alpha1 * hresult10[i] + alpha2 * hresult20[i];
    double sumf1 = alpha0 * hresult01[i] + alpha1 * hresult11[i] + alpha2 * hresult21[i];
    double sumf2 = alpha0 * hresult02[i] + alpha1 * hresult12[i] + alpha2 * hresult22[i];

    double dx = sumf0 + sumf1 + sumf2;
    alpha0 = sumf0 / dx;
    alpha1 = sumf1 / dx;
    alpha2 = sumf2 / dx;
    llk += log(dx);
  }
  return(-llk);
}


