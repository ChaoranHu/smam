// [[Rcpp::interfaces(r, cpp)]]

#include "movres.h"
#include "thmam_likelihood.h"

// some tools used in this script
NumericVector scale2vector(double x) { // convert scale to vector
  NumericVector x_vec(1, x);
  return(x_vec);
}


NumericMatrix vector2matrix(NumericVector x) { // convert vector to matrix
  NumericMatrix x_mat(1, x.size());
  x_mat(0, _) = x;

  return(x_mat);
}

double myProd(NumericVector x) { // prod in R
  NumericVector result = cumprod(x);
  return result[result.size()-1];
}

NumericMatrix generate_grid(int m, int dim) {
  NumericMatrix result(std::pow(static_cast <double> (m), static_cast <double> (dim)), dim+1);
  NumericVector h(dim+1, 1.0);
  for (int i = 0; i < result.nrow(); i++) {
    result(i, _) = h;
    h[dim - 1] = h[dim - 1] + 1;
    for (int j = dim - 1; j >= 0; j--) {
      if (h[j] > m) {
	h[j] = 1;
	h[j - 1] = h[j - 1] + 1;
      }
    }
  }
  return(result);
}

LogicalVector weak_equal (NumericVector x, NumericVector y) {
  NumericVector diff = x - y;
  diff = abs(diff);
  return(diff == 0.0);
}


// q function: Pr(Z, S(t)=j, \epsilon*_t = x_l | S(0) = i, \epsilon*_0 = y_k)
// theta: lam1, lam0, sigma, sig_err
// err_start, err_start_prob(no needed): \epsilon*_0 and corresponding probability
// err_end, err_end_prob:                \epsilon*_t and corresponding probability
// the length of all vector type parameters is equal to the dimension of data.


// [[Rcpp::export]]
double q10_mrme_approx(NumericVector z, double t, NumericVector theta,
		       NumericVector integrControl,
		       NumericVector err_start, NumericVector err_end,
		       NumericVector err_end_prob) {
  NumericVector h_w = z + err_start - err_end;
  NumericVector zero_cart(z.length());
  LogicalVector zero_ind = weak_equal(h_w, zero_cart);
  
  if (is_true(all(zero_ind))) {
    // Rcout << "I am zero !!!!!!!!!!!" << "\n";
    return(0.);
  } else {
    NumericMatrix h_w_mat = vector2matrix(h_w);
    NumericVector t_vec = scale2vector(t);
    NumericVector h_result = h10(h_w_mat, t_vec, theta[Range(0, 2)], integrControl);
    return(h_result[0] * myProd(err_end_prob)); 
  }
}


// [[Rcpp::export]]
double q01_mrme_approx(NumericVector z, double t, NumericVector theta,
		       NumericVector integrControl,
		       NumericVector err_start, NumericVector err_end,
		       NumericVector err_end_prob) {
  NumericVector h_w = z + err_start - err_end;
  NumericVector zero_cart(z.length());
  LogicalVector zero_ind = weak_equal(h_w, zero_cart);
  
  if (is_true(all(zero_ind))) {
    return(0.);
  } else {
    NumericMatrix h_w_mat = vector2matrix(h_w);
    NumericVector t_vec = scale2vector(t);
    NumericVector h_result = h01(h_w_mat, t_vec, theta[Range(0, 2)], integrControl);
    return(h_result[0] * myProd(err_end_prob));
  }
}


// [[Rcpp::export]]
double q00_mrme_approx(NumericVector z, double t, NumericVector theta,
		       NumericVector integrControl,
		       NumericVector err_start, NumericVector err_end,
		       NumericVector err_end_prob) {
  NumericVector h_w = z + err_start - err_end;
  NumericVector zero_cart(z.length());
  LogicalVector zero_ind = weak_equal(h_w, zero_cart);
  
  if (is_true(all(zero_ind))) {
    return(exp(-theta[1] * t) * myProd(err_end_prob));
  } else {
    NumericMatrix h_w_mat = vector2matrix(h_w);
    NumericVector t_vec = scale2vector(t);
    NumericVector h_result = h00(h_w_mat, t_vec, theta[Range(0, 2)], integrControl);
    return(h_result[0] * myProd(err_end_prob));
  }
}


// [[Rcpp::export]]
double q11_mrme_approx(NumericVector z, double t, NumericVector theta,
		       NumericVector integrControl,
		       NumericVector err_start, NumericVector err_end,
		       NumericVector err_end_prob) {
  NumericVector h_w = z + err_start - err_end;
  NumericVector zero_cart(z.length());
  LogicalVector zero_ind = weak_equal(h_w, zero_cart);
  
  if (is_true(all(zero_ind))) {
    return(0.);
  } else {
    NumericMatrix h_w_mat = vector2matrix(h_w);
    NumericVector t_vec = scale2vector(t);
    NumericVector h_result = h11(h_w_mat, t_vec, theta[Range(0, 2)], integrControl);
    return(h_result[0] * myProd(err_end_prob));
  }
}


// calculate the nllk of mrme model with approximated guassian error
// approx_norm_even & approx_norm_odd: the matrices specify the discrete
//           distribution used to approximate STANDARD normal.
// theta: lam1 lam0 sigma sig_err
// MENTION: Suppose discrete distribution is Pr(X=i)=p_i. Then, the error
//          is distributed as Pr(epsilon* = i*sig_err) = p_i. 


// [[Rcpp::export]]
double nllk_mrme_approx(NumericVector &theta, NumericMatrix &data,
			NumericVector &integrControl,
			NumericMatrix &approx_norm_even,
			NumericMatrix &approx_norm_odd) {
  if (is_true(any(theta <= 0))) return(NA_REAL);
  // if (theta[2] < theta[3]) return(NA_REAL);
  int n = data.nrow(), dim = data.ncol() - 1;
  int m_even = approx_norm_even.nrow(), m_odd = approx_norm_odd.nrow();
  NumericVector tt = data.column(0);
  NumericMatrix x  = data(Range(0, n-1), Range(1, dim));

  NumericVector approx_norm_even_value = approx_norm_even(_, 0) * theta[3];
  NumericVector approx_norm_odd_value  = approx_norm_odd(_, 0)  * theta[3];
  NumericVector approx_norm_even_prob  = approx_norm_even(_, 1);
  NumericVector approx_norm_odd_prob   = approx_norm_odd(_, 1);
  

  // initial forward variable matrix
  NumericMatrix fwd_mov_even(std::pow(static_cast <double> (m_even), static_cast <double> (dim)), dim+1);
  NumericMatrix fwd_res_even(std::pow(static_cast <double> (m_even), static_cast <double> (dim)), dim+1);
  NumericMatrix fwd_mov_odd(std::pow(static_cast <double> (m_odd), static_cast <double> (dim)), dim+1);
  NumericMatrix fwd_res_odd(std::pow(static_cast <double> (m_odd), static_cast <double> (dim)), dim+1);

  fwd_mov_even = generate_grid(m_even, dim);
  fwd_res_even = generate_grid(m_even, dim);
  fwd_mov_odd  = generate_grid(m_odd, dim);
  fwd_res_odd  = generate_grid(m_odd, dim);


  /*
    fwd_mov_even matrix holds the current forward variables
    that are moving part and even step (t0, t2, t4, ...)
    Suppose we have 2 dimensional data and the discrete
    distribution take m possible values, this matrix is
    m^2 by 3. It looks like following
    1stDimInd  2ndDimInd  fwdValue
    1          1         x
    1          2         x
    ...        ...         x
    m          m         x
  */

  
  double lam1 = theta[0], lam0 = theta[1];
  double pm = 1. / lam1 / (1. / lam1 + 1. / lam0), pr = 1. - pm;

  for (int i = 0; i < fwd_mov_even.nrow(); i++) {
    double cart1 = pm;
    for (int j = 0; j < dim; j++) {
      cart1 *= approx_norm_even(fwd_mov_even(i, j)-1, 1);
    }
    fwd_mov_even(i, dim) = cart1;
  }

  for (int i = 0; i < fwd_res_even.nrow(); i++) {
    double cart1 = pr;
    for (int j = 0; j < dim; j++) {
      cart1 *= approx_norm_even(fwd_res_even(i, j)-1, 1);
    }
    fwd_res_even(i, dim) = cart1;
  }
  // END initial forward variable matrix
  
  double dx = 0;
  double llk = 0;

  // PROCESS forward algorithm
  for (int k = 0; k < n; k++) { // k = row num of data

    NumericVector this_x = x(k, _);
    // Rcout << "k is :" << k << "\n";

    if (k % 2 == 0) {
    
      for (int i = 0; i < fwd_mov_odd.nrow(); i++) { // i = end point ind
	
	NumericVector end_ind = fwd_mov_odd(i, _);
	end_ind = end_ind[Range(0, dim-1)] - 1;
	NumericVector end_error = approx_norm_odd_value[end_ind];
	NumericVector end_prob  = approx_norm_odd_prob[end_ind];

	double cart2 = 0;
	for (int j = 0; j < fwd_mov_even.nrow(); j++) { // j = start point ind

	  NumericVector start_ind = fwd_mov_even(j, _);
	  start_ind = start_ind[Range(0, dim-1)] - 1;
	  NumericVector start_error = approx_norm_even_value[start_ind];
	
	  cart2 += fwd_mov_even(j, dim) * q11_mrme_approx(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);

	  cart2 += fwd_res_even(j, dim) * q01_mrme_approx(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  
	}

	fwd_mov_odd(i, dim) = cart2;
      }

      for (int i = 0; i < fwd_res_odd.nrow(); i++) { // i = end point ind

	NumericVector end_ind = fwd_res_odd(i, _);
	end_ind = end_ind[Range(0, dim-1)] - 1;
	NumericVector end_error = approx_norm_odd_value[end_ind];
	NumericVector end_prob  = approx_norm_odd_prob[end_ind];

	double cart2 = 0;
	for (int j = 0; j < fwd_mov_even.nrow(); j++) { // j = start point ind

	  NumericVector start_ind = fwd_mov_even(j, _);
	  start_ind = start_ind[Range(0, dim-1)] - 1;
	  NumericVector start_error = approx_norm_even_value[start_ind];

	  // Rcout << "i, j is" << i << j <<"\n";
	  // NumericVector test;
	  // test = this_x + start_error - end_error;
	  // Rcout << test << "\n";
	
	  cart2 += fwd_mov_even(j, dim) * q10_mrme_approx(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);

	  cart2 += fwd_res_even(j, dim) * q00_mrme_approx(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  
	}

	fwd_res_odd(i, dim) = cart2;
      }

      dx = sum(fwd_mov_odd(_, dim)) + sum(fwd_res_odd(_, dim));
      fwd_mov_odd(_, dim) = fwd_mov_odd(_, dim) / dx;
      fwd_res_odd(_, dim) = fwd_res_odd(_, dim) / dx;
      llk += log(dx);

    } else {
      
      for (int i = 0; i < fwd_mov_even.nrow(); i++) { // i = end point ind

	NumericVector end_ind = fwd_mov_even(i, _);
	end_ind = end_ind[Range(0, dim-1)] - 1;
	NumericVector end_error = approx_norm_even_value[end_ind];
	NumericVector end_prob  = approx_norm_even_prob[end_ind];

	double cart2 = 0;
	for (int j = 0; j < fwd_mov_odd.nrow(); j++) { // j = start point ind

	  NumericVector start_ind = fwd_mov_odd(j, _);
	  start_ind = start_ind[Range(0, dim-1)] - 1;
	  NumericVector start_error = approx_norm_odd_value[start_ind];
	
	  cart2 += fwd_mov_odd(j, dim) * q11_mrme_approx(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);

	  cart2 += fwd_res_odd(j, dim) * q01_mrme_approx(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  
	}

	fwd_mov_even(i, dim) = cart2;
      }

      for (int i = 0; i < fwd_res_even.nrow(); i++) { // i = end point ind

	NumericVector end_ind = fwd_res_even(i, _);
	end_ind = end_ind[Range(0, dim-1)] - 1;
	NumericVector end_error = approx_norm_even_value[end_ind];
	NumericVector end_prob  = approx_norm_even_prob[end_ind];

	double cart2 = 0;
	for (int j = 0; j < fwd_mov_odd.nrow(); j++) { // j = start point ind

	  NumericVector start_ind = fwd_mov_odd(j, _);
	  start_ind = start_ind[Range(0, dim-1)] - 1;
	  NumericVector start_error = approx_norm_odd_value[start_ind];

	  // Rcout << "i, j is" << i << j <<"\n";
	  // NumericVector test;
	  // test = this_x + start_error - end_error;
	  // Rcout << test << "\n";
	
	  cart2 += fwd_mov_odd(j, dim) * q10_mrme_approx(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);

	  cart2 += fwd_res_odd(j, dim) * q00_mrme_approx(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  
	}

	fwd_res_even(i, dim) = cart2;
      }

      dx = sum(fwd_mov_even(_, dim)) + sum(fwd_res_even(_, dim));
      fwd_mov_even(_, dim) = fwd_mov_even(_, dim) / dx;
      fwd_res_even(_, dim) = fwd_res_even(_, dim) / dx;
      llk += log(dx);
      
    }
  }

  return(-llk);
}



/******************************************

   calculate nllk of one dimentional mrme
   model with approximated error

*******************************************/

/******************************************
   TODO: use weak equal to check zero!!!
 ******************************************/




// q function for 1dim case

// [[Rcpp::export]]
double q10_mrme_approx_1dim(double z, double t, NumericVector theta,
			    NumericVector integrControl,
			    double err_start, double err_end,
			    double err_end_prob) {
  double h_w = z + err_start - err_end;
  if (h_w == 0.) {
    return(0.);
  } else {
    NumericMatrix h_w_mat = vector2matrix(scale2vector(h_w));
    NumericVector t_vec   = scale2vector(t);

    NumericVector h_result = h10(h_w_mat, t_vec, theta[Range(0, 2)],
				 integrControl);
    return(h_result[0] * err_end_prob);
  }
}

// [[Rcpp::export]]
double q01_mrme_approx_1dim(double z, double t, NumericVector theta,
			    NumericVector integrControl,
			    double err_start, double err_end,
			    double err_end_prob) {
  double h_w = z + err_start - err_end;
  if (h_w == 0.) {
    return(0.);
  } else {
    NumericMatrix h_w_mat = vector2matrix(scale2vector(h_w));
    NumericVector t_vec   = scale2vector(t);

    NumericVector h_result = h01(h_w_mat, t_vec, theta[Range(0, 2)],
				 integrControl);
    return(h_result[0] * err_end_prob);
  }
}

// [[Rcpp::export]]
double q00_mrme_approx_1dim(double z, double t, NumericVector theta,
			    NumericVector integrControl,
			    double err_start, double err_end,
			    double err_end_prob) {
  double h_w = z + err_start - err_end;
  if (h_w == 0.) {
    return(exp(-theta[1] * t) * err_end_prob);
  } else {
    NumericMatrix h_w_mat = vector2matrix(scale2vector(h_w));
    NumericVector t_vec   = scale2vector(t);

    NumericVector h_result = h00(h_w_mat, t_vec, theta[Range(0, 2)],
				 integrControl);
    return(h_result[0] * err_end_prob);
  }
}

// [[Rcpp::export]]
double q11_mrme_approx_1dim(double z, double t, NumericVector theta,
			    NumericVector integrControl,
			    double err_start, double err_end,
			    double err_end_prob) {
  double h_w = z + err_start - err_end;
  if (h_w == 0.) {
    return(0.);
  } else {
    NumericMatrix h_w_mat = vector2matrix(scale2vector(h_w));
    NumericVector t_vec   = scale2vector(t);

    NumericVector h_result = h11(h_w_mat, t_vec, theta[Range(0, 2)],
				 integrControl);
    return(h_result[0] * err_end_prob);
  }
}



// [[Rcpp::export]]
double nllk_mrme_approx_1dim(NumericVector &theta, NumericMatrix &data,
			     NumericVector &integrControl,
			     NumericMatrix &approx_norm_even,
			     NumericMatrix &approx_norm_odd) {
  if (is_true(any(theta <= 0))) return(NA_REAL);
  // if (theta[2] < theta[3]) return(NA_REAL);
  int n = data.nrow(); // dim = 1;
  int m_even = approx_norm_even.nrow(), m_odd = approx_norm_odd.nrow();
  NumericVector tt = data.column(0);
  NumericVector x  = data.column(1);

  NumericVector approx_norm_even_value = approx_norm_even(_, 0) * theta[3];
  NumericVector approx_norm_odd_value  = approx_norm_odd(_, 0)  * theta[3];
  NumericVector approx_norm_even_prob  = approx_norm_even(_, 1);
  NumericVector approx_norm_odd_prob   = approx_norm_odd(_, 1);

  // initial forward variable matrix
  NumericMatrix fwd_mov_even(m_even, 2);
  NumericMatrix fwd_res_even(m_even, 2);
  NumericMatrix fwd_mov_odd(m_odd, 2);
  NumericMatrix fwd_res_odd(m_odd, 2);

  fwd_mov_even(_, 0) = Range(1, m_even);
  fwd_res_even(_, 0) = Range(1, m_even);
  fwd_mov_odd(_, 0)  = Range(1, m_odd);
  fwd_res_odd(_, 0)  = Range(1, m_odd);

  double lam1 = theta[0], lam0 = theta[1];
  double pm = 1. / lam1 / (1. / lam1 + 1. / lam0), pr = 1. - pm;

  for (int i = 0; i < m_even; i++) {
    fwd_mov_even(i, 1) = pm * approx_norm_even_prob[i];
    fwd_res_even(i, 1) = pr * approx_norm_even_prob[i];
  }

  // END initial forward variable matrix

  double dx = 0, llk = 0;

  // PROCESS forward algorithm
  for (int k = 0; k < n; k++) { // k = row num of data

    // Rcout << "k is :" << k << "\n";
    
    double this_x = x[k];

    if (k % 2 == 0) {
      // Rcout << "even to odd" << "\n";
      for (int i = 0; i < m_odd; i++) { // i = end point ind

	double end_error = approx_norm_odd_value[i];
	double end_prob  = approx_norm_odd_prob[i];

	double cart2 = 0;
	for (int j = 0; j < m_even; j++) { // j = start point ind
	  
	  double start_error = approx_norm_even_value[j];
	  cart2 += fwd_mov_even(j, 1) * q11_mrme_approx_1dim(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  cart2 += fwd_res_even(j, 1) * q01_mrme_approx_1dim(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	}
	fwd_mov_odd(i, 1) = cart2;
	
      }

      for (int i = 0; i < m_odd; i++) { // i = end point ind

	double end_error = approx_norm_odd_value[i];
	double end_prob  = approx_norm_odd_prob[i];

	double cart2 = 0;
	for (int j = 0; j < m_even; j++) { // j = start point ind
	  
	  double start_error = approx_norm_even_value[j];
	  cart2 += fwd_mov_even(j, 1) * q10_mrme_approx_1dim(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  cart2 += fwd_res_even(j, 1) * q00_mrme_approx_1dim(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  
	}
	fwd_res_odd(i, 1) = cart2;
	
      }

      dx = sum(fwd_mov_odd(_, 1)) + sum(fwd_res_odd(_, 1));
      fwd_mov_odd(_, 1) = fwd_mov_odd(_, 1) / dx;
      fwd_res_odd(_, 1) = fwd_res_odd(_, 1) / dx;
      llk += log(dx);

    } else {

      // Rcout << "odd to even" << "\n";
      for (int i = 0; i < m_even; i++) { // i = end point ind

	double end_error = approx_norm_even_value[i];
	double end_prob  = approx_norm_even_prob[i];

	double cart2 = 0;
	for (int j = 0; j < m_odd; j++) { // j = start point ind

	  double start_error = approx_norm_odd_value[j];
	  cart2 += fwd_mov_odd(j, 1) * q11_mrme_approx_1dim(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  cart2 += fwd_res_odd(j, 1) * q01_mrme_approx_1dim(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  
	}
	fwd_mov_even(i, 1) = cart2;
      }


      for (int i = 0; i < m_even; i++) { // i = end point ind

	double end_error = approx_norm_even_value[i];
	double end_prob  = approx_norm_even_prob[i];

	double cart2 = 0;
	for (int j = 0; j < m_odd; j++) { // j = start point ind

	  double start_error = approx_norm_odd_value[j];
	  cart2 += fwd_mov_odd(j, 1) * q10_mrme_approx_1dim(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  cart2 += fwd_res_odd(j, 1) * q00_mrme_approx_1dim(this_x, tt[k], theta, integrControl, start_error, end_error, end_prob);
	  
	}
	fwd_res_even(i, 1) = cart2;
      }

      dx = sum(fwd_mov_even(_, 1)) + sum(fwd_res_even(_, 1));
      fwd_mov_even(_, 1) = fwd_mov_even(_, 1) / dx;
      fwd_res_even(_, 1) = fwd_res_even(_, 1) / dx;
      llk += log(dx);
    }
  }

  return(-llk);
}

/*

Check whether nllk_mrme_approx and nllk_mrme_approx_1dim are consistent.

smam:::nllk_mrme_approx_1dim(c(1,1,1,0.1), matrix(1:8, ncol = 2), c(1,1,1), cbind(c(-1, 0, 1), c(0.3, 0.4, 0.3)), cbind(c(-5, -2.5, 2.5, 5), c(0.25, 0.25, 0.25, 0.25)))

smam:::nllk_mrme_approx(c(1,1,1,0.1), matrix(1:8, ncol = 2), c(1,1,1), cbind(c(-1, 0, 1), c(0.3, 0.4, 0.3)), cbind(c(-5, -2.5, 2.5, 5), c(0.25, 0.25, 0.25, 0.25)))

*/
