#include <R.h>
#include <Rmath.h>


/* pmm(w, t) dw = Pr[M(t) in (w, w + dw), s(0) = 1, s(t) = 1] */
double pmm1(double w, double t, double lamM, double lamR) {
  double lr = lamR * (t - w), lm = lamM * w;
  double dens = 0.0, inc;
  int n = 1;
  while (1) {
    inc = dpois(n, lm, 0) * dpois(n - 1, lr, 0) * lamR;
/*  printf("n = %d, inc = %f\n", n, inc); */
    dens += inc;
    if (inc == 0) break;
    n++;
  }
  return(dens);
}


void pmm(double *w, double *t, double *lamM, double *lamR, int *wlen,
	 double *dens) {
  int i;
  for (i = 0; i < *wlen; i++) {
    dens[i] = pmm1(w[i], t[i], *lamM, *lamR);
  }
}

/* pmr(w, t) dw = Pr[M(t) in (w, w + dw), s(0) = 1, s(t) = 0] */
double pmr1(double w, double t, double lamM, double lamR) {
  double lr = lamR * (t - w), lm = lamM * w;  
  double dens = 0.0, inc;
  int n = 1;
  while (1) {
    inc = dpois(n, lm, 0) * dpois(n, lr, 0) * lamM;
    dens += inc;
    if (inc == 0) break;
    n++;
  }
  dens += lamM * exp(- lm - lr);
  return(dens);
}

void pmr(double *w, double *t, double *lamR, double *lamM, int *wlen,
	 double *dens) {
  int i;
  for (i = 0; i < *wlen; i++) {
    dens[i] = pmr1(w[i], t[i], *lamM, *lamR);
  }
}


/* density of resting for duration w during (0, t), starting from moving */
/* double stayDens1(double w, double t, double lamR, double lamM) { */
/*   double m = 0.0, inc; */
/*   double lr = lamR * (t - w), lm = lamM * w; */
/*   int j = 0; */
/*   double pj = dpois(j, lm, 0), pjp1; */
/*   while (1) { */
/*     pjp1 = dpois(j + 1, lm, 0); */
/*     inc = dpois(j, lr, 0) * ( lamR * pjp1 + lamM * pj ); */
/*     // m = m + inc; */
/*     m += inc; */
/*     // Rprintf("j = %d, m = %f\n", j, m); */
/*     if (inc == 0) break; */
/*     j = j + 1; */
/*     pj = pjp1; */
/*   } */
/*   return(m); */
/* } */

/* void stayDens(double *w, double *t, double *lamR, double *lamM, int *tlen, */
/* 	      double *m) { */
/*   int i; */
/*   for (i = 0; i < *tlen; i++) { */
/*     m[i] = stayDens1(w[i], *t, *lamR, *lamM); */
/*   } */
/* } */


/* Simulation code 
   Note the orders of mm and mr
*/
void staySim(int *n, double *s, double *mm, double *mr,
	     double *t) {
  /* starting from moving state, return simulated time length
     spent in moving state during time (0, s)
     1/mr: mean of exp time for resting
     1/mm: mean of exp time for moving
     t[i]: total moving time in (0, s) for t[i]
  */
  int i, ind;
  double tr, tm;
  for (i = 0; i < *n; i++) {
      tr = 0.0; tm = 0.0;
      while (1) {
	  tm += rexp(*mm); ind = 1; /* time in moving */
	  if (tm + tr > *s) break;
	  tr += rexp(*mr); ind = 2; /* time in resting */
	  if (tm + tr > *s) break;
      }
      if (ind == 1) t[i] = *s - tr; /* last period is moving */
      else t[i] = tm;  /* last period is resting */
  }
}
