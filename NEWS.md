# CHANGES IN smam VERSION 0.5.4

## MINOR CHANGES

* Modified `estVarMRME_Godambe` and `estVarMRME_pBootstrap`.

* Added 'f109' data.


# CHANGES IN smam VERSION 0.5.3, 0.5.2, 0.5.1

## MINOR CHANGES

* Fixed some bugs in solaris system. It is difficult to debug in solaris system. CRAN check packages
in solaris system but I cannot do it locally. You can use R-hub to check your own package under
solaris environment but R-hub does not support GSL, which is needed by `smam`. So, this is why
I submitted such many versions.


# CHANGES IN smam VERSION 0.5.0

## MAJOR CHANGES

* Added moving-resting model with measurement error is modeled as Gaussian noise as `fitMRME`.

* Added estimator variance function for `fitMRME` as `estVarMRME_Godambe` and `estVarMRME_pBootstrap`.



# CHANGES IN smam VERSION 0.4.0

## MAJOR CHANGES

* Added moving-resting-handling model with both full likelihood and composite likelihood.

* Provided option to use only part of dataset to fit BMME model, moving-resting model
and moving-resting-handling model (in function `fitBMME`, `fitMR`, and `fitMRH`).

* Provided tools for seasonal behavior analysis.

* Parallel code is implemented for moving-resting-handling model to speed code up.

* Added code for predicting state at given time point under MR model and MRH model (beta version).

## MINOR CHANGES

* Added vignette for quick start.

* Redocumented whole package.

* Renamed `rBmme -> rBMME`, `fitBmme -> fitBMME`, `rMovRes -> rMR` and `fitMovRes -> fitMR`.

* Provided option to show simulated states in `rMR` and `rMRH`.



# CHANGES IN smam VERSION 0.3.0

## MAJOR CHANGES

* Added full likelihood based hidden Markov model

* Full likelihood and composite likelihood are now done with Rcpp




# CHANGES IN smam VERSION 0.2.2

## MINOR CHANGES

* Added Depends: R (>= 3.0.0).

* Corrected double to void of C function pmr.




# CHANGES IN smam VERSION 0.2.1

## MAJOR CHANGES

* First public release.

