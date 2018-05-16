# CHANGES IN coga VERSION 0.2.2

## MAJOR CHANGES

* Rewrite function `dcoga2dim` and `pcoga2dim` to make it faster and more robust.

* Rewrite function `dcoga` to make it can handle wider ranges of parameters and x.

* Add new function `dcoga_approx` and `pcoga_approx` for approximation method, that can impove the speed of code under three or more variables case.

## MINOR CHANGES

* Remove `microbenchmark` from *Suggests* according to the requirement from CRAN.

* Move vignette to inst in order to pass the cran check with devel R.
