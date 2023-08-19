# README #


[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/smam)](https://cran.r-project.org/package=smam)
[![R-CMD-check](https://github.com/ChaoranHu/smam/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChaoranHu/smam/actions/workflows/R-CMD-check.yaml)

This README would normally document whatever steps are necessary to get this application up and running.

### What is this R package for? ###

* Statistical modeling of animal movements in R. Statistcal models in reference are provided by this package. They are moving-resting model with embedded Brownian motion^[1,2], moving-resting-handling model with embedded Brownian motion^[3], Brownian motion with measurement error^[4], moving-resting model with measurement error^[5], moving-moving model with two embedded Brownian motions. A quick guide for this package is given as package vignette.

* Version: 0.7.0

### How do I get set up? ###

* Run the following R code.

```
## install.packages("devtools")
devtools::install_github("ChaoranHu/smam")
```

or

```
install.packages("smam")
```

Note: The first way allows you to access to the most recent version of this package but it also requires development toolbox. This package includes Cpp and C codes, so you need a CPP compiler (for mac, you can use Xcode, which can be installed from Apple Store). Also, this package uses GNU GSL, so you need to install GNU GSL first. If you do not have these tools, you could use second way to install from CRAN.

* Run `library(smam)` to load in R.

### How do I use? ###

Please read the vignette of this package for more information, documentation and examples.

### Who do I talk to? ###

Chaoran Hu, <huchaoran.stat@gmail.com>


### Reference ###

[1] Yan, J., Chen, Y., Lawrence-Apfel, K., Ortega, I. M., Pozdnyakoc, V., Williams, S., and Meyer, T. (2014) A moving-resting process with an embedded Brownian motion for animal movements. Population Ecology. 56(2): 401--415. (doi:10.1007/s10144-013-0428-8)

[2] Pozdnyakov, V., Elbroch, L., Labarga, A., Meyer, T., and Yan, J. (2017) Discretely observed Brownian motion governed by telegraph process: estimation. Methodology and Computing in Applied Probability, 21: 907--920. (doi:10.1007/s11009-017-9547-6)

[3] Pozdnyakov, V., Elbroch, L.M., Hu, C. et al. On Estimation for Brownian Motion Governed by Telegraph Process with Multiple Off States. Methodol Comput Appl Probab, 22, 1275–1291 (2020). (doi:10.1007/s11009-020-09774-1)

[4] Pozdnyakov, V. , Meyer, T. , Wang, Y. and Yan, J. (2014), On modeling animal movements using Brownian motion with measurement error. Ecology, 95: 247-253. (doi:10.1890/13-0532.1)

[5] Hu, C., Elbroch, L.M., Meyer, T., Pozdnyakov, V. and Yan, J. (2021), Moving-resting process with measurement error in animal movement modeling. Methods in Ecology and Evolution. (doi:10.1111/2041-210X.13694)
