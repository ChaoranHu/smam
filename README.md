# README #

[![Build Status](https://travis-ci.org/ChaoranHu/smam.svg?branch=master)](https://travis-ci.org/ChaoranHu/smam) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/smam)](https://cran.r-project.org/package=smam)


This README would normally document whatever steps are necessary to get this application up and running.

### What is this R package for? ###

* Statistical modeling of animal movements in R. Three statistcal models are provided by this package. They are moving-resting model with embedded Brownian motion^[1,2], moving-resting-hunting model with embedded Brownian motion^[3] and Brownian motion with measurement error^[4]. A quick guide for this package is given as package vignette.

* Version: 0.4.0

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

Note: The first way need your computer have development tools, but can help you follow the newest version from my github. This package includes Cpp and C codes, so you need a CPP compiler (for mac, you can use Xcode, which can be installed from Apple Store). Also, this package use GNU GSL, so you also need to install GNU GSL first. If you do not have these tools, please use second way to install from CRAN.

* Run `library(smam)` to load in R.

### How do I use? ###

Please read the vignette of this package for more information, documentation and examples.

### Who do I talk to? ###

Chaoran Hu, <chaoran.hu@uconn.edu>


### Reference ###

[1] Yan, J., Chen, Y., Lawrence-Apfel, K., Ortega, I. M., Pozdnyakoc, V., Williams, S., and Meyer, T. (2014) A moving-resting process with an embedded Brownian motion for animal movements. Population Ecology. 56(2): 401--415. doi:10.1007/s10144-013-0428-8

[2] Pozdnyakov, V., Elbroch, L., Labarga, A., Meyer, T., and Yan, J. (2017) Discretely observed Brownian motion governed by telegraph process: estimation. Methodology and Computing in Applied Probability. doi:10.1007/s11009-017-9547-6.

[3] Pozdnyakov, V., Elbroch, L., Hu, C., Meyer, T., and Yan, J. (2018+) On estimation for Brownian motion governed by telegraph process with multiple off states. Under Review.

[4] Pozdnyakov, V. , Meyer, T. , Wang, Y. and Yan, J. (2014), On modeling animal movements using Brownian motion with measurement error. Ecology, 95: 247-253. doi:10.1890/13-0532.1
