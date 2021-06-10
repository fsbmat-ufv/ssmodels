
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ssmodels

<!-- badges: start -->
<!-- badges: end -->

In order to facilitate the adjustment of the sample selection models
existing in the literature, we created the ssmodels package. Our package
allows the adjustment of the classic Heckman model, and the estimation
of the parameters of this model via the maximum likelihood method and
two-step method, in addition to the adjustment of the Heckman-t models,
introduced in the literature in 2012 by @marchenko2012heckman and the
Heckman-Skew model introduced in the literature by @ogundimu2016sample
in 2016. We also implemented functions to adjust the generalized version
of the Heckman model, introduced by @bastosBarreto, that allows the
inclusion of covariables to the dispersion and correlation parameters
and a function to adjust the Heckman-BS model introduced in 2020 by
@bastos that uses the Birnbaum-Saunders distribution as a joint
distribution of the selection and primary regression variables.

## Installation

You can install the released version of ssmodels from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ssmodels")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fsbmat-ufv/ssmodels")
```
