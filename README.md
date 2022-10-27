
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ssmodels <a href='https://fsbmat-ufv.github.io/ssmodels/'><img src='man/figures/logo.png' align="right" height="80" style="margin:10px;" /></a>

# ssmodels

<!-- badges: start -->
<!-- badges: end -->
<p style="text-align: justify;">
In order to facilitate the adjustment of the sample selection models
existing in the literature, we created the ssmodels package. Our package
allows the adjustment of the classic Heckman model (Heckman (1976),
Heckman (1979)), and the estimation of the parameters of this model via
the maximum likelihood method and two-step method, in addition to the
adjustment of the Heckman-t models, introduced in the literature by
Marchenko and Genton (2012) and the Heckman-Skew model introduced in the
literature by Ogundimu and Hutton (2016). We also implemented functions
to adjust the generalized version of the Heckman model, introduced by
Bastos, Barreto-Souza, and Genton (2021), that allows the inclusion of
covariables to the dispersion and correlation parameters and a function
to adjust the Heckman-BS model introduced by Bastos and Barreto-Souza
(2020) that uses the Birnbaum-Saunders distribution as a joint
distribution of the selection and primary regression variables.
</p>

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

## Code of Conduct

`{ssmodels}` is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bastos" class="csl-entry">

Bastos, Fernando de Souza, and Wagner Barreto-Souza. 2020.
“Birnbaum–Saunders Sample Selection Model.” *Journal of Applied
Statistics*.

</div>

<div id="ref-bastosBarreto" class="csl-entry">

Bastos, Fernando de Souza, Wagner Barreto-Souza, and Marc G Genton.
2021. “A Generalized Heckman Model with Varying Sample Selection Bias
and Dispersion Parameters.” *Statistica Sinica*.

</div>

<div id="ref-heckman1976common" class="csl-entry">

Heckman, James J. 1976. “The Common Structure of Statistical Models of
Truncation, Sample Selection and Limited Dependent Variables and a
Simple Estimator for Such Models.” In *Annals of Economic and Social
Measurement, Volume 5, Number 4*, 475–92. NBER.

</div>

<div id="ref-heckman1979sample" class="csl-entry">

———. 1979. “Sample Selection Bias as a Specification Error.”
*Econometrica: Journal of the Econometric Society*, 153–61.

</div>

<div id="ref-marchenko2012heckman" class="csl-entry">

Marchenko, Yulia V, and Marc G Genton. 2012. “A Heckman Selection-t
Model.” *Journal of the American Statistical Association* 107 (497):
304–17.

</div>

<div id="ref-ogundimu2016sample" class="csl-entry">

Ogundimu, Emmanuel O, and Jane L Hutton. 2016. “A Sample Selection Model
with Skew-Normal Distribution.” *Scandinavian Journal of Statistics* 43
(1): 172–90.

</div>

</div>
