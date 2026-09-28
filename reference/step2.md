# Heckman's Two-Step Method

Estimates the parameters of the classical Heckman selection model using
the two-step method. The first step fits a probit model for the
selection equation. In the second step, the inverse Mills ratio (IMR) is
included as an additional regressor in the outcome equation.

## Usage

``` r
step2(YS, XS, YO, XO)
```

## Arguments

- YS:

  A binary vector indicating selection (`1` if observed, `0` otherwise).

- XS:

  A matrix of covariates for the selection equation.

- YO:

  A numeric vector representing the outcome variable of interest.

- XO:

  A matrix of covariates for the outcome equation.

## Value

A numeric vector containing the parameter estimates from the two-step
Heckman model:

- Coefficients of the selection equation (probit model).

- Coefficients of the outcome equation (excluding the IMR term).

- Estimated `sigma`.

- Estimated `rho`.

## Details

This function implements the two-step estimation procedure of the
classical Heckman model. In the first step, a probit model is estimated
to predict the selection indicator `YS` using the selection covariates
`XS`. The IMR is calculated from this model. In the second step, an
ordinary least squares (OLS) regression of the observed outcome `YO` on
`XO` and the IMR is performed for the uncensored observations
(`YS == 1`).

The function also calculates:

- `sigma`: The estimated standard deviation of the outcome equation's
  error term.

- `rho`: The estimated correlation between the error terms of the
  selection and outcome equations.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## Examples

``` r
data(MEPS2001)
attach(MEPS2001)
#> The following objects are masked from data:
#> 
#>     age, educ, income
#> The following objects are masked from data3:
#> 
#>     age, educ, income
#> The following objects are masked from nhanes:
#> 
#>     age, educ, income
#> The following objects are masked from PSID2:
#> 
#>     age, educ
#> The following objects are masked from Mroz87 (pos = 7):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 8):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
#> The following objects are masked from MEPS2001 (pos = 9):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
#> The following objects are masked from MEPS2001 (pos = 10):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
#> The following objects are masked from Mroz87 (pos = 11):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 12):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
#> The following objects are masked from MEPS2001 (pos = 13):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
#> The following objects are masked from MEPS2001 (pos = 14):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
YS <- dambexp
XS <- cbind(age, female, educ, blhisp, totchr, ins)
YO <- lnambx
XO <- cbind(age, female, educ, blhisp, totchr, ins, income)
step2(YS, XS, YO, XO)
#>        XSage     XSfemale       XSeduc     XSblhisp     XStotchr        XSins 
#>  0.041950311  0.629046530  0.035839783 -0.458985308  0.785774049  0.182828788 
#>        XOage     XOfemale       XOeduc     XOblhisp     XOtotchr        XOins 
#>  0.450637479  1.262172664  0.182049809 -0.649939053  1.134936458  0.247641171 
#>     XOincome        sigma          rho 
#> -0.002021763  2.620419828  1.639810343 
```
