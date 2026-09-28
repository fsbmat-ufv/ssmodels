# Two-Step Estimation of the Classic Heckman Model

Estimates the parameters of the classical Heckman sample selection model
using the two-step procedure.

## Usage

``` r
twostep(selection, outcome, data = sys.frame(sys.parent()))
```

## Arguments

- selection:

  A formula for the selection equation.

- outcome:

  A formula for the outcome equation.

- data:

  A data frame containing the variables.

## Value

A numeric vector containing:

- Estimated coefficients of the selection equation (Probit model),

- Estimated coefficients of the outcome equation (excluding IMR),

- Estimated standard deviation of the outcome errors (`phi`),

- Estimated correlation between the error terms (`cor`).

## Details

The two-step method first estimates a Probit model for the selection
equation, then fits an outcome equation that includes the Inverse Mills
Ratio (IMR) as an additional regressor to correct for sample selection
bias.

## References

There are no references for Rd macro `\insertAllCites` on this help
page. For details, see Heckman (1979) .

## Examples

``` r
data(MEPS2001)
attach(MEPS2001)
#> The following objects are masked from MEPS2001 (pos = 3):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
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
#> The following objects are masked from Mroz87 (pos = 8):
#> 
#>     age, educ
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
#> The following objects are masked from MEPS2001 (pos = 11):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
#> The following objects are masked from Mroz87 (pos = 12):
#> 
#>     age, educ
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
#> The following objects are masked from MEPS2001 (pos = 15):
#> 
#>     age, age2, agefem, ambexp, blhisp, dambexp, dhospexp, educ,
#>     fairpoor, female, ffs, good, hospexp, income, ins, instype,
#>     instype_s1, lambexp, lnambx, totchr, vgood, year01
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
twostep(selectEq, outcomeEq, data = MEPS2001)
#> XS(Intercept)         XSage      XSfemale        XSeduc      XSblhisp 
#>  -0.668643899   0.086814848   0.663505390   0.061883892  -0.365784312 
#>      XStotchr         XSins      XSincome XO(Intercept)         XOage 
#>   0.795747277   0.169106526   0.002677301   5.288927373   0.202466773 
#>      XOfemale        XOeduc      XOblhisp      XOtotchr         XOins 
#>   0.292133967   0.012388871  -0.182865733   0.500633176  -0.046509658 
#>       phi.IMR           IMR 
#>   1.290541875  -0.359316986 
```
