# Two-Step Method for Parameter Estimation of the Classical Heckman Model

Estimates the parameters of the classical Heckman sample selection model
using the two-step estimation method.

## Usage

``` r
HCinitial(selection, outcome, data = sys.frame(sys.parent()))
```

## Arguments

- selection:

  A formula specifying the selection equation.

- outcome:

  A formula specifying the outcome equation.

- data:

  A data frame containing the variables in the model.

## Value

A named numeric vector containing:

- Coefficients from the selection equation (probit model),

- Coefficients from the outcome equation (excluding the IMR),

- Estimated `sigma`,

- Estimated `rho`.

## Details

This function implements the two-step approach proposed by Heckman
(1979) to estimate the parameters of the classic sample selection model.
It is particularly useful for obtaining initial values for maximum
likelihood estimation (MLE).

In the first step, a probit model is fitted to the selection equation to
estimate the probability of selection. The second step involves
estimating a linear regression of the outcome equation for the observed
(selected) data, incorporating the inverse Mills ratio (IMR) as an
additional regressor to correct for sample selection bias.

The function also estimates:

- `sigma`: The standard deviation of the outcome equation's error term.

- `rho`: The correlation coefficient between the errors of the selection
  and outcome equations.

## References

James J Heckman (1979). “Sample selection bias as a specification
error.” *Econometrica: Journal of the econometric society*, 153–161.

## Examples

``` r
data(MEPS2001)
attach(MEPS2001)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
HCinitial(selectEq, outcomeEq, data = MEPS2001)
#> xs(Intercept)         xsage      xsfemale        xseduc      xsblhisp 
#>  -0.668643899   0.086814848   0.663505390   0.061883892  -0.365784312 
#>      xstotchr         xsins      xsincome xo(Intercept)         xoage 
#>   0.795747277   0.169106526   0.002677301   5.288927373   0.202466773 
#>      xofemale        xoeduc      xoblhisp      xototchr         xoins 
#>   0.292133967   0.012388871  -0.182865733   0.500633176  -0.046509658 
#>         sigma           rho 
#>   1.290541875  -0.359316986 
```
