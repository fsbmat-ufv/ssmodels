# Summary of Heckman-t Model

Prints a detailed summary of the parameter estimates and model fit
statistics for an object of class `HeckmantS`.

## Usage

``` r
# S3 method for class 'HeckmantS'
summary(object, ...)
```

## Arguments

- object:

  An object of class `HeckmantS`, containing the fitted model results.

- ...:

  Additional arguments (currently unused).

## Value

Prints to the console:

- Model fit statistics (log-likelihood, AIC, BIC, number of
  observations).

- Coefficient tables with standard errors and significance stars.

Invisibly returns `NULL`.

## Details

This method displays the maximum likelihood estimation results for the
Heckman sample selection model with Student's t-distributed errors. It
includes separate coefficient tables for:

- Selection equation (Probit model),

- Outcome equation,

- Error terms (including `sigma`, `rho`, and `df`).

Model fit statistics (log-likelihood, AIC, BIC, and number of
observations) are also provided for model evaluation.

## See also

[`HeckmantS`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmantS.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(MEPS2001)
attach(MEPS2001)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
model <- HeckmantS(selectEq, outcomeEq, data = MEPS2001, df = 12)
summary(model)
} # }
```
