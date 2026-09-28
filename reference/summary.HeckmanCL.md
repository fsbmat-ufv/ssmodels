# Summary of Classic Heckman Model

Prints a detailed summary of the parameter estimates and model fit
statistics for an object of class `HeckmanCL`.

## Usage

``` r
# S3 method for class 'HeckmanCL'
summary(object, ...)
```

## Arguments

- object:

  An object of class `HeckmanCL`, containing the fitted model results.

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
classical Heckman sample selection model. It includes separate
coefficient tables for:

- Selection equation (Probit model),

- Outcome equation,

- Error terms (`sigma` and `rho`).

Additionally, it reports the model fit statistics (log-likelihood, AIC,
BIC, and number of observations).

## See also

[`HeckmanCL`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmanCL.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(MEPS2001)
attach(MEPS2001)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
model <- HeckmanCL(selectEq, outcomeEq, data = MEPS2001)
summary(model)
} # }
```
