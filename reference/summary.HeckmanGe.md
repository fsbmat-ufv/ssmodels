# Summary of Generalized Heckman Model

Prints a detailed summary of the parameter estimates and model fit
statistics for an object of class `HeckmanGe`.

## Usage

``` r
# S3 method for class 'HeckmanGe'
summary(object, ...)
```

## Arguments

- object:

  An object of class `HeckmanGe`, containing the fitted model results.

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
generalized Heckman sample selection model. It includes separate
coefficient tables for:

- Selection equation (Probit model),

- Outcome equation,

- Dispersion (scale) model parameters,

- Correlation model parameters.

Model fit statistics (log-likelihood, AIC, BIC, and number of
observations) are also reported for interpretation and model assessment.

## See also

[`HeckmanGe`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmanGe.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(MEPS2001)
attach(MEPS2001)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
outcomeS <- ~ educ + income
outcomeC <- ~ blhisp + female
model <- HeckmanGe(selectEq, outcomeEq, outcomeS = outcomeS, outcomeC = outcomeC, data = MEPS2001)
summary(model)
} # }
```
