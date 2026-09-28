# Summary of Birnbaum-Saunders Heckman Model

Prints a detailed summary of the parameter estimates and model fit
statistics for an object of class `HeckmanBS`.

## Usage

``` r
# S3 method for class 'HeckmanBS'
summary(object, ...)
```

## Arguments

- object:

  An object of class `HeckmanBS`, containing the fitted model results.

- ...:

  Additional arguments (currently unused).

## Value

Prints to the console:

- Model fit statistics (log-likelihood, AIC, BIC, number of
  observations).

- Coefficient tables with standard errors and significance stars.

Invisibly returns `NULL`.

## Details

This method provides a summary of the maximum likelihood estimation
results for the Heckman sample selection model with Birnbaum-Saunders
errors. It includes separate coefficient tables for:

- Selection equation (Probit model),

- Outcome equation,

- Error terms (`sigma` and `rho`).

Model fit criteria such as the log-likelihood, AIC, and BIC are also
reported.

## See also

[`HeckmanBS`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmanBS.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(MEPS2001)
attach(MEPS2001)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- ambexp ~ age + female + educ + blhisp + totchr + ins
model <- HeckmanBS(selectEq, outcomeEq, data = MEPS2001)
summary(model)
} # }
```
