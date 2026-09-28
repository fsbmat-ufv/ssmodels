# Summary of Skew-Normal Heckman Model

Prints a detailed summary of the parameter estimates and model fit
statistics for an object of class `HeckmanSK`.

## Usage

``` r
# S3 method for class 'HeckmanSK'
summary(object, ...)
```

## Arguments

- object:

  An object of class `HeckmanSK`, containing the fitted model results.

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
Heckman sample selection model with Skew-Normal errors. It includes
separate coefficient tables for:

- Selection equation (Probit model),

- Outcome equation,

- Error terms (`sigma`, `rho`, and `lambda`).

Additionally, it reports model fit statistics such as the
log-likelihood, AIC, BIC, and the number of observations.

## See also

[`HeckmanSK`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmanSK.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(Mroz87)
attach(Mroz87)
selectEq <- lfp ~ huswage + kids5 + mtr + fatheduc + educ + city
outcomeEq <- log(wage) ~ educ + city
model <- HeckmanSK(selectEq, outcomeEq, data = Mroz87, lambda = -1.5)
summary(model)
} # }
```
