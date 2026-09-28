# Generalized Heckman Model Estimation

Fits a generalized Heckman sample selection model that allows for
heteroskedasticity in the outcome equation and correlation of the error
terms depending on covariates. The estimation is performed via Maximum
Likelihood using the BFGS algorithm.

## Usage

``` r
HeckmanGe(
  selection,
  outcome,
  outcomeS,
  outcomeC,
  data = sys.frame(sys.parent()),
  start = NULL
)
```

## Arguments

- selection:

  A formula specifying the selection equation.

- outcome:

  A formula specifying the outcome equation.

- outcomeS:

  A formula or matrix specifying covariates for the scale (variance)
  model.

- outcomeC:

  A formula or matrix specifying covariates for the correlation model.

- data:

  A data frame containing the variables in the model.

- start:

  An optional numeric vector with starting values for the optimization.

## Value

A list containing:

- `coefficients`: Named vector of estimated model parameters.

- `value`: Negative of the maximum log-likelihood.

- `loglik`: Maximum log-likelihood.

- `counts`: Number of gradient evaluations performed.

- `hessian`: Hessian matrix at the optimum.

- `fisher_infoHG`: Approximate Fisher information matrix.

- `prop_sigmaHG`: Standard errors for the parameter estimates.

- `level`: Levels of the selection variable.

- `nObs`: Number of observations in the dataset.

- `nParam`: Number of estimated parameters.

- `N0`: Number of censored (unobserved) observations.

- `N1`: Number of uncensored (observed) observations.

- `NXS`: Number of covariates in the selection equation.

- `NXO`: Number of covariates in the outcome equation.

- `df`: Degrees of freedom (observations minus parameters).

- `aic`: Akaike Information Criterion.

- `bic`: Bayesian Information Criterion.

- `initial.value`: Starting values used for optimization.

- `NE`: Number of parameters in the scale model.

- `NV`: Number of parameters in the correlation model.

## Details

This function extends the classical Heckman selection model by
incorporating models for the error term's variance (scale) and the
correlation between the selection and outcome equations. The scale model
(`outcomeS`) allows the error variance of the outcome equation to depend
on covariates, while the correlation model (`outcomeC`) allows the error
correlation to vary with covariates.

The optimization is initialized with default or user-supplied starting
values, and the results include robust standard errors derived from the
inverse of the observed Fisher information matrix.

## References

Fernando de Souza Bastos, Wagner Barreto-Souza, Marc G Genton (2022). “A
Generalized Heckman Model With Varying Sample Selection Bias and
Dispersion Parameters.” *Statistica Sinica*.

## Examples

``` r
if (FALSE) { # \dontrun{
data(MEPS2001)
attach(MEPS2001)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
outcomeS <- ~ educ + income
outcomeC <- ~ blhisp + female
HeckmanGe(selectEq, outcomeEq, outcomeS = outcomeS, outcomeC = outcomeC, data = MEPS2001)
} # }
```
