# Post-process Parameter Vector for Generalized Heckman Models

Internal helper function that assigns meaningful names to a vector of
estimated parameters and applies transformations to `sigma` and `rho`
parameters to obtain interpretable estimates.

## Usage

``` r
postprocess_theta(theta_par, NXS, NXO, NE, NV, XS, XO, outcomeS, outcomeC)
```

## Arguments

- theta_par:

  A numeric vector containing the estimated parameters.

- NXS:

  Integer. Number of covariates in the selection equation.

- NXO:

  Integer. Number of covariates in the outcome equation.

- NE:

  Integer. Number of covariates (or 1 for intercept-only) in the scale
  model.

- NV:

  Integer. Number of covariates (or 1 for intercept-only) in the
  correlation model.

- XS:

  Design matrix for the selection equation (used for naming).

- XO:

  Design matrix for the outcome equation (used for naming).

- outcomeS:

  Design matrix or variable for the scale (variance) model.

- outcomeC:

  Design matrix or variable for the correlation model.

## Value

A named numeric vector with:

- Transformed `sigma` and `rho` values (if needed),

- Meaningful names assigned to all parameters.

## Details

The parameter vector `theta_par` is structured as follows:

- The first `NXS` elements are coefficients for the selection equation.

- The next `NXO` elements are coefficients for the outcome equation.

- The next `NE` elements are coefficients (or log-sigma if `NE == 1`)
  for the scale model.

- The next `NV` elements are coefficients (or atanh(rho) if `NV == 1`)
  for the correlation model.

For interpretation, the function applies:

- [`exp()`](https://rdrr.io/r/base/Log.html) transformation for `sigma`
  if `NE == 1`.

- [`tanh()`](https://rdrr.io/r/base/Hyperbolic.html) transformation for
  `rho` if `NV == 1`.
