# Extract Model Components for Selection Models

This internal utility function extracts key components—such as model
frames, matrices, and response variables—from formulas and a data set.
It is used by models like `HeckmanCL`, `HeckmanGe`, `HeckmanSK`,
`HeckmanBS`, and `HeckmantS`. Additionally, it can handle covariate
matrices for modeling dispersion (`sigma`) and correlation (`rho`)
structures.

## Usage

``` r
extract_model_components(
  selection,
  outcome,
  data,
  outcomeS = NULL,
  outcomeC = NULL,
  drop.levels = TRUE
)
```

## Arguments

- selection:

  A formula for the selection equation.

- outcome:

  A formula for the outcome equation.

- data:

  A data frame containing all variables.

- outcomeS:

  Optional formula or matrix for the dispersion model (`sigma`).

- outcomeC:

  Optional formula or matrix for the correlation model (`rho`).

- drop.levels:

  Logical. If `TRUE`, drops unused factor levels.

## Value

A list with the following components:

- `XS`:

  Model matrix for the selection equation.

- `YS`:

  Response vector for the selection equation.

- `NXS`:

  Number of covariates in the selection model.

- `XO`:

  Model matrix for the outcome equation.

- `YO`:

  Response vector for the outcome equation.

- `NXO`:

  Number of covariates in the outcome model.

- `Msigma`:

  Matrix for the dispersion model (or `NULL` if not provided).

- `NE`:

  Number of covariates for the dispersion model (0 if not provided).

- `Mrho`:

  Matrix for the correlation model (or `NULL` if not provided).

- `NV`:

  Number of covariates for the correlation model (0 if not provided).

- `YSLevels`:

  Factor levels of the binary selection response.

## Details

If provided, `outcomeS` and `outcomeC` can be formulas or matrices for
modeling dispersion and correlation structures, respectively. The
function ensures that the selection equation response is binary.
