#' Extract model components from formulas
#'
#' @description
#' This function extracts the model frame, model matrix, and model response from formulas and a dataset.
#' It is used in models such as HeckmanCL, HeckmanGe, HeckmanSK, HeckmanBS, and HeckmantS.
#' Optionally, it handles covariate matrices for modeling dispersion and correlation parameters.
#'
#' @param selection A formula representing the selection equation.
#' @param outcome A formula representing the outcome equation.
#' @param data A data frame containing all variables used in the formulas.
#' @param outcomeS Optional matrix or formula of covariates for the dispersion (sigma) parameter.
#' @param outcomeC Optional matrix or formula of covariates for the correlation (rho) parameter.
#' @param drop.levels Logical; if TRUE, unused factor levels are dropped.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{XS}: model matrix for selection equation.
#'   \item \code{YS}: response vector for selection equation.
#'   \item \code{NXS}: number of covariates in the selection model.
#'   \item \code{XO}: model matrix for outcome equation.
#'   \item \code{YO}: response vector for outcome equation.
#'   \item \code{NXO}: number of covariates in the outcome model.
#'   \item \code{Msigma}: matrix for modeling dispersion (NULL if not specified).
#'   \item \code{NE}: number of covariates for the dispersion model.
#'   \item \code{Mrho}: matrix for modeling correlation (NULL if not specified).
#'   \item \code{NV}: number of covariates for the correlation model.
#'   \item \code{YSLevels}: factor levels of the selection response variable.
#' }
#'
#' @importFrom stats model.frame model.matrix model.response terms
#' @export
extract_model_components <- function(selection, outcome, data, outcomeS = NULL, outcomeC = NULL, drop.levels = TRUE) {

  # Build model frame for the selection equation
  mfS <- model.frame(
    formula = selection,
    data = data,
    drop.unused.levels = drop.levels,
    na.action = na.pass
  )
  mtS <- terms(mfS)
  XS <- model.matrix(mtS, mfS)
  YS <- model.response(mfS)
  NXS <- ncol(XS)

  # Ensure selection outcome has exactly two levels (binary response)
  if (length(levels(as.factor(YS))) != 2) {
    stop("The selection equation must have exactly two levels (e.g., 0 and 1).")
  }

  YSLevels <- levels(as.factor(YS))

  # Build model frame for the outcome equation
  mfO <- model.frame(
    formula = outcome,
    data = data,
    drop.unused.levels = drop.levels,
    na.action = na.pass
  )
  mtO <- terms(mfO)
  XO <- model.matrix(mtO, mfO)
  YO <- model.response(mfO)
  NXO <- ncol(XO)

  # Construct dispersion matrix (Msigma) if provided
  if (!is.null(outcomeS)) {
    if (length(outcomeS) == 1) {
      # Only one covariate: intercept-only model for sigma
      NE <- 1
      Msigma <- cbind(rep(1, nrow(XO)))
    } else {
      # Multiple covariates: use model matrix excluding intercept
      Msigma <- outcomeS
      NE <- dim(model.matrix(~outcomeS))[2] - 1
    }
  } else {
    Msigma <- NULL
    NE <- 0
  }

  # Construct correlation matrix (Mrho) if provided
  if (!is.null(outcomeC)) {
    if (length(outcomeC) == 1) {
      # Only one covariate: intercept-only model for rho
      NV <- 1
      Mrho <- cbind(rep(1, nrow(XO)))
    } else {
      # Multiple covariates: use model matrix excluding intercept
      Mrho <- outcomeC
      NV <- dim(model.matrix(~outcomeC))[2] - 1
    }
  } else {
    Mrho <- NULL
    NV <- 0
  }

  # Return all extracted and constructed components
  return(list(
    XS = XS, YS = YS, NXS = NXS,
    XO = XO, YO = YO, NXO = NXO,
    Msigma = Msigma, NE = NE,
    Mrho = Mrho, NV = NV,
    YSLevels = YSLevels
  ))
}
