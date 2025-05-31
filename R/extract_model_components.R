#' Extract Model Components for Selection Models
#'
#' @description
#' This internal utility function extracts key components—such as model frames, matrices,
#' and response variables—from formulas and a data set. It is used by models like
#' \code{HeckmanCL}, \code{HeckmanGe}, \code{HeckmanSK}, \code{HeckmanBS}, and \code{HeckmantS}.
#' Additionally, it can handle covariate matrices for modeling dispersion (\code{sigma}) and
#' correlation (\code{rho}) structures.
#'
#' @details
#' If provided, \code{outcomeS} and \code{outcomeC} can be formulas or matrices for modeling
#' dispersion and correlation structures, respectively. The function ensures that the
#' selection equation response is binary.
#'
#' @param selection A formula for the selection equation.
#' @param outcome A formula for the outcome equation.
#' @param data A data frame containing all variables.
#' @param outcomeS Optional formula or matrix for the dispersion model (\code{sigma}).
#' @param outcomeC Optional formula or matrix for the correlation model (\code{rho}).
#' @param drop.levels Logical. If \code{TRUE}, drops unused factor levels.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{XS}}{Model matrix for the selection equation.}
#'   \item{\code{YS}}{Response vector for the selection equation.}
#'   \item{\code{NXS}}{Number of covariates in the selection model.}
#'   \item{\code{XO}}{Model matrix for the outcome equation.}
#'   \item{\code{YO}}{Response vector for the outcome equation.}
#'   \item{\code{NXO}}{Number of covariates in the outcome model.}
#'   \item{\code{Msigma}}{Matrix for the dispersion model (or \code{NULL} if not provided).}
#'   \item{\code{NE}}{Number of covariates for the dispersion model (0 if not provided).}
#'   \item{\code{Mrho}}{Matrix for the correlation model (or \code{NULL} if not provided).}
#'   \item{\code{NV}}{Number of covariates for the correlation model (0 if not provided).}
#'   \item{\code{YSLevels}}{Factor levels of the binary selection response.}
#' }
#'
#' @importFrom stats model.frame model.matrix model.response terms
#' @keywords internal
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
