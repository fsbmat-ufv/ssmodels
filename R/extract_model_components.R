#' Extract model components from formulas
#'
#' @description
#' This function extracts the model.frame, model.matrix and model.response from a formula and data.
#' For models like HeckmanCL, HeckmanGe, HeckmanSK, HeckmanBS and HeckmantS, it also handles
#' optional matrices for dispersion (outcomeS) and correlation (outcomeC).
#'
#' @param selection A formula for the selection equation.
#' @param outcome A formula for the outcome equation.
#' @param data A data frame containing the variables in the formulas.
#' @param outcomeS Optional: covariates for the dispersion model.
#' @param outcomeC Optional: covariates for the correlation model.
#' @param drop.levels Logical. If TRUE, drops unused factor levels.
#'
#' @return A list with components for selection and outcome models, including optional matrices for dispersion and correlation.
#'
#' @importFrom stats model.frame model.matrix model.response terms
#' @export
extract_model_components <- function(selection, outcome, data, outcomeS = NULL, outcomeC = NULL, drop.levels = TRUE) {

  # Selection equation
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

  if (length(levels(as.factor(YS))) != 2) {
    stop("The selection equation must have exactly two levels (e.g., 0 and 1).")
  }

  YSLevels <- levels(as.factor(YS))

  # Outcome equation
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

  # Dispersion Matrix
  if (!is.null(outcomeS)) {
    if (length(outcomeS) == 1) {
      NE <- 1
      Msigma <- cbind(rep(1, nrow(XO)))
    } else {
      Msigma <- outcomeS
      NE <- dim(model.matrix(~outcomeS))[2] - 1
    }
  } else {
    Msigma <- NULL
    NE <- 0
  }

  # Correlation Matrix
  if (!is.null(outcomeC)) {
    if (length(outcomeC) == 1) {
      NV <- 1
      Mrho <- cbind(rep(1, nrow(XO)))
    } else {
      Mrho <- outcomeC
      NV <- dim(model.matrix(~outcomeC))[2] - 1
    }
  } else {
    Mrho <- NULL
    NV <- 0
  }

  return(list(
    XS = XS, YS = YS, NXS = NXS,
    XO = XO, YO = YO, NXO = NXO,
    Msigma = Msigma, NE = NE,
    Mrho = Mrho, NV = NV,
    YSLevels = YSLevels
  ))
}
