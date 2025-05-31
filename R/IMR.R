#' Inverse Mills Ratio (IMR) Calculation
#'
#' @description
#' Computes the column vector of the Inverse Mills Ratio (IMR) from
#' a Probit selection equation.
#'
#' @details
#' This function fits a Probit model to the provided selection equation
#' and returns the Inverse Mills Ratio (IMR) for each observation.
#' The IMR is useful for correcting sample selection bias in regression
#' models, following the classical Heckman approach.
#'
#' @param selection A formula specifying the selection equation.
#' @param data A data frame containing the variables in the model.
#'
#' @return
#' A numeric matrix with one column containing the IMR values
#' for each observation.
#'
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' IMR(selectEq, data = MEPS2001)
#'
#' @export
IMR <- function(selection, data = sys.frame(sys.parent())) {

  # Step 1: Build model frame from selection equation
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("selection", "data", "subset"), names(mf), 0)
  mfS <- mf[c(1, m)]
  mfS$drop.unused.levels <- TRUE
  mfS$na.action <- na.pass
  mfS[[1]] <- as.name("model.frame")
  names(mfS)[2] <- "formula"
  mfS <- eval(mfS, parent.frame())

  # Step 2: Extract design matrix and response vector
  mtS <- terms(mfS)
  XS <- model.matrix(mtS, mfS)
  YS <- model.response(mfS)

  # Step 3: Fit Probit model to selection equation
  fit1 <- glm(YS ~ XS - 1, family = binomial(link = "probit"))

  # Step 4: Compute Inverse Mills Ratio (IMR)
  imr <- dnorm(fit1$linear.predictors) / pnorm(fit1$linear.predictors)

  # Step 5: Return IMR as column matrix
  return(cbind(imr))
}
