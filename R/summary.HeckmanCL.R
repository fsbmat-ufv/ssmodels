#' Summary of Classic Heckman Model
#'
#' @description
#' Prints a detailed summary of the parameter estimates and model fit
#' statistics for an object of class \code{HeckmanCL}.
#'
#' @details
#' This method displays the maximum likelihood estimation results
#' for the classical Heckman sample selection model. It includes
#' separate coefficient tables for:
#' \itemize{
#'   \item Selection equation (Probit model),
#'   \item Outcome equation,
#'   \item Error terms (\code{sigma} and \code{rho}).
#' }
#' Additionally, it reports the model fit statistics (log-likelihood,
#' AIC, BIC, and number of observations).
#'
#' @param object An object of class \code{HeckmanCL}, containing the
#' fitted model results.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' Prints to the console:
#' \itemize{
#'   \item Model fit statistics (log-likelihood, AIC, BIC, number of observations).
#'   \item Coefficient tables with standard errors and significance stars.
#' }
#' Invisibly returns \code{NULL}.
#'
#' @examples
#' \dontrun{
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' model <- HeckmanCL(selectEq, outcomeEq, data = MEPS2001)
#' summary(model)
#' }
#'
#' @seealso \code{\link{HeckmanCL}}
#'
#' @export
summary.HeckmanCL <- function(object, ...) {

  # Extract components from the HeckmanCL object
  fisher_infoHC <- object$fisher_infoHC       # Fisher information matrix (not used here)
  prop_sigmaHC  <- object$prop_sigmaHC        # Standard errors of the estimates
  coeffsHC      <- object$coefficients        # Estimated coefficients
  counts        <- object$counts              # Number of iterations in optimization
  value         <- object$value               # Log-likelihood at convergence
  loglik        <- object$loglik              # Log-likelihood function
  NObs          <- object$NObs                # Total number of observations
  nParam        <- object$nParam              # Number of free parameters
  df            <- object$df                  # Degrees of freedom
  NXS           <- object$NXS                 # Number of selection equation covariates
  NXO           <- object$NXO                 # Number of outcome equation covariates
  N0            <- object$N0                  # Number of censored observations
  N1            <- object$N1                  # Number of uncensored observations
  aic           <- object$aic                 # Akaike Information Criterion
  bic           <- object$bic                 # Bayesian Information Criterion

  # Build coefficient tables for each block of the model
  tb  <- miscTools::coefTable(coeffsHC, prop_sigmaHC, df = df)

  tb1 <- miscTools::coefTable(
    coeffsHC[1:NXS],
    prop_sigmaHC[1:NXS],
    df = df
  )

  tb2 <- miscTools::coefTable(
    coeffsHC[(NXS + 1):(NXS + NXO)],
    prop_sigmaHC[(NXS + 1):(NXS + NXO)],
    df = df
  )

  tb3 <- miscTools::coefTable(
    coeffsHC[(NXS + NXO + 1):(NXS + NXO + 2)],
    prop_sigmaHC[(NXS + NXO + 1):(NXS + NXO + 2)],
    df = df
  )

  # Print model summary
  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("          Classic Heckman Model (Package: ssmodels)           \n")
  cat("--------------------------------------------------------------\n")
  cat("--------------------------------------------------------------\n")
  cat("Maximum Likelihood estimation \n")
  cat("optim function with method BFGS - iterations number:", counts, "\n")
  cat("Log-Likelihood:", value, "\n")
  cat("AIC:", aic, "BIC:", bic, "\n")
  cat("Number of observations:", NObs, "(", N0, "censored and", N1, "observed", ")", "\n")
  cat(nParam, "free parameters", "( df =", df, ")", "\n")
  cat("--------------------------------------------------------------\n")

  # Print selection equation coefficients
  cat("Probit selection equation:\n")
  printCoefmat(tb1, signif.stars = TRUE, signif.legend = FALSE, digits = 4)
  cat("--------------------------------------------------------------\n")

  # Print outcome equation coefficients
  cat("Outcome equation:\n")
  printCoefmat(tb2, signif.stars = TRUE, signif.legend = FALSE, digits = 4)
  cat("--------------------------------------------------------------\n")

  # Print error term estimates (sigma and rho)
  cat("Error terms:\n")
  printCoefmat(tb3, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
  cat("--------------------------------------------------------------\n")
}
