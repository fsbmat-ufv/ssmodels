#' Summary of Heckman-t Model
#'
#' @description
#' Prints a detailed summary of the parameter estimates and model fit
#' statistics for an object of class \code{HeckmantS}.
#'
#' @details
#' This method displays the maximum likelihood estimation results
#' for the Heckman sample selection model with Student's t-distributed errors.
#' It includes separate coefficient tables for:
#' \itemize{
#'   \item Selection equation (Probit model),
#'   \item Outcome equation,
#'   \item Error terms (including \code{sigma}, \code{rho}, and \code{df}).
#' }
#' Model fit statistics (log-likelihood, AIC, BIC, and number of observations)
#' are also provided for model evaluation.
#'
#' @param object An object of class \code{HeckmantS}, containing the
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
#' model <- HeckmantS(selectEq, outcomeEq, data = MEPS2001, df = 12)
#' summary(model)
#' }
#'
#' @seealso \code{\link{HeckmantS}}
#'
#' @export
summary.HeckmantS <- function(object, ...) {
  # Extract components from the HeckmantS object
  fisher_infotS <- object$fisher_infotS       # Fisher information matrix (not printed)
  prop_sigmatS  <- object$prop_sigmatS        # Standard errors of the estimates
  coeffstS      <- object$coefficients        # Estimated coefficients
  counts        <- object$counts              # Number of iterations in optimization
  value         <- object$value               # Log-likelihood at convergence
  loglik        <- object$loglik              # Log-likelihood function
  NObs          <- object$NObs                # Total number of observations
  nParam        <- object$nParam              # Number of free parameters
  df            <- object$df                  # Degrees of freedom
  NXS           <- object$NXS                 # Number of covariates in selection equation
  NXO           <- object$NXO                 # Number of covariates in outcome equation
  N0            <- object$N0                  # Number of censored observations
  N1            <- object$N1                  # Number of observed (uncensored) observations
  aic           <- object$aic                 # Akaike Information Criterion
  bic           <- object$bic                 # Bayesian Information Criterion

  # Generate coefficient tables
  tb1 <- miscTools::coefTable(coeffstS[1:NXS], prop_sigmatS[1:NXS], df = df)                         # Selection equation
  tb2 <- miscTools::coefTable(coeffstS[(NXS + 1):(NXS + NXO)], prop_sigmatS[(NXS + 1):(NXS + NXO)], df = df)  # Outcome equation
  tb3 <- miscTools::coefTable(coeffstS[(NXS + NXO + 1):(NXS + NXO + 3)], prop_sigmatS[(NXS + NXO + 1):(NXS + NXO + 3)], df = df)  # Error and df terms
  tb  <- miscTools::coefTable(coeffstS, prop_sigmatS, df = df)  # Full model (not printed but kept)

  # Print summary output
  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("      t-Student Heckman Model (Package: ssmodels)             \n")
  cat("--------------------------------------------------------------\n")
  cat("--------------------------------------------------------------\n")
  cat("Maximum Likelihood estimation \n")
  cat("optim function with method BFGS - iterations number:", counts, "\n")
  cat("Log-Likelihood:", value, "\n")
  cat("AIC:", aic, "BIC:", bic, "\n")
  cat("Number of observations:", NObs, "(", N0, "censored and", N1, "observed )\n")
  cat(nParam, "free parameters ( df =", df, ")\n")
  cat("--------------------------------------------------------------\n")

  cat("Probit selection equation:\n")
  printCoefmat(tb1, signif.stars = TRUE, signif.legend = FALSE, digits = 4)
  cat("--------------------------------------------------------------\n")

  cat("Outcome equation:\n")
  printCoefmat(tb2, signif.stars = TRUE, signif.legend = FALSE, digits = 4)
  cat("--------------------------------------------------------------\n")

  cat("Error terms:\n")
  printCoefmat(tb3, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
  cat("--------------------------------------------------------------\n")
}
