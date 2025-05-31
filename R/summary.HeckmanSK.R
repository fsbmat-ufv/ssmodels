#' Summary of Skew-Normal Heckman Model
#'
#' @description
#' Prints a detailed summary of the parameter estimates and model fit
#' statistics for an object of class \code{HeckmanSK}.
#'
#' @details
#' This method displays the maximum likelihood estimation results
#' for the Heckman sample selection model with Skew-Normal errors.
#' It includes separate coefficient tables for:
#' \itemize{
#'   \item Selection equation (Probit model),
#'   \item Outcome equation,
#'   \item Error terms (\code{sigma}, \code{rho}, and \code{lambda}).
#' }
#' Additionally, it reports model fit statistics such as the log-likelihood,
#' AIC, BIC, and the number of observations.
#'
#' @param object An object of class \code{HeckmanSK}, containing the
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
#' data(Mroz87)
#' attach(Mroz87)
#' selectEq <- lfp ~ huswage + kids5 + mtr + fatheduc + educ + city
#' outcomeEq <- log(wage) ~ educ + city
#' model <- HeckmanSK(selectEq, outcomeEq, data = Mroz87, lambda = -1.5)
#' summary(model)
#' }
#'
#' @seealso \code{\link{HeckmanSK}}
#'
#' @export
summary.HeckmanSK <- function(object, ...) {
  # Extract components from the HeckmanSK object
  fisher_infoSK <- object$fisher_infoSK       # Fisher information matrix (not printed)
  prop_sigmaSK  <- object$prop_sigmaSK        # Standard errors of the estimates
  coeffsSK      <- object$coefficients        # Estimated coefficients
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
  tb1 <- miscTools::coefTable(coeffsSK[1:NXS], prop_sigmaSK[1:NXS], df = df)                        # Selection equation
  tb2 <- miscTools::coefTable(coeffsSK[(NXS + 1):(NXS + NXO)], prop_sigmaSK[(NXS + 1):(NXS + NXO)], df = df)  # Outcome equation
  tb3 <- miscTools::coefTable(coeffsSK[(NXS + NXO + 1):(NXS + NXO + 3)], prop_sigmaSK[(NXS + NXO + 1):(NXS + NXO + 3)], df = df)  # Error components
  tb  <- miscTools::coefTable(coeffsSK, prop_sigmaSK, df = df)  # Full model (unused in print but preserved)

  # Print summary output
  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("      Skew Normal Heckman Model (Package: ssmodels)           \n")
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
