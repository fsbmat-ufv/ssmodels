#' Summary of Birnbaum-Saunders Heckman Model
#'
#' @description
#' Prints a detailed summary of the parameter estimates and model fit
#' statistics for an object of class \code{HeckmanBS}.
#'
#' @details
#' This method provides a summary of the maximum likelihood estimation results
#' for the Heckman sample selection model with Birnbaum-Saunders errors.
#' It includes separate coefficient tables for:
#' \itemize{
#'   \item Selection equation (Probit model),
#'   \item Outcome equation,
#'   \item Error terms (\code{sigma} and \code{rho}).
#' }
#' Model fit criteria such as the log-likelihood, AIC, and BIC are also reported.
#'
#' @param object An object of class \code{HeckmanBS}, containing
#' the fitted model results.
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
#' outcomeEq <- ambexp ~ age + female + educ + blhisp + totchr + ins
#' model <- HeckmanBS(selectEq, outcomeEq, data = MEPS2001)
#' summary(model)
#' }
#'
#' @seealso \code{\link{HeckmanBS}}
#'
#' @export
summary.HeckmanBS <- function(object, ...) {

  # Extract components from the HeckmanBS object
  fisher_infoBS <- object$fisher_infoBS       # Fisher information matrix (not used here)
  prop_sigmaBS  <- object$prop_sigmaBS        # Standard errors of the estimates
  coeffsBS      <- object$coefficients        # Estimated coefficients
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
  tb  <- miscTools::coefTable(coeffsBS, prop_sigmaBS, df = df)

  tb1 <- miscTools::coefTable(
    coeffsBS[1:NXS],
    prop_sigmaBS[1:NXS],
    df = df
  )

  tb2 <- miscTools::coefTable(
    coeffsBS[(NXS + 1):(NXS + NXO)],
    prop_sigmaBS[(NXS + 1):(NXS + NXO)],
    df = df
  )

  tb3 <- miscTools::coefTable(
    coeffsBS[(NXS + NXO + 1):(NXS + NXO + 2)],
    prop_sigmaBS[(NXS + NXO + 1):(NXS + NXO + 2)],
    df = df
  )

  # Print model summary
  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("   Birnbaum-Saunders Heckman Model (Package: ssmodels)        \n")
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
