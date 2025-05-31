#' Summary of Generalized Heckman Model
#'
#' @description
#' Prints a detailed summary of the parameter estimates and model fit
#' statistics for an object of class \code{HeckmanGe}.
#'
#' @details
#' This method displays the maximum likelihood estimation results
#' for the generalized Heckman sample selection model. It includes
#' separate coefficient tables for:
#' \itemize{
#'   \item Selection equation (Probit model),
#'   \item Outcome equation,
#'   \item Dispersion (scale) model parameters,
#'   \item Correlation model parameters.
#' }
#' Model fit statistics (log-likelihood, AIC, BIC, and number of observations)
#' are also reported for interpretation and model assessment.
#'
#' @param object An object of class \code{HeckmanGe}, containing the
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
#' outcomeS <- ~ educ + income
#' outcomeC <- ~ blhisp + female
#' model <- HeckmanGe(selectEq, outcomeEq, outcomeS = outcomeS, outcomeC = outcomeC, data = MEPS2001)
#' summary(model)
#' }
#'
#' @seealso \code{\link{HeckmanGe}}
#'
#' @export
summary.HeckmanGe <- function(object, ...) {
  # Extract components from the HeckmanGe object
  fisher_infoHG <- object$fisher_infoHG       # Fisher information matrix (not printed)
  prop_sigmaHG  <- object$prop_sigmaHG        # Standard errors of the estimates
  coeffsHG      <- object$coefficients        # Estimated coefficients
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
  NE            <- object$NE                  # Number of covariates in scale model
  NV            <- object$NV                  # Number of covariates in correlation model

  # Compute coefficient tables depending on the structure of the model
  if (NE == 1 & NV == 1) {
    tb1 <- miscTools::coefTable(coeffsHG[1:NXS], prop_sigmaHG[1:NXS], df = df)
    tb2 <- miscTools::coefTable(coeffsHG[(NXS + 1):(NXS + NXO)], prop_sigmaHG[(NXS + 1):(NXS + NXO)], df = df)
    tb3 <- miscTools::coefTable(exp(coeffsHG[NXS + NXO + NE]), prop_sigmaHG[NXS + NXO + NE], df = df)
    tb4 <- miscTools::coefTable(tanh(coeffsHG[NXS + NXO + NE + NV]), prop_sigmaHG[NXS + NXO + NE + NV], df = df)

  } else if (NE != 1 & NV == 1) {
    tb1 <- miscTools::coefTable(coeffsHG[1:NXS], prop_sigmaHG[1:NXS], df = df)
    tb2 <- miscTools::coefTable(coeffsHG[(NXS + 1):(NXS + NXO)], prop_sigmaHG[(NXS + 1):(NXS + NXO)], df = df)
    tb3 <- miscTools::coefTable(coeffsHG[(NXS + NXO + 1):(NXS + NXO + NE)], prop_sigmaHG[(NXS + NXO + 1):(NXS + NXO + NE)], df = df)
    tb4 <- miscTools::coefTable(tanh(coeffsHG[NXS + NXO + NE + NV]), prop_sigmaHG[NXS + NXO + NE + NV], df = df)

  } else if (NE == 1 & NV != 1) {
    tb1 <- miscTools::coefTable(coeffsHG[1:NXS], prop_sigmaHG[1:NXS], df = df)
    tb2 <- miscTools::coefTable(coeffsHG[(NXS + 1):(NXS + NXO)], prop_sigmaHG[(NXS + 1):(NXS + NXO)], df = df)
    tb3 <- miscTools::coefTable(exp(coeffsHG[NXS + NXO + NE]), prop_sigmaHG[NXS + NXO + NE], df = df)
    tb4 <- miscTools::coefTable(coeffsHG[(NXS + NXO + NE + 1):(NXS + NXO + NE + NV)], prop_sigmaHG[(NXS + NXO + NE + 1):(NXS + NXO + NE + NV)], df = df)

  } else {
    tb1 <- miscTools::coefTable(coeffsHG[1:NXS], prop_sigmaHG[1:NXS], df = df)
    tb2 <- miscTools::coefTable(coeffsHG[(NXS + 1):(NXS + NXO)], prop_sigmaHG[(NXS + 1):(NXS + NXO)], df = df)
    tb3 <- miscTools::coefTable(coeffsHG[(NXS + NXO + 1):(NXS + NXO + NE)], prop_sigmaHG[(NXS + NXO + 1):(NXS + NXO + NE)], df = df)
    tb4 <- miscTools::coefTable(coeffsHG[(NXS + NXO + NE + 1):(NXS + NXO + NE + NV)], prop_sigmaHG[(NXS + NXO + NE + 1):(NXS + NXO + NE + NV)], df = df)
  }

  # Print summary output
  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("       Generalized Heckman Model (Package: ssmodels)          \n")
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

  cat("Dispersion terms:\n")
  printCoefmat(tb3, signif.stars = TRUE, signif.legend = FALSE, digits = 4)
  cat("--------------------------------------------------------------\n")

  cat("Correlation terms:\n")
  printCoefmat(tb4, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
  cat("--------------------------------------------------------------\n")
}
