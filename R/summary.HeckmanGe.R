#' Summary of Generalized Heckman Model
#'
#' @description
#' Print estimates of the parameters of the Generalized Heckman model using Maximum Likelihood Estimation.
#'
#' @param object An object of class \code{HeckmanGe}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' Prints the summary output including coefficient tables and model fit statistics.
#'
#' @export summary.HeckmanGe
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
