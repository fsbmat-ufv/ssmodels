#' Summary of Heckman-ts Model
#'
#' @return
#' Print estimates of the parameters of the Heckman-ts model
#' @param object HeckmantS class object.
#' @param ... others functions.
#' @export summary.HeckmantS
#' @export
summary.HeckmantS   <- function(object, ... ){

  fisher_infotS <- object$fisher_infotS
  prop_sigmatS  <- object$prop_sigmatS
  coeffstS      <- object$coefficients
  counts        <- object$counts
  value         <- object$value
  loglik        <- object$loglik
  NObs          <- object$NObs
  nParam        <- object$nParam
  df            <- object$df
  NXS           <- object$NXS
  NXO           <- object$NXO
  N0            <- object$N0
  N1            <- object$N1
  aic           <- object$aic
  bic           <- object$bic

  tb <- miscTools::coefTable(coeffstS,
    prop_sigmatS,
    df = df)

  tb1 <- miscTools::coefTable(coeffstS[1:NXS],
    prop_sigmatS[1:NXS],
    df = df)

  tb2 <- miscTools::coefTable(coeffstS[(NXS + 1):(NXS + NXO)],
    prop_sigmatS[(NXS + 1):(NXS + NXO)],
    df = df)

  tb3 <- miscTools::coefTable(coeffstS[(NXS + NXO + 1):(NXS + NXO + 3)],
    prop_sigmatS[(NXS + NXO + 1):(NXS + NXO + 3)],
    df = df)

  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("      t-Student Heckman Model (Package: ssmodels)             \n")
  cat("--------------------------------------------------------------\n")
  cat("--------------------------------------------------------------\n")
  cat("Maximum Likelihood estimation \n")
  cat("optim function with method BFGS-iterations numbers:", counts,
    "\n")
  cat("Log-Likelihood:", value, "\n")
  cat("AIC:", aic, "BIC:", bic, "\n")
  cat("Number of observations:", NObs, "(", N0, "censored and",
    N1, "observed", ")", "\n")
  cat(nParam, "free parameters", "(", "df=", df, ")", "\n")
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
