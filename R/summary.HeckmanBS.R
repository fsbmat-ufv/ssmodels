#' Summary of Birnbaum-Saunders Heckman Model
#'
#' @return
#' Print estimates of the parameters of the Heckman-BS model
#' @param object HeckmanBS class object.
#' @param ... others functions.
#' @export summary.HeckmanBS
#' @export
summary.HeckmanBS   <- function (object, ...)
{
  fisher_infoBS <- object$fisher_infoBS
  prop_sigmaBS <- object$prop_sigmaBS
  coeffsBS <- object$coefficients
  counts <- object$counts
  value <- object$value
  loglik <- object$loglik
  NObs <- object$NObs
  nParam <- object$nParam
  df <- object$df
  NXS <- object$NXS
  NXO <- object$NXO
  N0 <- object$N0
  N1 <- object$N1
  aic <- object$aic
  bic <- object$bic
  tb <- miscTools::coefTable(coeffsBS, prop_sigmaBS, df = df)
  tb1 <- miscTools::coefTable(coeffsBS[1:NXS], prop_sigmaBS[1:NXS],
                              df = df)
  tb2 <- miscTools::coefTable(coeffsBS[(NXS + 1):(NXS + NXO)],
                              prop_sigmaBS[(NXS + 1):(NXS + NXO)], df = df)
  tb3 <- miscTools::coefTable(coeffsBS[(NXS + NXO + 1):(NXS +
                                                          NXO + 2)], prop_sigmaBS[(NXS + NXO + 1):(NXS + NXO +
                                                                                                     2)], df = df)
  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("   Birnbaum-Saunders Heckman Model (Package: ssmodels)        \n")
  cat("--------------------------------------------------------------\n")
  cat("--------------------------------------------------------------\n")
  cat("Maximum Likelihood estimation \n")
  cat("optim function with method BFGS-iterations numbers:",
      counts, "\n")
  cat("Log-Likelihood:", value, "\n")
  cat("AIC:", aic, "BIC:", bic, "\n")
  cat("Number of observations:", NObs, "(", N0, "censored and",
      N1, "observed", ")", "\n")
  cat(nParam, "free parameters", "(", "df=", df, ")", "\n")
  cat("--------------------------------------------------------------\n")
  cat("Probit selection equation:\n")
  printCoefmat(tb1, signif.stars = TRUE, signif.legend = FALSE,
               digits = 4)
  cat("--------------------------------------------------------------\n")
  cat("Outcome equation:\n")
  printCoefmat(tb2, signif.stars = TRUE, signif.legend = FALSE,
               digits = 4)
  cat("--------------------------------------------------------------\n")
  cat("Error terms:\n")
  printCoefmat(tb3, signif.stars = TRUE, signif.legend = TRUE,
               digits = 4)
  cat("--------------------------------------------------------------\n")
}
