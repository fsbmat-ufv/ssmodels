#' Summary of Classic Heckman Model
#'
#' Print estimates of the parameters of the Classic Heckman model
#' @param object HeckmanCL class object.
#' @param ... others functions.
#' @export summary.HeckmanCL
#' @export
summary.HeckmanCL   <- function(object, ... ){

      fisher_infoHC <- object$fisher_infoHC
      prop_sigmaHC  <- object$prop_sigmaHC
      coeffsHC      <- object$coefficients
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

 tb <- miscTools::coefTable(coeffsHC,
   prop_sigmaHC,
   df = df)

tb1 <- miscTools::coefTable(coeffsHC[1:NXS],
  prop_sigmaHC[1:NXS],
  df = df)

tb2 <- miscTools::coefTable(coeffsHC[(NXS + 1):(NXS + NXO)],
    prop_sigmaHC[(NXS + 1):(NXS + NXO)],
  df = df)

tb3 <- miscTools::coefTable(coeffsHC[(NXS + NXO + 1):(NXS + NXO + 2)],
    prop_sigmaHC[(NXS + NXO + 1):(NXS + NXO + 2)],
  df = df)

cat("\n")
cat("--------------------------------------------------------------\n")
cat("          Classic Heckman Model (Package: ssmodels)           \n")
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
