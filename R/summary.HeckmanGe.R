#' Summary of Generalized Heckman Model
#'
#' Print eEstimates of the parameters of the Generalized Heckman model
#' @param object HeckmanGe class object.
#' @param ... others functions.
#' @export summary.HeckmanGe
#' @export
summary.HeckmanGe   <- function(object, ... ){
  fisher_infoHG <- object$fisher_infoHG
  prop_sigmaHG  <- object$prop_sigmaHG
  coeffsHG      <- object$coefficients
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
  NE <- object$NE
  NV <- object$NV

if(NE == 1 & NV == 1){
   tb <- miscTools::coefTable(coeffsHG,
     prop_sigmaHG,
     df = df)

  tb1 <- miscTools::coefTable(coeffsHG[1:NXS],
    prop_sigmaHG[1:NXS],
    df = df)

  tb2 <- miscTools::coefTable(coeffsHG[(NXS + 1):(NXS + NXO)],
    prop_sigmaHG[(NXS + 1):(NXS + NXO)],
    df = df)

  tb3 <- miscTools::coefTable(exp(coeffsHG[(NXS + NXO + NE)]),
    prop_sigmaHG[(NXS + NXO + NE)],
    df = df)

  tb4 <- miscTools::coefTable(tanh(coeffsHG[(NXS + NXO + NE + NV)]), prop_sigmaHG[(NXS +
          NXO + NE + NV)], df = df)

} else{
      if(NE != 1 & NV == 1) {
          tb <- miscTools::coefTable(coeffsHG,
            prop_sigmaHG,
            df = df)

          tb1 <- miscTools::coefTable(coeffsHG[1:NXS],
            prop_sigmaHG[1:NXS],
            df = df)

          tb2 <- miscTools::coefTable(coeffsHG[(NXS + 1):(NXS + NXO)],
            prop_sigmaHG[(NXS + 1):(NXS + NXO)],
            df = df)

          tb3 <- miscTools::coefTable(coeffsHG[(NXS + NXO + 1):(NXS + NXO + NE)],
            prop_sigmaHG[(NXS + NXO + 1):(NXS + NXO + NE)],
            df = df)

          tb4 <- miscTools::coefTable(tanh(coeffsHG[(NXS + NXO + NE + NV)]),
            prop_sigmaHG[(NXS + NXO + NE + NV)],
            df = df)

      } else{
          if(NE == 1 & NV != 1) {
              tb <- miscTools::coefTable(coeffsHG,
                prop_sigmaHG,
                df = df)

              tb1 <- miscTools::coefTable(coeffsHG[1:NXS], prop_sigmaHG[1:NXS],
                df = df)

              tb2 <- miscTools::coefTable(coeffsHG[(NXS + 1):(NXS + NXO)],
                prop_sigmaHG[(NXS + 1):(NXS + NXO)],
                df = df)

              tb3 <- miscTools::coefTable(exp(coeffsHG[(NXS + NXO + NE)]),
                prop_sigmaHG[(NXS + NXO + NE)],
                df = df)

              tb4 <- miscTools::coefTable(coeffsHG[(NXS + NXO + NE + 1):(NXS + NXO + NE + NV)],
                prop_sigmaHG[(NXS + NXO + NE + 1):(NXS + NXO + NE + NV)],
                df = df)

          } else{
              tb <- miscTools::coefTable(coeffsHG,
                prop_sigmaHG,
                df = df)

              tb1 <- miscTools::coefTable(coeffsHG[1:NXS],
                prop_sigmaHG[1:NXS],
                df = df)

              tb2 <- miscTools::coefTable(coeffsHG[(NXS + 1):(NXS + NXO)],
                prop_sigmaHG[(NXS + 1):(NXS + NXO)],
                df = df)

              tb3 <- miscTools::coefTable(coeffsHG[(NXS + NXO + 1):(NXS + NXO + NE)],
                prop_sigmaHG[(NXS + NXO + 1):(NXS + NXO + NE)],
                df = df)

              tb4 <- miscTools::coefTable(coeffsHG[(NXS + NXO + NE + 1):(NXS + NXO + NE + NV)],
                prop_sigmaHG[(NXS + NXO + NE + 1):(NXS + NXO + NE + NV)],
                df = df)
          }
      }
  }


  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("       Generalized Heckman Model (Package: ssmodels)          \n")
  cat("--------------------------------------------------------------\n")
  cat("--------------------------------------------------------------\n")
  cat("Maximum Likelihood estimation \n")
  cat("optim function with method BFGS-iterations numbers:", counts, "\n")
  cat("Log-Likelihood:", value, "\n")
  cat("AIC:", aic, "BIC:", bic, "\n")
  cat("Number of observations:", NObs, "(", N0, "censored and", N1, "observed", ")", "\n")
  cat(nParam, "free parameters", "(", "df=", df, ")", "\n")
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
