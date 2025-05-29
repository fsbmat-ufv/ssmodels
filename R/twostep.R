#' Heckman's two-step method
#'
#' @description
#' Estimate model parameters via two-step method.
#'
#' @param selection A formula for the selection equation.
#' @param outcome A formula for the outcome equation.
#' @param data A data frame containing the variables used in the model.
#'
#' @return
#' Returns a numeric vector with the parameter estimates of the Classical Heckman model using the two-step method.
#' For more information, see \insertCite{heckman1979sample;textual}{ssmodels}
#'
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' twostep(selectEq, outcomeEq)
#'
#' @importFrom Rdpack reprompt
#' @references {
#' \insertAllCited{}
#' }
#' @export
twostep <- function(selection, outcome, data = sys.frame(sys.parent())) {
  # Extract selection components
  sel <- extract_model_components(selection, data)
  YS <- sel$y                     # Selection response
  XS <- sel$X                     # Selection design matrix
  NXS <- ncol(XS)                 # Number of selection covariates

  # Extract outcome components
  out <- extract_model_components(outcome, data)
  YO <- out$y                     # Outcome response
  XO <- out$X                     # Outcome design matrix
  NXO <- ncol(XO)                 # Number of outcome covariates

  # Step 1: First-stage estimation (Probit model)
  fit1 <- glm(YS ~ XS - 1, family = binomial(link = "probit"))

  # Compute Inverse Mills Ratio
  IMR <- dnorm(fit1$linear.predictors) / pnorm(fit1$linear.predictors)

  # Step 2: Second-stage regression (OLS with IMR)
  xMat <- cbind(XO, IMR)
  fit2 <- lm(YO ~ XO + IMR - 1, subset = YS == 1)
  xMat <- data.frame(xMat)

  # Compute delta for variance correction
  delta <- xMat$IMR * (xMat$IMR + fit1$fitted.values)

  # Estimate residual variance
  Var <- (sum((fit2$residuals)^2) +
            (coef(fit2)[length(coef(fit2))]^2 * sum(delta[YS == 1]))) / sum(YS)

  # Estimate correlation between error terms
  cor <- coef(fit2)[length(coef(fit2))] / sqrt(Var)

  # Return vector of parameters
  beta <- c(coef(fit1), coef(fit2)[-length(coef(fit2))], phi = sqrt(Var), cor)
  return(beta)
}
