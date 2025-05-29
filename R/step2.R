#' Heckman's two-step method
#'
#' @description
#' Estimate model parameters via two-step method.
#'
#' @return
#' Returns a numerical vector with the parameter estimates of the Classical
#' Heckman model via a two-step method. For more information see
#' \insertCite{heckman1979sample;textual}{ssmodels}
#'
#' @param YS Selection vector.
#' @param XS Selection matrix.
#' @param YO Outcome vector.
#' @param XO Covariate matrix for the outcome equation.
#'
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' YS <- dambexp
#' XS <- cbind(age, female, educ, blhisp, totchr, ins)
#' YO <- lnambx
#' XO <- cbind(age, female, educ, blhisp, totchr, ins, income)
#' step2(YS, XS, YO, XO)
#'
#' @importFrom Rdpack reprompt
#' @references {
#' \insertAllCited{}
#' }
#' @export
step2 <- function(YS, XS, YO, XO) {
  # Step 1: Estimate selection equation using probit model
  fit1 <- glm(YS ~ XS - 1, family = binomial(link = "probit"))

  # Step 2: Compute Inverse Mills Ratio (IMR)
  IMR <- dnorm(fit1$linear.predictors) / pnorm(fit1$linear.predictors)

  # Step 3: Estimate outcome equation with IMR as additional covariate
  xMat <- cbind(XO, IMR)
  fit2 <- lm(YO ~ XO + IMR - 1, subset = YS == 1)

  # Step 4: Create a data frame for convenience
  xMat <- data.frame(xMat)

  # Step 5: Compute delta values (used for variance correction)
  delta <- xMat$IMR * (xMat$IMR + fit1$fitted.values)

  # Step 6: Estimate corrected variance of the outcome
  Var <- (sum((fit2$residuals)^2) +
            (coef(fit2)[length(coef(fit2))]^2) * sum(delta[YS == 1])) / sum(YS)

  # Step 7: Estimate correlation coefficient (rho) between selection and outcome equations
  cor <- coef(fit2)[length(coef(fit2))] / sqrt(Var)
  names(cor) <- "rho"

  # Step 8: Compute standard deviation (sigma) of the outcome equation
  sigma <- sqrt(Var)
  names(sigma) <- "sigma"

  # Step 9: Combine all estimated parameters into a single vector
  beta <- c(coef(fit1), coef(fit2)[-length(coef(fit2))], sigma, cor)

  return(beta)
}
