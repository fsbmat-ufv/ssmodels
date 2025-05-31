#' Heckman's Two-Step Method
#'
#' @description
#' Estimates the parameters of the classical Heckman selection model using the two-step method.
#' The first step fits a probit model for the selection equation. In the second step, the inverse Mills ratio (IMR)
#' is included as an additional regressor in the outcome equation.
#'
#' @details
#' This function implements the two-step estimation procedure of the classical Heckman model.
#' In the first step, a probit model is estimated to predict the selection indicator \code{YS} using
#' the selection covariates \code{XS}. The IMR is calculated from this model.
#' In the second step, an ordinary least squares (OLS) regression of the observed outcome \code{YO} on
#' \code{XO} and the IMR is performed for the uncensored observations (\code{YS == 1}).
#'
#' The function also calculates:
#' \itemize{
#'   \item \code{sigma}: The estimated standard deviation of the outcome equation's error term.
#'   \item \code{rho}: The estimated correlation between the error terms of the selection and outcome equations.
#' }
#'
#' @param YS A binary vector indicating selection (\code{1} if observed, \code{0} otherwise).
#' @param XS A matrix of covariates for the selection equation.
#' @param YO A numeric vector representing the outcome variable of interest.
#' @param XO A matrix of covariates for the outcome equation.
#'
#' @return
#' A numeric vector containing the parameter estimates from the two-step Heckman model:
#' \itemize{
#'   \item Coefficients of the selection equation (probit model).
#'   \item Coefficients of the outcome equation (excluding the IMR term).
#'   \item Estimated \code{sigma}.
#'   \item Estimated \code{rho}.
#' }
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
#' @references
#' \insertAllCited{}
#'
#' @export
step2 <- function(YS, XS, YO, XO) {
    # Two-step estimation
    fit1 <- glm(YS ~ XS - 1, family = binomial(link = "probit"))
    IMR <- dnorm(fit1$linear.predictors)/pnorm(fit1$linear.predictors)
    xMat <- cbind(XO, IMR)
    fit2 <- lm(YO ~ XO + IMR - 1, subset = YS == 1)
    xMat <- data.frame(xMat)
    # Geracao de valores da nova covariavel delta
    delta <- (xMat$IMR) * (xMat$IMR + fit1$fitted.values)
    # Calculo da variancia de Y1
    Var <- (sum((fit2$residuals)^2) + (((coef(fit2)[length(coef(fit2))])^2) * sum(delta[YS ==
        1])))/(sum(YS))
    # Calculo da correlacao entre Y1 e Y2
    cor <- (coef(fit2)[length(coef(fit2))])/sqrt(Var)
    names(cor) <- "rho"
    sigma <- sqrt(Var)
    names(sigma) <- "sigma"
    #Chute inicial para optim do modelo BS
    beta <- c(coef(fit1), coef(fit2)[-length(coef(fit2))], sigma, cor)
    return(beta)
}
