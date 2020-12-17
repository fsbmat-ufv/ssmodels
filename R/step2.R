#' Heckman's two-step method
#'
#' Estimate model parameters via two-step method
#'
#' @param YS Selection vector.
#' @param XS Selection Matrix.
#' @param YO Interest vector.
#' @param XO Matrix of the equation of interest.
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
