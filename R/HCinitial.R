#' Two-Step Method for Parameter Estimation of the Heckman Model
#'
#' @description
#' Estimates the parameters of the classic
#' Heckman model via the two-step method.
#'
#' @return
#'
#' Returns a numerical vector with estimates
#' of the parameters of the classical Heckman
#' model using the two-step method
#'
#' @details
#' Generally, the two-step method is very
#' useful for finding initial values for
#' the Likelihood Estimation method. In
#' first step performs a probit analysis
#' on a selection equation. The second
#' step analyzes an outcome equation based
#' on the first-step binary probit model.
#'
#'
#' @param selection Selection equation.
#' @param outcome Primary Regression Equation.
#' @param data Database.
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' HCinitial(selectEq,outcomeEq, data = MEPS2001)
#'
#' @export
HCinitial <- function(selection, outcome, data = sys.frame(sys.parent())) {
    # Extrair matriz do modelo e matriz das equações de seleção e regressão Matriz
    # de seleção
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("selection", "data", "subset"), names(mf), 0)
    mfs <- mf[c(1, m)]
    mfs$drop.unused.levels <- TRUE
    mfs$na.action <- na.pass
    mfs[[1]] <- as.name("model.frame")
    names(mfs)[2] <- "formula"
    # model.frame requires the parameter to be 'formula'
    mfs <- eval(mfs, parent.frame())
    mts <- terms(mfs)
    xs <- model.matrix(mts, mfs)
    ys <- model.response(mfs)
    yslevels <- levels(as.factor(ys))
    ############### matriz de regressão
    m <- match(c("outcome", "data", "subset", "weights", "offset"), names(mf), 0)
    mfo <- mf[c(1, m)]
    mfo$na.action <- na.pass
    mfo$drop.unused.levels <- TRUE
    mfo$na.action <- na.pass
    mfo[[1]] <- as.name("model.frame")
    names(mfo)[2] <- "formula"
    mfo <- eval(mfo, parent.frame())
    mto <- attr(mfo, "terms")
    xo <- model.matrix(mto, mfo)
    yo <- model.response(mfo)
    # two-step estimation
    fit1 <- glm(ys ~ xs - 1, family = binomial(link = "probit"))
    imr <- dnorm(fit1$linear.predictors)/pnorm(fit1$linear.predictors)
    xmat <- cbind(xo, imr)
    fit2 <- lm(yo ~ xo + imr - 1, subset = ys == 1)
    xmat <- data.frame(xmat)
    # geracao de valores da nova covariavel delta
    delta <- (xmat$imr) * (xmat$imr + fit1$fitted.values)
    # calculo da variancia de y1
    var <- (sum((fit2$residuals)^2) + (((coef(fit2)[length(coef(fit2))])^2) * sum(delta[ys ==
        1])))/(sum(ys))
    # calculo da correlacao entre y1 e y2
    cor <- (coef(fit2)[length(coef(fit2))])/sqrt(var)
    names(cor) <- "rho"
    sigma <- sqrt(var)
    names(sigma) <- "sigma"
    # estimativas do modelo de heckman classico via metodo de dois passos
    result <- c(coef(fit1), coef(fit2)[-length(coef(fit2))], sigma, cor)
    return(result)
}
