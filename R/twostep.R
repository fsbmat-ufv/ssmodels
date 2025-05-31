#' Two-Step Estimation of the Classic Heckman Model
#'
#' @description
#' Estimates the parameters of the classical Heckman sample selection model
#' using the two-step procedure.
#'
#' @details
#' The two-step method first estimates a Probit model for the selection equation,
#' then fits an outcome equation that includes the Inverse Mills Ratio (IMR)
#' as an additional regressor to correct for sample selection bias.
#'
#' @param selection A formula for the selection equation.
#' @param outcome A formula for the outcome equation.
#' @param data A data frame containing the variables.
#'
#' @return
#' A numeric vector containing:
#' \itemize{
#'   \item Estimated coefficients of the selection equation (Probit model),
#'   \item Estimated coefficients of the outcome equation (excluding IMR),
#'   \item Estimated standard deviation of the outcome errors (\code{phi}),
#'   \item Estimated correlation between the error terms (\code{cor}).
#' }
#'
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' twostep(selectEq, outcomeEq, data = MEPS2001)
#'
#' @references
#' \insertAllCited{}
#' For details, see \insertCite{heckman1979sample;textual}{ssmodels}.
#'
#' @importFrom Rdpack reprompt
#' @export
twostep <- function(selection, outcome, data = sys.frame(sys.parent())) {
    ############################################################################################################################################## Extrair matriz do modelo e matriz das equações de seleção e regressão # Matriz
    ############################################################################################################################################## de seleção
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("selection", "data", "subset"), names(mf), 0)
    mfS <- mf[c(1, m)]
    mfS$drop.unused.levels <- TRUE
    mfS$na.action <- na.pass
    mfS[[1]] <- as.name("model.frame")
    names(mfS)[2] <- "formula"
    # model.frame requires the parameter to be 'formula'
    mfS <- eval(mfS, parent.frame())
    mtS <- terms(mfS)
    XS <- model.matrix(mtS, mfS)
    NXS <- ncol(XS)
    YS <- model.response(mfS)
    YSLevels <- levels(as.factor(YS))
    ############################################################################################################################################## Matriz de regressão #
    m <- match(c("outcome", "data", "subset", "weights", "offset"), names(mf), 0)
    mfO <- mf[c(1, m)]
    mfO$na.action <- na.pass
    mfO$drop.unused.levels <- TRUE
    mfO$na.action <- na.pass
    mfO[[1]] <- as.name("model.frame")
    names(mfO)[2] <- "formula"
    mfO <- eval(mfO, parent.frame())
    mtO <- attr(mfO, "terms")
    XO <- model.matrix(mtO, mfO)
    NXO <- ncol(XO)
    YO <- model.response(mfO)
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
    # Chute inicial para optim do modelo BS
    beta <- c(coef(fit1), coef(fit2)[-length(coef(fit2))], phi = sqrt(Var), cor)
    return(beta)
}
