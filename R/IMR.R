#' Inverse Mills Ratio
#'
#' @description
#' Column vector of the inverse ratio of Mills
#' @return
#' Return column vector of the inverse ratio of Mills
#' @param selection Selection equation.
#' @param data Database.
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' IMR(selectEq)
#' @export
IMR <- function(selection, data = sys.frame(sys.parent())) {
    ############################################################################################################################################## Extrair matriz do modelo e matriz da equação de seleção # Matriz de seleção #
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("selection", "data", "subset"), names(mf), 0)
    mfS <- mf[c(1, m)]
    mfS$drop.unused.levels <- TRUE
    mfS$na.action <- na.pass
    mfS[[1]] <- as.name("model.frame")
    names(mfS)[2] <- "formula"
    mfS <- eval(mfS, parent.frame())
    mtS <- terms(mfS)
    XS <- model.matrix(mtS, mfS)
    NXS <- ncol(XS)
    YS <- model.response(mfS)
    YSLevels <- levels(as.factor(YS))
    ############################################################################################################################################## Ajuste do modelo Probit #
    fit1 <- glm(YS ~ XS - 1, family = binomial(link = "probit"))
    ############################################################################################################################################## Cálculo da razão Inversa de Mills #
    IMR <- dnorm(fit1$linear.predictors)/pnorm(fit1$linear.predictors)
    return(cbind(IMR))
}
