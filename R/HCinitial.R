#' Two-Step Method for Parameter Estimation of the Heckman Model
#'
#' @description
#' Estimates the parameters of the classical Heckman model using the two-step method.
#'
#' @details
#' This method is especially useful for obtaining initial values for maximum likelihood estimation.
#' The first step fits a Probit model to the selection equation.
#' The second step fits a linear regression to the outcome equation, adjusting for selection bias using the inverse Mills ratio.
#'
#' @param selection A formula for the selection equation.
#' @param outcome A formula for the outcome equation.
#' @param data A data frame containing the variables.
#'
#' @return A numeric vector containing the estimated parameters of the classic Heckman model:
#' \itemize{
#'   \item Coefficients from the selection equation (Probit model),
#'   \item Coefficients from the outcome equation (excluding the IMR),
#'   \item Estimated standard deviation of the error term (\code{sigma}),
#'   \item Estimated correlation between error terms (\code{rho}).
#' }
#'
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' HCinitial(selectEq, outcomeEq, data = MEPS2001)
#'
#' @export
HCinitial <- function(selection, outcome, data = sys.frame(sys.parent())) {

  # Step 1: Construct model frame for selection equation
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("selection", "data", "subset"), names(mf), 0)
  mfs <- mf[c(1, m)]
  mfs$drop.unused.levels <- TRUE
  mfs$na.action <- na.pass
  mfs[[1]] <- as.name("model.frame")
  names(mfs)[2] <- "formula"  # required argument name
  mfs <- eval(mfs, parent.frame())
  mts <- terms(mfs)
  xs <- model.matrix(mts, mfs)
  ys <- model.response(mfs)
  yslevels <- levels(as.factor(ys))

  # Step 2: Construct model frame for outcome equation
  m <- match(c("outcome", "data", "subset", "weights", "offset"), names(mf), 0)
  mfo <- mf[c(1, m)]
  mfo$drop.unused.levels <- TRUE
  mfo$na.action <- na.pass
  mfo[[1]] <- as.name("model.frame")
  names(mfo)[2] <- "formula"
  mfo <- eval(mfo, parent.frame())
  mto <- attr(mfo, "terms")
  xo <- model.matrix(mto, mfo)
  yo <- model.response(mfo)

  # Step 3: Fit Probit model for the selection equation
  fit1 <- glm(ys ~ xs - 1, family = binomial(link = "probit"))

  # Step 4: Compute Inverse Mills Ratio (IMR)
  imr <- dnorm(fit1$linear.predictors) / pnorm(fit1$linear.predictors)

  # Step 5: Fit outcome regression using IMR as additional covariate
  xmat <- cbind(xo, imr)
  fit2 <- lm(yo ~ xo + imr - 1, subset = ys == 1)

  # Step 6: Estimate additional quantities for sigma and rho
  xmat <- data.frame(xmat)
  delta <- xmat$imr * (xmat$imr + fit1$fitted.values)
  resid_squared_sum <- sum((fit2$residuals)^2)
  imr_coef <- coef(fit2)[length(coef(fit2))]
  var <- (resid_squared_sum + (imr_coef^2) * sum(delta[ys == 1])) / sum(ys)
  sigma <- sqrt(var)
  names(sigma) <- "sigma"

  # Step 7: Estimate rho
  rho <- imr_coef / sigma
  names(rho) <- "rho"

  # Step 8: Return parameter estimates
  result <- c(coef(fit1), coef(fit2)[-length(coef(fit2))], sigma, rho)
  return(result)
}
