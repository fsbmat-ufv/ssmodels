#' Two-Step Method for Parameter Estimation of the Classical Heckman Model
#'
#' @description
#' Estimates the parameters of the classical Heckman sample selection model using the two-step estimation method.
#'
#' @details
#' This function implements the two-step approach proposed by Heckman (1979) to estimate the parameters
#' of the classic sample selection model. It is particularly useful for obtaining initial values
#' for maximum likelihood estimation (MLE).
#'
#' In the first step, a probit model is fitted to the selection equation to estimate the probability of selection.
#' The second step involves estimating a linear regression of the outcome equation for the observed (selected) data,
#' incorporating the inverse Mills ratio (IMR) as an additional regressor to correct for sample selection bias.
#'
#' The function also estimates:
#' \itemize{
#'   \item \code{sigma}: The standard deviation of the outcome equation's error term.
#'   \item \code{rho}: The correlation coefficient between the errors of the selection and outcome equations.
#' }
#'
#' @param selection A formula specifying the selection equation.
#' @param outcome A formula specifying the outcome equation.
#' @param data A data frame containing the variables in the model.
#'
#' @return A named numeric vector containing:
#' \itemize{
#'   \item Coefficients from the selection equation (probit model),
#'   \item Coefficients from the outcome equation (excluding the IMR),
#'   \item Estimated \code{sigma},
#'   \item Estimated \code{rho}.
#' }
#'
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' HCinitial(selectEq, outcomeEq, data = MEPS2001)
#'
#' @references
#' \insertRef{heckman1979sample}{ssmodels}
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
