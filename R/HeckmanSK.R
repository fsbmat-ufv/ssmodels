#' Skew-Normal Sample Selection Model Fit Function
#'
#' @description
#' Fits a sample selection model based on the Skew-Normal distribution
#' using Maximum Likelihood Estimation (MLE). This model allows for
#' asymmetry in the distribution of the outcome variable's error term,
#' addressing potential skewness.
#'
#' @details
#' The function implements MLE for a sample selection model where
#' the outcome equation's errors follow a Skew-Normal distribution,
#' as proposed in \insertCite{ogundimu2016sample;textual}{ssmodels}.
#' The optimization is performed via the BFGS algorithm.
#'
#' The results include estimates for:
#' \itemize{
#'   \item Selection equation coefficients.
#'   \item Outcome equation coefficients.
#'   \item Standard deviation of the error term (\code{sigma}).
#'   \item Correlation between the selection and outcome errors (\code{rho}).
#'   \item Skewness parameter (\code{lambda}).
#'   \item Robust standard errors from the Fisher information matrix.
#' }
#'
#' @param selection A formula specifying the selection equation.
#' @param outcome A formula specifying the outcome equation.
#' @param lambda Initial start value for the skewness parameter (\code{lambda}).
#' @param start Optional numeric vector of initial parameter values.
#' @param data A data frame containing the variables.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{coefficients}: Named vector of estimated model parameters.
#'   \item \code{value}: The (negative) log-likelihood at convergence.
#'   \item \code{loglik}: The maximum log-likelihood.
#'   \item \code{counts}: Number of gradient evaluations.
#'   \item \code{hessian}: Hessian matrix at the optimum.
#'   \item \code{fisher_infoSK}: Approximate Fisher information matrix.
#'   \item \code{prop_sigmaSK}: Standard errors for the estimates.
#'   \item \code{level}: Levels of the selection variable.
#'   \item \code{nObs}: Number of observations.
#'   \item \code{nParam}: Number of model parameters.
#'   \item \code{N0}: Number of censored (unobserved) observations.
#'   \item \code{N1}: Number of observed (uncensored) observations.
#'   \item \code{NXS}: Number of covariates in the selection equation.
#'   \item \code{NXO}: Number of covariates in the outcome equation.
#'   \item \code{df}: Degrees of freedom (observations minus parameters).
#'   \item \code{aic}: Akaike Information Criterion.
#'   \item \code{bic}: Bayesian Information Criterion.
#'   \item \code{initial.value}: Initial parameter values used.
#' }
#'
#' @examples
#' data("Mroz87")
#' attach(Mroz87)
#' selectEq <- lfp ~ huswage + kids5 + mtr + fatheduc + educ + city
#' outcomeEq <- log(wage) ~ educ + city
#' HeckmanSK(selectEq, outcomeEq, data = Mroz87, lambda = -1.5)
#'
#' @references
#' \insertRef{ogundimu2016sample}{ssmodels}
#'
#' @importFrom Rdpack reprompt
#' @export
HeckmanSK <- function(selection, outcome, data = sys.frame(sys.parent()), lambda, start = NULL) {

  # Extract design matrices and outcomes from the selection and outcome models
  components <- extract_model_components(selection = selection, outcome = outcome, data = data)
  XS  <- components$XS       # Design matrix for selection equation
  YS  <- components$YS       # Response variable for selection equation
  NXS <- components$NXS      # Number of selection parameters
  XO  <- components$XO       # Design matrix for outcome equation
  YO  <- components$YO       # Response variable for outcome equation
  NXO <- components$NXO      # Number of outcome parameters
  YSLevels <- components$YSLevels  # Levels of the selection variable

  # Replace missing values in the outcome with zeros
  YO[is.na(YO)] <- 0

  #############################
  # Log-likelihood function
  #############################
  loglik_SK <- function(start) {
    # Define parameter indices
    NXS <- dim(model.matrix(~XS))[2] - 1
    NXO <- dim(model.matrix(~XO))[2] - 1
    istartS <- 1:NXS
    istartO <- seq(tail(istartS, 1) + 1, length = NXO)
    isigma  <- tail(istartO, 1) + 1
    irho    <- tail(isigma, 1) + 1
    ilamb1  <- tail(irho, 1) + 1

    # Extract parameters
    g      <- start[istartS]
    b      <- start[istartO]
    sigma  <- start[isigma]
    if (sigma < 0) return(NA)
    rho    <- start[irho]
    if ((rho < -1) || (rho > 1)) return(NA)
    lamb1  <- start[ilamb1]

    # Linear predictors
    XS.g <- XS %*% g
    XO.b <- XO %*% b
    u2   <- YO - XO.b
    z    <- u2 / sigma
    r    <- sqrt(1 - rho^2)
    u    <- 1 + lamb1^2 - (lamb1 * rho)^2
    lstar <- -lamb1 * rho / sqrt(u)

    # Compute selection log-likelihood
    h <- sn::psn(-XS.g, 0, 1, -lstar)

    # Log-likelihood vector
    ll <- ifelse(YS == 0,
                 log(h),
                 log(2 / sigma) + log(dnorm(z)) + log(pnorm(lamb1 * z)) +
                   log(pnorm((XS.g + rho * z) / r)))

    return(sum(ll))
  }

  #############################
  # Gradient of the log-likelihood
  #############################
  gradlik_SK <- function(start) {
    # Redefine dimensions and parameter counts
    NXS <- dim(model.matrix(~XS))[2] - 1
    NXO <- dim(model.matrix(~XO))[2] - 1
    nObs <- length(YS)
    nParam <- NXS + NXO + 3

    # Subset by selection
    XS0 <- XS[YS == 0, , drop = FALSE]
    XS1 <- XS[YS == 1, , drop = FALSE]
    XO1 <- XO[YS == 1, , drop = FALSE]
    YO[is.na(YO)] <- 0
    YO1 <- YO[YS == 1]

    N0 <- sum(YS == 0)
    N1 <- sum(YS == 1)
    w0 <- rep(1, N0)
    w1 <- rep(1, N1)

    # Assign parameter indices
    istartS <- 1:NXS
    istartO <- seq(tail(istartS, 1) + 1, length = NXO)
    isigma  <- tail(istartO, 1) + 1
    irho    <- tail(isigma, 1) + 1
    ilamb1  <- tail(irho, 1) + 1

    # Extract parameters
    g      <- start[istartS]
    b      <- start[istartO]
    sigma  <- start[isigma]
    if (sigma < 0) return(matrix(NA, nObs, nParam))
    rho    <- start[irho]
    if ((rho < -1) || (rho > 1)) return(matrix(NA, nObs, nParam))
    lamb1  <- start[ilamb1]

    # Linear predictors
    XS0.g <- as.numeric(XS0 %*% g)
    XS1.g <- as.numeric(XS1 %*% g)
    XO1.b <- as.numeric(XO1 %*% b)
    u2    <- YO1 - XO1.b
    z     <- u2 / sigma
    r     <- 1 - rho^2
    u     <- 1 + lamb1^2 - (lamb1 * rho)^2
    lstar <- -lamb1 * rho / sqrt(u)

    # Derivatives of lstar w.r.t rho and lambda
    dlrho <- (-rho * u + lamb1 * rho * (lamb1 - lamb1 * rho^2)) / u^(3/2)
    dllam <- (-rho / sqrt(u)) + lamb1 * rho * (lamb1 - lamb1 * rho^2) / (sqrt(u) * u)

    # Intermediate quantities
    omeg <- (XS1.g + rho * z) / sqrt(r)
    K1   <- dnorm(omeg) / pnorm(omeg)
    h    <- function(t) sn::psn(-t, 0, 1, -lstar)
    h2   <- function(t) sn::psn(-t, 0, 1, lamb1 * rho / sqrt(u))
    K2   <- dnorm(-XS0.g) * pnorm(lstar * XS0.g) / h(XS0.g)
    eta  <- lamb1 * z
    K3   <- dnorm(eta) / pnorm(eta)
    K4   <- dnorm(XS0.g * sqrt(1 + lamb1^2) / sqrt(u)) / h(XS0.g)

    # Initialize gradient matrix
    gradient <- matrix(0, nObs, nParam)

    # Compute gradients
    gradient[YS == 0, istartS] <- w0 * XS0 * (-2 * K2)
    gradient[YS == 1, istartS] <- w1 * XS1 * (K1 / sqrt(r))
    gradient[YS == 1, istartO] <- w1 * XO1 * ((z / sigma) - (lamb1 / sigma) * K3 - (rho / (sigma * sqrt(r))) * K1)
    gradient[YS == 1, isigma]  <- w1 * (-1 / sigma + (z^2) / sigma - (lamb1 * K3 * z) / sigma - (rho * K1 * z) / (sigma * sqrt(r)))
    gradient[YS == 0, irho]    <- w0 * (-(2 / sqrt(2 * pi)) * (lamb1 * (1 + lamb1^2) / (sqrt(u) * (u + (lamb1 * rho)^2))) *
                                          dnorm(sqrt(1 + (lamb1 * rho / sqrt(u))^2) * XS0.g) / h2(XS0.g))
    gradient[YS == 1, irho]    <- w1 * (K1 * (rho * XS1.g + z)) / (sqrt(r)^3)
    gradient[YS == 0, ilamb1]  <- w0 * dllam * sqrt(2 / pi) * (1 / (1 + lstar^2)) *
      dnorm(sqrt(1 + lstar^2) * XS0.g) / h(XS0.g)
    gradient[YS == 1, ilamb1]  <- w1 * (K3 * z)

    return(colSums(gradient))
  }

  #############################
  # Initial values
  #############################
  if (is.null(start)) {
    message("Start not provided using default start values.")
    start <- c(step2(YS, XS, YO, XO), lambda)
  }

  #############################
  # Optimization via BFGS
  #############################
  theta_SK <- optim(start, loglik_SK, gradlik_SK, method = "BFGS", hessian = TRUE,
                    control = list(fnscale = -1))

  #############################
  # Organize and return results
  #############################
  names(theta_SK$par) <- c(colnames(XS), colnames(XO), "sigma", "rho", "lambda")

  result <- list(
    coefficients   = theta_SK$par,
    value          = theta_SK$value,
    loglik         = -theta_SK$value,
    counts         = theta_SK$counts[2],
    hessian        = theta_SK$hessian,
    fisher_infoSK  = solve(-theta_SK$hessian),
    prop_sigmaSK   = sqrt(diag(solve(-theta_SK$hessian))),
    level          = YSLevels,
    nObs           = length(YS),
    nParam         = length(start),
    N0             = sum(YS == 0),
    N1             = sum(YS == 1),
    NXS            = ncol(XS),
    NXO            = ncol(XO),
    df             = length(YS) - length(start),
    aic            = -2 * theta_SK$value + 2 * length(start),
    bic            = -2 * theta_SK$value + length(start) * log(length(YS)),
    initial.value  = start
  )

  class(result) <- c("HeckmanSK", class(theta_SK))
  return(result)
}
