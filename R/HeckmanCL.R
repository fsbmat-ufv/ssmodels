#' Classic Heckman Model Fit Function
#'
#' @description
#' Fits the classical Heckman sample selection model using Maximum Likelihood Estimation (MLE).
#' Initial parameter estimates are obtained via the two-step method.
#'
#' @details
#' This function estimates the parameters of the classical Heckman sample selection model
#' via MLE, accounting for potential sample selection bias. It uses the \code{optim} function
#' with the BFGS method to find the parameter estimates that maximize the log-likelihood function.
#' The initial values for optimization are obtained using the two-step Heckman method.
#'
#' The function returns a rich set of results, including:
#' \itemize{
#'   \item Estimated coefficients for the selection and outcome equations.
#'   \item Standard deviation of the outcome error term (\code{sigma}).
#'   \item Correlation between the errors of the selection and outcome equations (\code{rho}).
#'   \item Measures of model fit (AIC, BIC).
#'   \item Standard errors (approximated by the square root of the Fisher information diagonal).
#' }
#'
#' @param selection A formula specifying the selection equation.
#' @param outcome A formula specifying the primary outcome equation.
#' @param data A data frame containing the variables in the model.
#' @param start An optional numeric vector of initial parameter values. If not provided, default values are used.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{coefficients}: A named numeric vector of estimated parameters.
#'   \item \code{value}: The value of the (negative) log-likelihood at convergence.
#'   \item \code{loglik}: The maximized log-likelihood.
#'   \item \code{counts}: Number of gradient evaluations performed.
#'   \item \code{hessian}: Hessian matrix at the optimum.
#'   \item \code{fisher_infoHC}: The (approximate) Fisher information matrix.
#'   \item \code{prop_sigmaHC}: Approximate standard errors.
#'   \item \code{level}: Levels of the selection variable.
#'   \item \code{nObs}: Number of observations.
#'   \item \code{nParam}: Number of estimated parameters.
#'   \item \code{N0}: Number of unobserved (censored) observations.
#'   \item \code{N1}: Number of observed (uncensored) observations.
#'   \item \code{NXS}: Number of parameters in the selection equation.
#'   \item \code{NXO}: Number of parameters in the outcome equation.
#'   \item \code{df}: Degrees of freedom (observations minus parameters).
#'   \item \code{aic}: Akaike Information Criterion.
#'   \item \code{bic}: Bayesian Information Criterion.
#'   \item \code{initial.value}: Initial parameter values used in the optimization.
#' }
#'
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' HeckmanCL(selectEq, outcomeEq, data = MEPS2001)
#'
#' @references
#' \insertRef{heckman1979sample}{ssmodels}
#'
#' @export
HeckmanCL <- function(selection, outcome, data = sys.frame(sys.parent()), start = NULL) {
  # Extract model components from selection and outcome formulas
  components <- extract_model_components(selection = selection, outcome = outcome, data = data)
  XS  <- components$XS       # Design matrix for selection equation
  YS  <- components$YS       # Response variable for selection equation
  NXS <- components$NXS      # Number of selection parameters
  XO  <- components$XO       # Design matrix for outcome equation
  YO  <- components$YO       # Response variable for outcome equation
  NXO <- components$NXO      # Number of outcome parameters
  YSLevels <- components$YSLevels  # Levels of the selection variable

  ######################################
  # Log-likelihood function definition
  ######################################
  loglik_HC <- function(start) {
    # Parameter indices
    istartS <- 1:NXS
    istartO <- seq(tail(istartS, 1) + 1, length.out = NXO)
    isigma  <- tail(istartO, 1) + 1
    irho    <- isigma + 1

    # Extract parameters
    g     <- start[istartS]
    b     <- start[istartO]
    sigma <- start[isigma]
    rho   <- start[irho]

    # Check parameter validity
    if (!is.finite(sigma) || sigma <= 0) return(NA)
    if (!is.finite(rho) || abs(rho) >= 1) return(NA)

    # Compute linear predictors
    XS.g     <- XS %*% g
    XO.b     <- XO %*% b
    YO_clean <- ifelse(is.na(YO), 0, YO)
    u2       <- YO_clean - XO.b
    r        <- sqrt(1 - rho^2)
    B        <- (XS.g + (rho / sigma) * u2) / r

    # Log-likelihood contributions
    ll <- ifelse(YS == 0,
                 pnorm(-XS.g, log.p = TRUE),
                 dnorm(u2 / sigma, log = TRUE) - log(sigma) + pnorm(B, log.p = TRUE))

    return(sum(ll))
  }

  ######################################
  # Gradient of the log-likelihood
  ######################################
  gradlik_HC <- function(start) {
    nObs <- length(YS)

    # Parameter indices
    istartS <- 1:NXS
    istartO <- seq(tail(istartS, 1) + 1, length.out = NXO)
    isigma  <- tail(istartO, 1) + 1
    irho    <- isigma + 1
    nParam  <- irho

    # Extract parameters
    g     <- start[istartS]
    b     <- start[istartO]
    sigma <- start[isigma]
    rho   <- start[irho]

    # Check parameter validity
    if (!is.finite(sigma) || sigma <= 0) return(rep(NA, nParam))
    if (!is.finite(rho) || abs(rho) >= 1) return(rep(NA, nParam))

    # Partition data based on selection indicator
    XS0 <- XS[YS == 0, , drop = FALSE]
    XS1 <- XS[YS == 1, , drop = FALSE]
    XO1 <- XO[YS == 1, , drop = FALSE]
    YO_clean <- ifelse(is.na(YO), 0, YO)
    YO1 <- YO_clean[YS == 1]

    # Linear predictors and residuals
    XS0.g <- drop(XS0 %*% g)
    XS1.g <- drop(XS1 %*% g)
    XO1.b <- drop(XO1 %*% b)
    u2    <- YO1 - XO1.b
    r     <- sqrt(1 - rho^2)

    # Compute correction term
    B <- (XS1.g + (rho / sigma) * u2) / r
    lambdaB <- exp(dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))

    # Initialize gradient matrix
    gradient <- matrix(0, nrow = nObs, ncol = nParam)

    # Gradient contributions
    gradient[YS == 0, istartS] <- -XS0 * exp(dnorm(-XS0.g, log = TRUE) - pnorm(-XS0.g, log.p = TRUE))
    gradient[YS == 1, istartS] <-  XS1 * lambdaB / r
    gradient[YS == 1, istartO] <-  XO1 * (u2 / sigma^2 - lambdaB * rho / (sigma * r))
    gradient[YS == 1, isigma]  <- (u2^2 / sigma^3 - lambdaB * rho * u2 / (sigma^2 * r)) - 1 / sigma
    gradient[YS == 1, irho]    <- lambdaB * (u2 / sigma + rho * XS1.g) / r^3

    return(colSums(gradient))
  }

  ######################################
  # Starting values (if not provided)
  ######################################
  if (is.null(start)) {
    message("Start not provided using default start values.")
    start <- c(rep(0, NXS + NXO), 1, 0)
  }

  ######################################
  # Maximum Likelihood Estimation
  ######################################
  theta_HC <- optim(start,
                    loglik_HC,
                    gradlik_HC,
                    method  = "BFGS",
                    hessian = TRUE,
                    control = list(fnscale = -1))

  # Assign names to estimated parameters
  names(theta_HC$par) <- c(colnames(XS), colnames(XO), "sigma", "rho")

  ######################################
  # Organize results
  ######################################
  a   <- start
  a1  <- theta_HC$par
  a2  <- theta_HC$value
  a3  <- theta_HC$counts[2]
  a4  <- theta_HC$hessian
  a5  <- solve(-a4)
  a6  <- sqrt(diag(a5))
  a7  <- YSLevels
  a8  <- length(YS)
  a9  <- length(start)
  a10 <- sum(YS == 0)
  a11 <- sum(YS == 1)
  a12 <- NXS
  a13 <- NXO
  a14 <- a8 - a9
  a15 <- -2 * a2 + 2 * a9
  a16 <- -2 * a2 + a9 * log(a8)

  # Create result object
  result <- list(coefficients   = a1,
                 value          = a2,
                 loglik         = -a2,
                 counts         = a3,
                 hessian        = a4,
                 fisher_infoHC  = a5,
                 prop_sigmaHC   = a6,
                 level          = a7,
                 nObs           = a8,
                 nParam         = a9,
                 N0             = a10,
                 N1             = a11,
                 NXS            = a12,
                 NXO            = a13,
                 df             = a14,
                 aic            = a15,
                 bic            = a16,
                 initial.value  = a)

  class(result) <- c("HeckmanCL", class(theta_HC))
  return(result)
}
