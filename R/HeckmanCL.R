#' Classic Heckman Model fit Function
#'
#' @description
#' Estimates the parameters of the classic Heckman model
#' via Maximum Likelihood method. The initial start is obtained
#' via the two-step method.
#'
#' @return
#' Returns a list with the following components.
#'
#' Coefficients: Returns a numerical vector with the best estimated values
#' of the model parameters;
#'
#' Value: The value of function to be minimized (or maximized) corresponding
#' to par.
#'
#' loglik: Negative of value. Minimum (or maximum) of the likelihood function
#' calculated from the estimated coefficients.
#'
#' counts: Component of the Optim function. A two-element integer vector
#' giving the number of calls to fn and gr respectively. This excludes
#' those calls needed to compute the Hessian, if requested, and any calls
#' to fn to compute a finite-difference approximation to the gradient.
#'
#' hessian: Component of the Optim function, with pre-defined option
#' hessian=TRUE. A symmetric matrix giving an estimate of the Hessian
#' at the solution found. Note that this is the Hessian of the unconstrained
#' problem even if the box constraints are active.
#'
#' fisher_infoHC: Fisher information matrix
#'
#' prop_sigmaHC: Square root of the Fisher information matrix diagonal
#'
#' level: Selection variable levels
#'
#' nObs: Numeric value representing the size of the database
#'
#' nParam: Numerical value representing the number of model parameters
#'
#' N0: Numerical value representing the number of unobserved entries
#'
#' N1: Numerical value representing the number of complete entries
#'
#' NXS: Numerical value representing the number of parameters of the
#' selection model
#'
#' NXO: Numerical value representing the number of parameters of the
#' regression model
#'
#' df: Numerical value that represents the difference between the size
#' of the response vector of the selection equation and the number of
#' model parameters
#'
#' aic: Numerical value representing Akaike's information criterion.
#'
#' bic: Numerical value representing Schwarz's Bayesian Criterion
#'
#' initial.value: Numerical vector that represents the input values
#' (Initial Values) used in the parameter estimation.
#'
#' @param selection Selection equation.
#' @param outcome Primary Regression Equation.
#' @param start initial values.
#' @param data Database.
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' HeckmanCL(selectEq, outcomeEq, data = MEPS2001)
#' @export HeckmanCL
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
