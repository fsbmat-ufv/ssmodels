#' Heckman-t Model fit Function
#'
#' @description
#' Estimates the parameters of the Heckman-t model
#'
#' @return
#'
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
#' fisher_infotS: Fisher information matrix
#'
#' prop_sigmatS: Square root of the Fisher information matrix diagonal
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
#' @details
#' The HeckmantS() function fits the Sample Selection Model
#' based on the Student's t distribution. For more information see
#' \insertCite{marchenko2012heckman;textual}{ssmodels}
#'
#' @param selection Selection equation.
#' @param outcome Primary Regression Equation.
#' @param df Initial start to the degree of freedom.
#' @param start initial values.
#' @param data Database.
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' HeckmantS(selectEq, outcomeEq, data = MEPS2001, df=12)
#' @importFrom Rdpack reprompt
#' @references {
#' \insertAllCited{}
#' }
#' @export HeckmantS
#' @export
HeckmantS <- function(selection, outcome, data = sys.frame(sys.parent()), df, start = NULL) {
  # Step 1: Extract components from selection and outcome formulas
  components <- extract_model_components(selection = selection, outcome = outcome, data = data)
  XS  <- components$XS
  YS  <- components$YS
  NXS <- components$NXS
  XO  <- components$XO
  YO  <- components$YO
  NXO <- components$NXO
  YSLevels <- components$YSLevels

  # Step 2: Define the log-likelihood function for the Heckman-t model
  loglik_tS <- function(start) {
    # Compute the log-likelihood based on the Student's t-distribution
    # using the current parameter estimates.

    NXS <- dim(model.matrix(~XS))[2] - 1  # Number of columns in XS minus intercept
    NXO <- dim(model.matrix(~XO))[2] - 1  # Number of columns in XO minus intercept

    # Define parameter indices
    istartS <- 1:NXS
    istartO <- seq(tail(istartS, 1) + 1, length = NXO)
    isigma  <- tail(istartO, 1) + 1
    irho    <- tail(isigma, 1) + 1
    iv      <- tail(irho, 1) + 1

    # Extract parameter values
    g     <- start[istartS]
    b     <- start[istartO]
    sigma <- start[isigma]

    # Validate sigma
    if (sigma < 0) return(NA)

    rho <- start[irho]
    if ((rho < -1) || (rho > 1)) return(NA)

    v <- start[iv]

    # Compute linear predictors
    XS.g <- XS %*% g
    XO.b <- XO %*% b
    u2   <- YO - XO.b
    z    <- u2 / sigma
    r    <- sqrt(1 - rho^2)
    Q    <- sqrt((v + 1) / (v + z^2))
    eta  <- Q * ((rho * z + XS.g) / r)

    # Compute log-density constant
    gam <- log(gamma((v + 1) / 2)) - log(gamma(v / 2)) - 0.5 * log(pi) -
      0.5 * log(v) - log(sigma)

    # Compute log-likelihood contributions
    ll <- ifelse(YS == 0,
                 pt(-XS.g, v, log.p = TRUE),
                 gam - ((v + 1)/2) * log(1 + (z^2 / v)) + pt(eta, v + 1, log.p = TRUE))

    # Return total log-likelihood
    return(sum(ll))
  }
    ############## Gradient #
  # Step 3: Define the gradient of the log-likelihood function
  gradlik_tS <- function(start) {
    # Compute the analytical gradient of the log-likelihood
    # with respect to all model parameters.

    NXS <- dim(model.matrix(~XS))[2] - 1  # Number of columns in XS (excluding intercept)
    NXO <- dim(model.matrix(~XO))[2] - 1  # Number of columns in XO (excluding intercept)
    nObs <- length(YS)                   # Total number of observations
    NO <- length(YS[YS > 0])             # Number of non-zero observations
    nParam <- NXS + NXO + 3              # Total number of parameters
    N0 <- sum(YS == 0)                   # Number of censored observations
    N1 <- sum(YS == 1)                   # Number of uncensored observations

    w <- rep(1, N0 + N1)                 # Weight vector for all observations
    w0 <- rep(1, N0)                     # Weight vector for censored observations
    w1 <- rep(1, N1)                     # Weight vector for uncensored observations

    # Indices for parameter vector
    istartS <- 1:NXS
    istartO <- seq(tail(istartS, 1) + 1, length = NXO)
    isigma <- tail(istartO, 1) + 1
    irho <- tail(isigma, 1) + 1
    iv <- tail(irho, 1) + 1

    # Extract selection and outcome coefficients
    g <- start[istartS]
    b <- start[istartO]

    # Extract and transform scale (sigma), correlation (rho), and degrees of freedom (v)
    chuteS = start[isigma]
    lns = log(chuteS)
    sigma = exp(lns)

    chuteR <- start[irho]
    tau = log((1 + chuteR)/(1 - chuteR))/2
    rho = (exp(2 * tau) - 1)/(exp(2 * tau) + 1)

    chuteV = start[iv]
    lndf = log(chuteV)
    v = exp(lndf)

    # Subsets for censored (YS == 0) and uncensored (YS == 1) cases
    XS0 <- XS[YS == 0, , drop = FALSE]
    XS1 <- XS[YS == 1, , drop = FALSE]
    YO[is.na(YO)] <- 0
    YO1 <- YO[YS == 1]
    XO1 <- XO[YS == 1, , drop = FALSE]

    # Linear predictors
    XS0.g <- as.numeric((XS0) %*% g)
    XS1.g <- as.numeric((XS1) %*% g)
    XO1.b <- as.numeric((XO1) %*% b)

    # Residuals and transformations
    u2O <- YO1 - XO1.b
    z0 <- u2O / sigma
    r <- sqrt(1 - rho^2)
    Qv <- ((v + 1) / (v + z0^2))^(1/2)
    Ar <- 1 / sqrt(1 - rho^2)
    Arr <- rho * Ar
    QSI <- (Arr * z0) + Ar * XS1.g
    eta <- Qv * QSI
    dr = 4 * exp(2 * tau) / ((exp(2 * tau) + 1)^2)

    # Gradient with respect to v for censored cases
    f1 <- function(x) {
      ff <- pt(-XS0.g, x, log.p = TRUE)
      return(ff)
    }
    gv2 <- numDeriv::grad(f1, rep(1, length(XS0.g)) * v)

    # Gradient with respect to v for uncensored cases
    f2 <- function(x) {
      teste = (((x + 1) / (x + z0^2))^(1/2)) * QSI
      ff <- pt(teste, x + 1, log.p = TRUE)
      return(ff)
    }
    gv <- numDeriv::grad(f2, rep(1, length(z0)) * v)

    # Initialize gradient matrix
    gradient <- matrix(0, nObs, nParam)

    # Gradient for selection equation parameters (g)
    gradient[YS == 0, istartS] <- -w0 * (XS0) * (dt(-XS0.g, v) / pt(-XS0.g, v))
    gradient[YS == 1, istartS] <- w1 * (XS1) * Ar * Qv * (dt(eta, v + 1) / pt(eta, v + 1))

    # Gradient for outcome equation parameters (b)
    gradient[YS == 1, istartO] <- w1 * (XO1) * (Qv / sigma) *
      ((Qv * z0) + ((QSI * ((v + z0^2)^(-1)) * z0 - Arr) * (dt(eta, v + 1) / pt(eta, v + 1))))

    # Gradient for sigma
    gradient[YS == 1, isigma] <- w1 * (-1 + (Qv * z0)^2 +
                                         (dt(eta, v + 1) / pt(eta, v + 1)) * ((Qv * z0) * (QSI * ((v + z0^2)^(-1)) * z0 - Arr)))

    # Gradient for rho
    gradient[YS == 1, irho] <- w1 * (dt(eta, v + 1) / pt(eta, v + 1)) *
      Qv * (Ar^3) * (z0 + rho * XS1.g) * dr

    # Gradient for v (degrees of freedom)
    gradient[YS == 1, iv] <- w1 * ((((1/2) * digamma((v + 1)/2) - (1/2) * digamma(v/2)) -
                                      (1 / (2 * v))) - ((1/2) * log(1 + ((z0^2) / v))) +
                                     ((Qv * z0)^2) / (2 * v) + gv)

    gradient[YS == 0, iv] <- w0 * gv2

    # Return sum of gradients for each parameter
    return(colSums(gradient))
  }


  # Step 4: Define initial values if not provided
  if (is.null(start)) {
    message("Start not provided using default start values.")
    start <- c(rep(0, NXS + NXO), 1, 0, log(df))
  }

  # Step 5: Maximize the log-likelihood using BFGS optimization
  theta_tS <- optim(start, loglik_tS, gradlik_tS, method = "BFGS", hessian = T,
                    control = list(fnscale = -1))

  # Step 6: Assign names to estimated parameters
  names(theta_tS$par) <- c(colnames(XS), colnames(XO), "sigma", "rho", "df")

  # Step 7: Construct output list with all relevant model components
  result <- list(
    coefficients    = theta_tS$par,              # Estimated coefficients
    value           = theta_tS$value,            # Optim objective function value
    loglik          = -theta_tS$value,           # Log-likelihood
    counts          = theta_tS$counts[2],        # Number of iterations
    hessian         = theta_tS$hessian,          # Hessian matrix
    fisher_infotS   = solve(-theta_tS$hessian),  # Fisher Information Matrix
    prop_sigmatS    = sqrt(diag(solve(-theta_tS$hessian))), # Standard errors
    level           = YSLevels,                  # Selection levels
    nObs            = length(YS),                # Number of observations
    nParam          = length(start),             # Number of parameters
    N0              = sum(YS == 0),              # Number of censored obs
    N1              = sum(YS == 1),              # Number of uncensored obs
    NXS             = NXS,                       # Number of selection regressors
    NXO             = NXO,                       # Number of outcome regressors
    df              = length(YS) - length(start),# Degrees of freedom
    aic             = -2 * theta_tS$value + 2 * length(start), # AIC
    bic             = -2 * theta_tS$value + length(start) * log(length(YS)), # BIC
    initial.value   = start                      # Initial values used
  )
  class(result) <- c("HeckmantS", class(theta_tS))
  return(result)
}

