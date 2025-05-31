#' Heckman-t Model Fit Function
#'
#' @description
#' Fits a sample selection model based on the Student's t-distribution,
#' extending the classical Heckman model to account for heavy-tailed error terms.
#' The estimation is performed via Maximum Likelihood using the BFGS algorithm.
#'
#' @details
#' The function implements the Heckman sample selection model using
#' the Student's t-distribution for the error terms, as proposed by
#' \insertCite{marchenko2012heckman;textual}{ssmodels}. This extension
#' allows for robustness against outliers and heavy-tailed distributions.
#' Initial parameter values can be specified by the user or default to standard starting values.
#'
#' @param selection A formula specifying the selection equation.
#' @param outcome A formula specifying the outcome equation.
#' @param df Initial value for the degrees of freedom parameter of the t-distribution.
#' @param start Optional numeric vector of initial parameter values.
#' @param data A data frame containing the variables in the model.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{coefficients}: Named vector of estimated model parameters.
#'   \item \code{value}: Negative of the maximum log-likelihood.
#'   \item \code{loglik}: Maximum log-likelihood.
#'   \item \code{counts}: Number of gradient evaluations performed.
#'   \item \code{hessian}: Hessian matrix at the optimum.
#'   \item \code{fisher_infotS}: Approximate Fisher information matrix.
#'   \item \code{prop_sigmatS}: Standard errors for the parameter estimates.
#'   \item \code{level}: Levels of the selection variable.
#'   \item \code{nObs}: Number of observations.
#'   \item \code{nParam}: Number of model parameters.
#'   \item \code{N0}: Number of censored (unobserved) observations.
#'   \item \code{N1}: Number of uncensored (observed) observations.
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
#' HeckmantS(selectEq, outcomeEq, data = MEPS2001, df = 12)
#'
#' @references
#' \insertRef{marchenko2012heckman}{ssmodels}
#'
#' @importFrom Rdpack reprompt
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
  ################### Likelihood #
  loglik_tS <- function(start) {
    NXS <- dim(model.matrix(~XS))[2] - 1  #Numero de colunas de XS+1
    NXO <- dim(model.matrix(~XO))[2] - 1  #Numero de colunas de XO+1
    ## parameter indices
    istartS <- 1:NXS
    istartO <- seq(tail(istartS, 1) + 1, length = NXO)
    isigma <- tail(istartO, 1) + 1
    irho <- tail(isigma, 1) + 1
    iv <- tail(irho, 1) + 1
    g <- start[istartS]
    b <- start[istartO]
    sigma <- start[isigma]
    if (sigma < 0)
      return(NA)
    rho <- start[irho]
    if ((rho < -1) || (rho > 1))
      return(NA)
    v <- start[iv]
    XS.g <- (XS) %*% g
    XO.b <- (XO) %*% b
    u2 <- YO - XO.b
    z <- u2/sigma
    r <- sqrt(1 - rho^2)
    Q <- ((v + 1)/(v + (z^2)))^(1/2)
    eta <- Q * ((rho * z + XS.g)/r)
    gam <- log(gamma((v + 1)/2)) - log(gamma(v/2)) - 0.5 * log(pi) - 0.5 * log(v) -
      log(sigma)
    ll <- ifelse(YS == 0, (pt(-XS.g, v, log.p = TRUE)), (gam - ((v + 1)/2) *
                                                           log(1 + ((z^2)/v)) + pt(eta, v + 1, log.p = TRUE)))
    return(sum(ll))
  }
  ############## Gradient #
  gradlik_tS <- function(start) {
    NXS <- dim(model.matrix(~XS))[2] - 1  #Numero de colunas de XS+1
    NXO <- dim(model.matrix(~XO))[2] - 1  #Numero de colunas de XO+1
    nObs <- length(YS)
    NO <- length(YS[YS > 0])
    nParam <- NXS + NXO + 3  #Total of parameters
    N0 <- sum(YS == 0)
    N1 <- sum(YS == 1)

    w <- rep(1, N0 + N1)
    w0 <- rep(1, N0)
    w1 <- rep(1, N1)

    ## parameter indices
    istartS <- 1:NXS
    istartO <- seq(tail(istartS, 1) + 1, length = NXO)
    isigma <- tail(istartO, 1) + 1
    irho <- tail(isigma, 1) + 1
    iv <- tail(irho, 1) + 1
    g <- start[istartS]
    b <- start[istartO]

    chuteS = start[isigma]
    lns = log(chuteS)
    sigma = exp(lns)
    chuteR <- start[irho]
    tau = log((1 + chuteR)/(1 - chuteR))/2
    rho = (exp(2 * tau) - 1)/(exp(2 * tau) + 1)
    chuteV = start[iv]
    lndf = log(chuteV)
    v = exp(lndf)

    XS0 <- XS[YS == 0, , drop = FALSE]
    XS1 <- XS[YS == 1, , drop = FALSE]
    YO[is.na(YO)] <- 0
    YO1 <- YO[YS == 1]
    XO1 <- XO[YS == 1, , drop = FALSE]
    XS0.g <- as.numeric((XS0) %*% g)
    XS1.g <- as.numeric((XS1) %*% g)
    XO1.b <- as.numeric((XO1) %*% b)
    # u2 <- YO1 - XO1.b u2S1 <- YO1 - XS1.g
    u2O <- YO1 - XO1.b
    # u2 <- YO - XO.b
    z0 <- u2O/sigma
    r <- sqrt(1 - rho^2)
    Qv <- ((v + 1)/(v + (z0)^2))^(1/2)
    Ar <- 1/sqrt(1 - (rho^2))
    Arr <- rho * Ar
    QSI <- (Arr * z0) + Ar * XS1.g
    eta <- Qv * QSI
    dr = 4 * exp(2 * tau)/((exp(2 * tau) + 1)^2)

    # tau <- XS1.g/r myenv <- new.env() assign('v', v, envir = myenv)
    # #assign('XS0.g',XS0.g,envir = myenv) gv2 <-numericDeriv(quote(pt(-XS0.g,
    # v,log.p=TRUE)), c('v'), myenv)

    f1 <- function(x) {
      ff <- pt(-XS0.g, x, log.p = TRUE)
      return(ff)
    }

    gv2 <- numDeriv::grad(f1, rep(1, length(XS0.g)) * v)
    # sum(gv2) myenv2 <- new.env() assign('v', v, envir = myenv2) assign( 'z0',z0,
    # envir = myenv2) assign( 'QSI',QSI, envir = myenv2) f <-
    # quote(pt((((v+1)/(v+(z0)^2))^(1/2))*QSI,v+1,log.p=TRUE)) gv <- numericDeriv(f,
    # c('v'), myenv2)

    f2 <- function(x) {
      teste = (((x + 1)/(x + (z0)^2))^(1/2)) * QSI
      ff <- pt(teste, (x + 1), log.p = TRUE)
      return(ff)
    }

    gv <- numDeriv::grad(f2, rep(1, length(z0)) * v)

    gradient <- matrix(0, nObs, nParam)
    gradient[YS == 0, istartS] <- -w0 * (XS0) * (dt(-XS0.g, v)/pt(-XS0.g, v))
    gradient[YS == 1, istartS] <- w1 * (XS1) * Ar * Qv * (dt(eta, v + 1)/pt(eta,
                                                                            v + 1))
    gradient[YS == 1, istartO] <- w1 * (XO1) * (Qv/sigma) * ((Qv * z0) + (((QSI *
                                                                              ((v + (z0^2))^(-1))) * z0 - Arr) * (dt(eta, v + 1)/pt(eta, v + 1))))
    gradient[YS == 1, isigma] <- w1 * (-1 + (Qv * z0)^2 + (dt(eta, v + 1)/pt(eta,
                                                                             v + 1)) * ((Qv * (z0)) * (QSI * ((v + (z0^2))^(-1)) * (z0) - Arr)))
    gradient[YS == 1, irho] <- w1 * (dt(eta, v + 1)/pt(eta, v + 1)) * Qv * (Ar^(3)) *
      (z0 + rho * XS1.g) * dr
    # gradient[YS == 1, iv] <- w1 *
    # (((1/2)*v*(digamma((v+1)/2)-digamma(v/2))-(1/2))-((1/2)*v*log(1+((z0^2)/v)))+
    # (((Qv*(z0))^(2))/(2))+v*gv) gradient[YS == 0, iv] <- w0*v*gv2
    # colSums(gradient)
    gradient[YS == 1, iv] <- w1 * ((((1/2) * digamma((v + 1)/2) - (1/2) * digamma(v/2)) -
                                      (1/(2 * v))) - ((1/2) * log(1 + ((z0^2)/v))) + (((Qv * (z0))^(2))/(2 *
                                                                                                           v)) + gv)
    gradient[YS == 0, iv] <- w0 * gv2
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

