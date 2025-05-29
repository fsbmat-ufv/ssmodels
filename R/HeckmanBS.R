#' Heckman BS Model fit Function
#'
#' @description
#' Estimates the parameters of the Heckman-BS model
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
#' fisher_infoBS: Fisher information matrix
#'
#' prop_sigmaBS: Square root of the Fisher information matrix diagonal
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
#' The HeckmanBS() function fits the Sample Selection Model
#' based on the Birnbaum Saunders bivariate distribution,
#' it has the same number of parameters as the classical
#' Heckman model. For more information see
#' \insertCite{bastos;textual}{ssmodels}
#'
#'
#' @param selection Selection equation.
#' @param outcome Primary Regression Equation.
#' @param data Database.
#' @param start initial values.
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeBS <- ambexp ~ age + female + educ + blhisp + totchr + ins
#' HeckmanBS(selectEq, outcomeBS, data = MEPS2001)
#' @importFrom Rdpack reprompt
#' @importFrom stats model.frame
#' @references {
#' \insertAllCited{}
#' }
#' @export HeckmanBS
#' @export
HeckmanBS <- function(selection, outcome, data = sys.frame(sys.parent()), start = NULL) {

  components <- extract_model_components(selection = selection, outcome = outcome, data = data)
  XS  <- components$XS
  YS  <- components$YS
  NXS <- components$NXS
  XO  <- components$XO
  YO  <- components$YO
  NXO <- components$NXO

  if (is.null(start)) {
    message("Start not provided using default start values.")
    start <- c(rep(0, NXS + NXO), 1, 0)
  }

  # Auxiliary transformation: maps rho_star ∈ ℝ to rho ∈ (-1, 1)
  inv_logit_2 <- function(x) 2 / (1 + exp(-x)) - 1
  loglik_BS <- function(par) {
    # Check for invalid parameters
    if (any(!is.finite(par))) return(NA)

    # Parameter indexing
    igamma <- 1:NXS
    ibeta  <- seq(tail(igamma, 1) + 1, length.out = NXO)
    iphi1  <- tail(ibeta, 1) + 1
    irho   <- iphi1 + 1

    # Parameter extraction
    gamma <- par[igamma]
    beta  <- par[ibeta]
    phi1  <- par[iphi1]
    if (!is.finite(phi1) || phi1 < 1e-3 || phi1 > 100) return(NA)
    rho_star <- par[irho]
    rho <- inv_logit_2(rho_star)
    phi2 <- 1  # Fixed value for identifiability

    # Prepare data subsets
    XS0 <- XS[YS == 0, , drop = FALSE]
    XS1 <- XS[YS == 1, , drop = FALSE]
    YO_clean <- ifelse(is.na(YO), 0, YO)
    YO1 <- YO_clean[YS == 1]
    XO1 <- XO[YS == 1, , drop = FALSE]

    # Linear predictors
    XS0.g <- exp(XS0 %*% gamma)
    XS1.g <- exp(XS1 %*% gamma)
    XO1.b <- exp(XO1 %*% beta)

    # Compute likelihood terms
    term0 <- sqrt(YO1 * (phi1 + 1)/(phi1 * XO1.b)) - sqrt((phi1 * XO1.b)/(YO1 * (phi1 + 1)))
    term1 <- exp((-phi1 / 4) * term0^2)
    term2 <- sqrt((phi1 + 1)/(phi1 * XO1.b * YO1)) + sqrt((phi1 * XO1.b)/((phi1 + 1) * YO1^3))
    term3 <- (1 / (2 * sqrt(2 * pi))) * sqrt(phi1 / 2)
    temp_rho2 <- pmax(1 - rho^2, 1e-6)
    term4 <- sqrt((phi2 + 1)/(2 * XS1.g * temp_rho2))
    term5 <- (phi2 * XS1.g)/(phi2 + 1) - 1
    term6 <- rho * sqrt(phi1 / (2 * temp_rho2))
    integrand <- term4 * term5 + term6 * term0
    term7 <- pnorm(integrand, log.p = TRUE)

    term8 <- sqrt(phi2 / 2) * (sqrt((phi2 + 1)/(phi2 * XS0.g)) - sqrt((phi2 * XS0.g)/(phi2 + 1)))
    FT2 <- pnorm(term8, log.p = TRUE)

    # Stabilize log computations
    term1[term1 < 1e-300] <- 1e-300
    logterm1 <- log(term1)
    logterm2 <- log(term2)
    if (any(!is.finite(logterm1)) || any(!is.finite(logterm2))) return(NA)

    # Final log-likelihood
    ll <- sum(logterm1 + logterm2 + log(term3) + term7) + sum(FT2)
    return(ll)
  }

  ################################################################################

  gradlik_BS <- function(par) {
    # Parameter indexing
    igamma <- 1:NXS
    ibeta  <- seq(tail(igamma, 1) + 1, length.out = NXO)
    iphi1  <- tail(ibeta, 1) + 1
    irho   <- iphi1 + 1

    # Parameter extraction
    gamma <- par[igamma]
    beta  <- par[ibeta]
    phi1  <- par[iphi1]
    rho_star <- par[irho]
    rho <- inv_logit_2(rho_star)
    phi2 <- 1  # Fixed value

    # Prepare data subsets
    XS0 <- XS[YS == 0, , drop = FALSE]
    XS1 <- XS[YS == 1, , drop = FALSE]
    YO_clean <- ifelse(is.na(YO), 0, YO)
    YO1 <- YO_clean[YS == 1]
    XO1 <- XO[YS == 1, , drop = FALSE]

    # Linear predictors
    XS0.g <- exp(as.numeric(XS0 %*% gamma))
    XS1.g <- exp(as.numeric(XS1 %*% gamma))
    XO1.b <- exp(as.numeric(XO1 %*% beta))

    # Auxiliary terms
    temp_rho2 <- pmax(1 - rho^2, 1e-6)
    term0 <- sqrt(YO1 * (phi1 + 1)/(phi1 * XO1.b)) - sqrt((phi1 * XO1.b)/(YO1 * (phi1 + 1)))
    term1 <- exp((-phi1 / 4) * term0^2)
    term2 <- sqrt((phi1 + 1)/(phi1 * XO1.b * YO1)) + sqrt((phi1 * XO1.b)/((phi1 + 1) * YO1^3))
    term3 <- (1 / (2 * sqrt(2 * pi))) * sqrt(phi1 / 2)
    term4 <- sqrt((phi2 + 1)/(2 * XS1.g * temp_rho2))
    term5 <- (phi2 * XS1.g)/(phi2 + 1) - 1
    term6 <- rho * sqrt(phi1 / (2 * temp_rho2))
    integrand <- term4 * term5 + term6 * term0
    term8 <- sqrt(phi2 / 2) * (sqrt((phi2 + 1)/(phi2 * XS0.g)) - sqrt((phi2 * XS0.g)/(phi2 + 1)))

    # Ratios of PDF to CDF
    lambda_I  <- exp(dnorm(integrand, log = TRUE) - pnorm(integrand, log.p = TRUE))
    lambda_T8 <- exp(dnorm(term8, log = TRUE) - pnorm(term8, log.p = TRUE))

    # Partial derivatives
    term9  <- (-1 / 2) * (sqrt(YO1 * (phi1 + 1)/(phi1 * XO1.b)) + sqrt((phi1 * XO1.b)/(YO1 * (phi1 + 1)))) * model.matrix(~XO1 - 1)
    term10 <- (1 / 2) * (sqrt((XO1.b * phi1)/(YO1^3 * (phi1 + 1))) - sqrt((phi1 + 1)/(YO1 * phi1 * XO1.b))) * model.matrix(~XO1 - 1)
    term11 <- (-1 / 2) * sqrt(phi2 / 2) * (sqrt((phi2 + 1)/(phi2 * XS0.g)) + sqrt((phi2 * XS0.g)/(phi2 + 1))) * model.matrix(~XS0 - 1)
    term12 <- (-1 / 2) * (sqrt(YO1/(phi1^3 * XO1.b * (phi1 + 1))) + sqrt(XO1.b/(YO1 * (phi1 + 1)^3 * phi1)))
    term13 <- (1 / 2) * (sqrt(XO1.b / (phi1 * (phi1 + 1)^3 * YO1^3)) - sqrt(1 / (phi1^3 * (phi1 + 1) * XO1.b * YO1)))
    term14 <- rho / (2 * sqrt(2 * phi1 * temp_rho2))

    # Derivatives w.r.t. parameters
    df1b     <- (-phi1 / 2) * term0 * term9
    df2b     <- term10 / term2
    df4b     <- lambda_I * term9 * term6
    df4g     <- lambda_I * (1/4) * ((XS1.g + 2) / sqrt(XS1.g * temp_rho2)) * model.matrix(~XS1 - 1)
    df5g     <- lambda_T8 * term11
    df1phi1  <- - (term0^2) / 4 - (phi1 / 2) * term0 * term12
    df2phi1  <- term13 / term2
    df3phi1  <- 1 / (8 * term3 * sqrt(pi * phi1))
    df4phi1  <- lambda_I * (term0 * term14 + term6 * term12)

    # Derivative w.r.t. rho using chain rule
    dt5_drho <- if (abs(rho) < 1e-8) sqrt(phi1 / 2) else sqrt(phi1 / (2 * temp_rho2)) + rho^2 * sqrt(phi1) / ((temp_rho2)^(3/2) * sqrt(2))
    d_rho    <- (4 * exp(-rho_star)) / (1 + exp(-rho_star))^2
    df4rho_raw <- 0.5 * lambda_I * (term5 * term4 * (rho / temp_rho2) + dt5_drho * term0)
    gr_rho_star <- df4rho_raw * d_rho

    # Assemble gradient matrix
    gradient <- matrix(0, nrow = length(YS), ncol = length(par))
    gradient[YS == 0, igamma] <- df5g
    gradient[YS == 1, igamma] <- df4g
    gradient[YS == 1, ibeta]  <- df1b + df2b + df4b
    gradient[YS == 1, iphi1]  <- df1phi1 + df2phi1 + df3phi1 + df4phi1
    gradient[YS == 1, irho]   <- gr_rho_star

    return(colSums(gradient))
  }

  theta_BS <- optim(start, loglik_BS, gradlik_BS, method = "BFGS", hessian = TRUE,
                    control = list(fnscale = -1, maxit = 1000))
  names(theta_BS$par) <- c(colnames(XS), colnames(XO), "sigma", "rho")

  irho <- length(theta_BS$par)
  a1 <- theta_BS$par
  a4 <- theta_BS$hessian
  a5 <- tryCatch(solve(-a4), error = function(e) matrix(NA, nrow = length(a1), ncol = length(a1)))
  a6 <- suppressWarnings(sqrt(diag(a5)))
  a6[!is.finite(a6)] <- NA
  rho_star <- a1[irho]
  rho_hat <- 2 / (1 + exp(-rho_star)) - 1
  g_prime <- 4 * exp(-rho_star) / (1 + exp(-rho_star))^2
  se_rho_star <- a6[irho]
  se_rho <- abs(g_prime) * se_rho_star
  a1[irho] <- rho_hat
  a6[irho] <- se_rho

  result <- list(
    coefficients     = a1,
    value            = theta_BS$value,
    loglik           = -theta_BS$value,
    counts           = theta_BS$counts[2],
    hessian          = a4,
    fisher_infoBS    = a5,
    prop_sigmaBS     = a6,
    level            = levels(as.factor(YS)),
    nObs             = length(YS),
    nParam           = length(start),
    N0               = sum(YS == 0),
    N1               = sum(YS == 1),
    NXS              = NXS,
    NXO              = NXO,
    df               = length(YS) - length(start),
    aic              = -2 * theta_BS$value + 2 * length(start),
    bic              = -2 * theta_BS$value + length(start) * log(length(YS)),
    initial.value    = start
  )
  class(result) <- c("HeckmanBS")
  return(result)
}
