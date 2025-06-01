#' Generalized Heckman Model Estimation
#'
#' @description
#' Fits a generalized Heckman sample selection model that allows for heteroskedasticity
#' in the outcome equation and correlation of the error terms depending on covariates.
#' The estimation is performed via Maximum Likelihood using the BFGS algorithm.
#'
#' @details
#' This function extends the classical Heckman selection model by incorporating models
#' for the error term's variance (scale) and the correlation between the selection and outcome equations.
#' The scale model (\code{outcomeS}) allows the error variance of the outcome equation
#' to depend on covariates, while the correlation model (\code{outcomeC}) allows
#' the error correlation to vary with covariates.
#'
#' The optimization is initialized with default or user-supplied starting values,
#' and the results include robust standard errors derived from the inverse of the observed
#' Fisher information matrix.
#'
#' @param selection A formula specifying the selection equation.
#' @param outcome A formula specifying the outcome equation.
#' @param outcomeS A formula or matrix specifying covariates for the scale (variance) model.
#' @param outcomeC A formula or matrix specifying covariates for the correlation model.
#' @param data A data frame containing the variables in the model.
#' @param start An optional numeric vector with starting values for the optimization.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{coefficients}: Named vector of estimated model parameters.
#'   \item \code{value}: Negative of the maximum log-likelihood.
#'   \item \code{loglik}: Maximum log-likelihood.
#'   \item \code{counts}: Number of gradient evaluations performed.
#'   \item \code{hessian}: Hessian matrix at the optimum.
#'   \item \code{fisher_infoHG}: Approximate Fisher information matrix.
#'   \item \code{prop_sigmaHG}: Standard errors for the parameter estimates.
#'   \item \code{level}: Levels of the selection variable.
#'   \item \code{nObs}: Number of observations in the dataset.
#'   \item \code{nParam}: Number of estimated parameters.
#'   \item \code{N0}: Number of censored (unobserved) observations.
#'   \item \code{N1}: Number of uncensored (observed) observations.
#'   \item \code{NXS}: Number of covariates in the selection equation.
#'   \item \code{NXO}: Number of covariates in the outcome equation.
#'   \item \code{df}: Degrees of freedom (observations minus parameters).
#'   \item \code{aic}: Akaike Information Criterion.
#'   \item \code{bic}: Bayesian Information Criterion.
#'   \item \code{initial.value}: Starting values used for optimization.
#'   \item \code{NE}: Number of parameters in the scale model.
#'   \item \code{NV}: Number of parameters in the correlation model.
#' }
#'
#' @examples
#' \dontrun{
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' outcomeS <- ~ educ + income
#' outcomeC <- ~ blhisp + female
#' HeckmanGe(selectEq, outcomeEq, outcomeS = outcomeS, outcomeC = outcomeC, data = MEPS2001)
#' }
#'
#' @references
#' \insertRef{bastosBarreto}{ssmodels}
#'
#' @importFrom Rdpack reprompt
#' @export
HeckmanGe <- function(selection, outcome, outcomeS, outcomeC, data = sys.frame(sys.parent()), start = NULL) {
  ##############################################################################
  # Extract model matrix and response vectors from the selection and outcome equations
  ##############################################################################
  components <- extract_model_components(selection = selection,
                                         outcome = outcome,
                                         outcomeS = outcomeS,
                                         outcomeC = outcomeC,
                                         data = data)
  XS  <- components$XS       # Model matrix for the selection equation
  YS  <- components$YS       # Binary response for the selection equation
  NXS <- components$NXS      # Number of selection covariates
  XO  <- components$XO       # Model matrix for the outcome equation
  YO  <- components$YO       # Response for the outcome equation
  NXO <- components$NXO      # Number of outcome covariates
  Msigma <- components$Msigma  # Covariates for the scale model
  NE <- components$NE          # Number of dispersion parameters
  Mrho <- components$Mrho      # Covariates for the correlation model
  NV <- components$NV          # Number of correlation parameters
  YSLevels <- components$YSLevels  # Levels of the selection variable

  ##############################################################################
  # Log-likelihood function for the generalized Heckman model
  ##############################################################################
  loglik_gen <- function(start) {
    NXS <- ncol(model.matrix(~XS)) - 1
    NXO <- ncol(model.matrix(~XO)) - 1

    NE <- if (length(outcomeS) == 1) 1 else ncol(model.matrix(~outcomeS))
    NV <- if (length(outcomeC) == 1) 1 else ncol(model.matrix(~outcomeC))

    Msigma <- if (NE == 1) matrix(1, nrow(XO), 1) else outcomeS
    Mrho   <- if (NV == 1) matrix(1, nrow(XO), 1) else outcomeC

    ## parameter indices
    istartS <- 1:NXS
    istartO <- (max(istartS) + 1):(max(istartS) + NXO)
    ilambda <- (max(istartO) + 1):(max(istartO) + NE)
    ikappa  <- (max(ilambda) + 1):(max(ilambda) + NV)

    g <- start[istartS]
    b <- start[istartO]
    lambda <- start[ilambda]
    kappa <- start[ikappa]

    mu2 <- XS %*% g
    mu1 <- XO %*% b

    if (NE == 1) {
      sigma <- exp(model.matrix(~1, data = data.frame(rep(1, nrow(XO)))) %*%
                     lambda)
    } else {
      sigma <- exp(model.matrix(~Msigma) %*% lambda)
    }

    if (NV == 1) {
      rho <- tanh(model.matrix(~1, data = data.frame(rep(1, nrow(XO)))) %*%
                    kappa)
    } else {
      rho <- tanh(model.matrix(~Mrho) %*% kappa)
    }



    z1 <- YO - mu1
    z <- z1 * (1/sigma)
    r <- sqrt(1 - rho^2)
    A_rho <- 1/r
    A_rrho <- rho/r
    zeta <- mu2 * A_rho + z * A_rrho
    ll <- ifelse(YS == 0, (pnorm(-mu2, log.p = TRUE)), dnorm(z, log = TRUE) -
                   log(sigma) + (pnorm(zeta, log.p = TRUE)))
    sum(ll)
  }

  ##############################################################################
  # Analytical gradient of the log-likelihood
  ##############################################################################
    gradlik_gen <- function(start) {
      NXS <- ncol(model.matrix(~XS)) - 1
      NXO <- ncol(model.matrix(~XO)) - 1

      NE <- if (length(outcomeS) == 1) 1 else ncol(model.matrix(~outcomeS))
      NV <- if (length(outcomeC) == 1) 1 else ncol(model.matrix(~outcomeC))

      Msigma <- if (NE == 1) matrix(1, nrow(XO), 1) else outcomeS
      Mrho   <- if (NV == 1) matrix(1, nrow(XO), 1) else outcomeC

      nObs <- length(YS)
      nParam <- NXS + NXO + NE + NV

        XS0 <- XS[YS == 0, , drop = FALSE]
        XS1 <- XS[YS == 1, , drop = FALSE]
        YO[is.na(YO)] <- 0
        YO1 <- YO[YS == 1]
        XO1 <- XO[YS == 1, , drop = FALSE]
        ES0 <- Msigma[YS == 0, , drop = FALSE]
        ES1 <- Msigma[YS == 1, , drop = FALSE]
        VS0 <- Mrho[YS == 0, , drop = FALSE]
        VS1 <- Mrho[YS == 1, , drop = FALSE]
        N0 <- sum(YS == 0)
        N1 <- sum(YS == 1)

        M1 <- rep(1, N0 + N1)
        u1 <- rep(1, N1)  #YS=1
        u2 <- rep(1, N0)  #YS=0

        ## parameter indices
        istartS <- 1:NXS
        istartO <- (max(istartS) + 1):(max(istartS) + NXO)
        ilambda <- (max(istartO) + 1):(max(istartO) + NE)
        ikappa  <- (max(ilambda) + 1):(max(ilambda) + NV)


        g <- start[istartS]
        b <- start[istartO]
        lambda <- start[ilambda]
        kappa <- start[ikappa]

        mu20 <- as.numeric(XS0 %*% g)
        mu21 <- as.numeric(XS1 %*% g)
        mu11 <- as.numeric(XO1 %*% b)

        sigma1 <- as.numeric(exp(model.matrix(if (NE == 1) ~ES1 - 1 else ~ES1) %*% lambda))
        rho1   <- as.numeric(tanh(model.matrix(if (NV == 1) ~VS1 - 1 else ~VS1) %*% kappa))

        z <- (YO1 - mu11) / sigma1
        r <- sqrt(1 - rho1^2)
        A_rho <- 1 / r
        A_rrho <- rho1 / r

        zeta <- drop(mu21) * A_rho + z * A_rrho
        MZeta <- exp(dnorm(zeta, log = TRUE) - pnorm(zeta, log.p = TRUE))
        Mmu2 <- exp(dnorm(-mu20, log = TRUE) - pnorm(-mu20, log.p = TRUE))
        gradient <- matrix(0, nObs, nParam)
        gradient[YS == 0, istartS] <- -u2 * XS0 * Mmu2
        gradient[YS == 1, istartS] <- u1 * XS1 * MZeta * A_rho
        gradient[YS == 1, istartO] <- u1 * XO1 * (z - MZeta * A_rrho) / sigma1

        form_sigma <- model.matrix(if (NE == 1) ~ES1 - 1 else ~ES1)
        form_rho   <- model.matrix(if (NV == 1) ~VS1 - 1 else ~VS1)
        sech_sq <- pracma::sech(form_rho %*% kappa)^2

        gradient[YS == 1, ilambda] <- u1 * form_sigma * (z^2 - 1 - MZeta * z * A_rrho)
        gradient[YS == 1, ikappa]  <- u1 * form_rho * MZeta * A_rho * sech_sq *
          (drop(mu21) * rho1 * A_rho^2 + z * (1 + A_rrho^2))
        colSums(gradient)
    }

  ##############################################################################
  # Starting values
  ##############################################################################
  if (is.null(start)) {
    message("Start not provided using default start values.")
    start <- c(rep(0, NXS + NXO), 1, 0)
  }

  # Auxiliary vectors for lambda and kappa parameters
  lambda_start <- rep(1, NE)
  kappa_start  <- rep(0, NV)

  # Starting values for selection (g), outcome (b), scale (sigma) and correlation (rho)
  g_start     <- start[1:NXS]
  b_start     <- start[(NXS + 1):(NXS + NXO)]
  sigma_start <- start[NXS + NXO + 1]  # log(sigma) or intercept
  rho_start   <- start[NXS + NXO + 2]  # atanh(rho) or intercept

  # Combine all initial values
  start <- c(
    g_start,
    b_start,
    if (NE == 1) sigma_start else c(sigma_start, lambda_start),
    if (NV == 1) rho_start   else c(rho_start, kappa_start)
  )

  # Assign names to parameter vector
  start_names <- c(
    colnames(XS),
    colnames(XO),
    if (NE == 1) "sigma" else c("interceptS", colnames(outcomeS)),
    if (NV == 1) "correlation" else c("interceptC", colnames(outcomeC))
  )
  names(start) <- start_names

  ##############################################################################
  # Optimization using BFGS method
  ##############################################################################
  theta_HG <- optim(start, loglik_gen, gradlik_gen, method = "BFGS", hessian = TRUE,
                    control = list(fnscale = -1))

  # Post-process estimated parameters: apply transformations and assign names
  theta_HG$par <- postprocess_theta(theta_HG$par, NXS, NXO, NE, NV, XS, XO, outcomeS, outcomeC)

  ##############################################################################
  # Model diagnostics and information criteria
  ##############################################################################
  nObs    <- length(YS)               # Number of observations
  nParam  <- length(start)            # Total number of parameters
  N0      <- sum(YS == 0)             # Number of censored observations
  N1      <- sum(YS == 1)             # Number of observed (selected) observations
  df      <- nObs - nParam            # Degrees of freedom
  aic     <- -2 * theta_HG$value + 2 * nParam
  bic     <- -2 * theta_HG$value + nParam * log(nObs)

  # Recalculate NE and NV in case of transformation
  NE <- if (length(outcomeS) == 1) 1 else ncol(model.matrix(~outcomeS))
  NV <- if (length(outcomeC) == 1) 1 else ncol(model.matrix(~outcomeC))

  ##############################################################################
  # Fisher information and standard errors
  ##############################################################################
  fisher_infoHG <- tryCatch(solve(-theta_HG$hessian), error = function(e) matrix(NA, nParam, nParam))
  prop_sigmaHG  <- sqrt(diag(fisher_infoHG))

  ##############################################################################
  # Output object
  ##############################################################################
  result <- list(
    coefficients   = theta_HG$par,
    value          = theta_HG$value,
    loglik         = theta_HG$value,
    counts         = theta_HG$counts[2],
    hessian        = theta_HG$hessian,
    fisher_infoHG  = fisher_infoHG,
    prop_sigmaHG   = prop_sigmaHG,
    level          = YSLevels,
    nObs           = nObs,
    nParam         = nParam,
    N0             = N0,
    N1             = N1,
    NXS            = ncol(XS),
    NXO            = ncol(XO),
    df             = df,
    aic            = aic,
    bic            = bic,
    initial.value  = start,
    NE             = NE,
    NV             = NV
  )

  # Assign class to the result object
  class(result) <- c("HeckmanGe", class(theta_HG))
  return(result)
}
