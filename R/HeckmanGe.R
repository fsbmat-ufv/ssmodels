#' Function for fit of the Generalized
#' Heckman Model
#'
#' @description
#' Estimates the parameters of the Generalized
#' Heckman model
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
#' fisher_infoHG: Fisher information matrix
#'
#' prop_sigmaHG: Square root of the Fisher information matrix diagonal
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
#' NE: Numerical value that represents the number
#' of parameters related to the covariates fitted
#' to the dispersion parameter considering the
#' constant parameter.
#'
#' NV: Numerical value that represents the number
#' of parameters related to the covariates fitted
#' to the correlation parameter considering the
#' constant parameter.
#'
#' @details
#' The HeckmanGe() function fits a generalization of the Heckman sample
#' selection model, allowing sample selection bias and dispersion parameters
#' to depend on covariates. For more information, see
#' \insertCite{bastosBarreto;textual}{ssmodels}
#'
#' @param selection Selection equation.
#' @param outcome Primary Regression Equation.
#' @param outcomeS Matrix with Covariates for fit of the Dispersion Parameter.
#' @param outcomeC Matrix with Covariates for fit of the Correlation Parameter.
#' @param start initial values.
#' @param data Database.
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' outcomeS <- cbind(age,female,totchr,ins)
#' outcomeC <- 1
#' HeckmanGe(selectEq, outcomeEq,outcomeS, outcomeC, data = MEPS2001)
#' @importFrom Rdpack reprompt
#' @references {
#' \insertAllCited{}
#' }
#' @export HeckmanGe
#' @export
HeckmanGe <- function(selection, outcome, outcomeS, outcomeC, data = sys.frame(sys.parent()), start = NULL) {
    ##############################################################################
    # Extract model matrix and matrix from selection and regression equations
    ##############################################################################
  components <- extract_model_components(selection = selection,
                                         outcome = outcome,
                                         outcomeS = outcomeS,
                                         outcomeC = outcomeC,
                                         data = data)
  XS  <- components$XS
  YS  <- components$YS
  NXS <- components$NXS
  XO  <- components$XO
  YO  <- components$YO
  NXO <- components$NXO
  Msigma <- components$Msigma
  NE <- components$NE
  Mrho <- components$Mrho
  NV <- components$NV
  YSLevels <- components$YSLevels



    ############################ Likelihood #
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

    #######Gradient
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

    ############################## Start #
    if (is.null(start)) {
      start <- step2(YS, XS, YO, XO)
    }

    # Criação dos vetores auxiliares para lambda e kappa
    lambda_start <- rep(1, NE)
    kappa_start  <- rep(0, NV)

    # Parte comum: parâmetros da equação de seleção (g) e de resultado (b)
    g_start <- start[1:NXS]
    b_start <- start[(NXS + 1):(NXS + NXO)]
    sigma_start <- start[NXS + NXO + 1]  # log(sigma) ou intercepto
    rho_start   <- start[NXS + NXO + 2]  # atanh(rho) ou intercepto

    # Composição do vetor final de start
    start <- c(
      g_start,
      b_start,
      if (NE == 1) sigma_start else c(sigma_start, lambda_start),
      if (NV == 1) rho_start   else c(rho_start, kappa_start)
    )

    # Nomeação dos parâmetros
    start_names <- c(
      colnames(XS),
      colnames(XO),
      if (NE == 1) "sigma" else c("interceptS", colnames(outcomeS)),
      if (NV == 1) "correlation" else c("interceptC", colnames(outcomeC))
    )

    names(start) <- start_names


    ###################### optim function
    theta_HG <- optim(start, loglik_gen, gradlik_gen, method = "BFGS", hessian = T,
                      control = list(fnscale = -1))

    theta_HG$par <- postprocess_theta(theta_HG$par, NXS, NXO, NE, NV, XS, XO, outcomeS, outcomeC)



    # Inferência e diagnóstico
    nObs    <- length(YS)
    nParam  <- length(start)
    N0      <- sum(YS == 0)
    N1      <- sum(YS == 1)
    df      <- nObs - nParam
    aic     <- -2 * theta_HG$value + 2 * nParam
    bic     <- -2 * theta_HG$value + nParam * log(nObs)

    # Verifica número de parâmetros nos termos adicionais
    NE <- if (length(outcomeS) == 1) 1 else ncol(model.matrix(~outcomeS))
    NV <- if (length(outcomeC) == 1) 1 else ncol(model.matrix(~outcomeC))

    # Informações derivadas
    fisher_infoHG <- tryCatch(solve(-theta_HG$hessian), error = function(e) matrix(NA, nParam, nParam))
    prop_sigmaHG  <- sqrt(diag(fisher_infoHG))

    # Output final
    result <- list(
      coefficients   = theta_HG$par,
      value          = theta_HG$value,
      loglik         = -theta_HG$value,
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

    class(result) <- c("HeckmanGe", class(theta_HG))
    return(result)

}
