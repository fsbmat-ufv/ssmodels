#' Function for fit of the Generalized
#' Heckman Model
#'
#' @description
#' Estimates the parameters of the Generalized
#' Heckman model
#'
#' #' @return
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
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("selection", "data", "subset"), names(mf), 0)
    mfS <- mf[c(1, m)]
    mfS$drop.unused.levels <- TRUE
    mfS$na.action <- na.pass
    mfS[[1]] <- as.name("model.frame")
    names(mfS)[2] <- "formula"
    # model.frame requires the parameter to be formula
    mfS <- eval(mfS, parent.frame())
    mtS <- terms(mfS)
    XS <- model.matrix(mtS, mfS)
    NXS <- ncol(XS)
    YS <- model.response(mfS)
    YSLevels <- levels(as.factor(YS))
    ############################## Regression Matrix #
    m <- match(c("outcome", "data", "subset", "weights", "offset"), names(mf), 0)
    mfO <- mf[c(1, m)]
    mfO$na.action <- na.pass
    mfO$drop.unused.levels <- TRUE
    mfO$na.action <- na.pass
    mfO[[1]] <- as.name("model.frame")
    names(mfO)[2] <- "formula"
    mfO <- eval(mfO, parent.frame())
    mtO <- attr(mfO, "terms")
    XO <- model.matrix(mtO, mfO)
    NXO <- ncol(XO)
    YO <- model.response(mfO)

    #################### Dispersion Matrix #
    E <- outcomeS
    # NE <- ncol(E)
    if (length(E) == 1) {
        NE <- 1
        Msigma <- cbind(rep(1, nrow(XO)))
    } else {
        NE <- dim(model.matrix(~E))[2] - 1
        Msigma <- E
    }
    ################## Correlation Matrix #
    V <- outcomeC
    if (length(V) == 1) {
        NV <- 1
        Mrho <- cbind(rep(1, nrow(XO)))
    } else {
        NV <- dim(model.matrix(~V))[2] - 1
        Mrho <- V
    }

    ############################ Likelihood #
    loglik_gen <- function(start) {
        NXS <- dim(model.matrix(~XS))[2] - 1
        NXO <- dim(model.matrix(~XO))[2] - 1

        if (length(E) == 1) {
            NE <- 1
            Msigma <- cbind(rep(1, nrow(XO)))
        } else {
            NE <- dim(model.matrix(~E))[2]
            Msigma <- E
        }
        if (length(V) == 1) {
            NV <- 1
            Mrho <- cbind(rep(1, nrow(XO)))
        } else {
            NV <- dim(model.matrix(~V))[2]
            Mrho <- V
        }

        ## parameter indices
        istartS <- 1:NXS
        istartO <- seq(tail(istartS, 1) + 1, length = NXO)
        ilambda <- seq(tail(istartO, 1) + 1, length = NE)
        ikappa <- seq(tail(ilambda, 1) + 1, length = NV)

        g <- start[istartS]
        b <- start[istartO]
        lambda <- start[ilambda]
        # if(sigma < 0) return(NA)
        kappa <- start[ikappa]
        # if( ( rho < -1) || ( rho > 1)) return(NA)

        mu2 <- XS %*% g
        mu1 <- XO %*% b
        # sigma <- exp(model.matrix(~E) %*% lambda) rho <- tanh(model.matrix(~1,data =
        # data.frame(rep(1,nrow(XO))))%*% kappa)

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
        NXS <- dim(model.matrix(~XS))[2] - 1
        NXO <- dim(model.matrix(~XO))[2] - 1

        if (length(E) == 1) {
            NE <- 1
            Msigma <- cbind(rep(1, nrow(XO)))
        } else {
            NE <- dim(model.matrix(~E))[2]
            Msigma <- E
        }
        if (length(V) == 1) {
            NV <- 1
            Mrho <- cbind(rep(1, nrow(XO)))
        } else {
            NV <- dim(model.matrix(~V))[2]
            Mrho <- V
        }
        nObs <- length(YS)
        NO <- length(YS[YS > 0])
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
        istartO <- seq(tail(istartS, 1) + 1, length = NXO)
        ilambda <- seq(tail(istartO, 1) + 1, length = NE)
        ikappa <- seq(tail(ilambda, 1) + 1, length = NV)

        # if(sigma < 0) return(NA) ikappa <- tail(ilambda, 1) + 1

        g <- start[istartS]
        b <- start[istartO]
        lambda <- start[ilambda]
        # if(sigma < 0) return(matrix(NA, nObs, nParam))
        kappa <- start[ikappa]
        # if( ( rho < -1) || ( rho > 1)) return(matrix(NA, nObs, nParam))
        mu20 <- as.numeric(XS0 %*% g)
        mu21 <- as.numeric(XS1 %*% g)
        mu11 <- as.numeric(XO1 %*% b)
        # sigma0 <- exp(as.numeric(model.matrix(~ES0) %*% lambda)) sigma1 <-
        # exp(as.numeric(model.matrix(~ES1) %*% lambda)) rho0 <-
        # tanh(as.numeric(model.matrix(~VS0) %*% kappa)) rho1 <-
        # tanh(as.numeric(model.matrix(~VS1) %*% kappa))
        if (NE == 1) {
            sigma0 <- exp(as.numeric(model.matrix(~ES0 - 1) %*% lambda))
            sigma1 <- exp(as.numeric(model.matrix(~ES1 - 1) %*% lambda))
        } else {
            sigma0 <- exp(as.numeric(model.matrix(~ES0) %*% lambda))
            sigma1 <- exp(as.numeric(model.matrix(~ES1) %*% lambda))
        }

        if (NV == 1) {
            rho0 <- tanh(as.numeric(model.matrix(~VS0 - 1) %*% kappa))
            rho1 <- tanh(as.numeric(model.matrix(~VS1 - 1) %*% kappa))
        } else {
            rho0 <- tanh(as.numeric(model.matrix(~VS0) %*% kappa))
            rho1 <- tanh(as.numeric(model.matrix(~VS1) %*% kappa))
        }
        z1 <- YO1 - mu11
        z <- z1/sigma1
        r <- sqrt(1 - rho1^2)
        A_rho <- 1/r
        A_rrho <- rho1/r
        # B <- (XS1.g + rho/sigma*u2)/r
        zeta <- (mu21 * A_rho + z * A_rrho)
        MZeta <- exp(dnorm(zeta, log = TRUE) - pnorm(zeta, log.p = TRUE))
        Mmu2 <- exp(dnorm(-mu20, log = TRUE) - pnorm(-mu20, log.p = TRUE))
        Q_rho <- mu21 * rho1 * ((A_rho)^2) + z * (1 + A_rrho^2)
        Q_rrho <- mu21 * (1 + 2 * A_rrho^2) + 2 * z * rho1 * (1 + A_rrho^2)

        gradient <- matrix(0, nObs, nParam)
        gradient[YS == 0, istartS] <- -u2 * XS0 * Mmu2
        gradient[YS == 1, istartS] <- u1 * XS1 * MZeta * A_rho
        gradient[YS == 1, istartO] <- u1 * XO1 * (z - MZeta * A_rrho) * (1/sigma1)
        if (NE == 1) {
            gradient[YS == 1, ilambda] <- u1 * model.matrix(~ES1 - 1) * (z^2 - 1 -
                                                                             MZeta * z * A_rrho)
        } else {
            gradient[YS == 1, ilambda] <- u1 * model.matrix(~ES1) * (z^2 - 1 - MZeta *
                                                                         z * A_rrho)
        }
        if (NV == 1) {
            gradient[YS == 1, ikappa] <- u1 * model.matrix(~VS1 - 1) * MZeta * A_rho *
                ((pracma::sech(as.numeric(model.matrix(~VS1 - 1) %*% kappa)))^2) *
                (mu21 * rho1 * (A_rho^2) + z * (1 + (A_rrho^2)))
        } else {
            gradient[YS == 1, ikappa] <- u1 * model.matrix(~VS1) * MZeta * A_rho *
                ((pracma::sech(as.numeric(model.matrix(~VS1) %*% kappa)))^2) * (mu21 *
                                                                                    rho1 * (A_rho^2) + z * (1 + (A_rrho^2)))
        }
        colSums(gradient)
    }

    ##############################Start#
    if (is.null(start))
        start <- step2(YS, XS, YO, XO)
    ####################################
    ilambda <- rep(1, NE)
    ikappa <- rep(0, NV)
    if (length(V) == 1 & length(E) == 1) {
        start <- c(start[(1:NXS)], start[((NXS + 1):(NXS + NXO))], start[(NXS + NXO +
                                                                              1)], start[(NXS + NXO + 2)])
        names(start) <- c(colnames(XS), colnames(XO), "sigma", colnames(E), "correlation")
    } else {
        if (length(V) == 1 & length(E) != 1) {
            start <- c(start[(1:NXS)], start[((NXS + 1):(NXS + NXO))], start[(NXS +
                                                                                  NXO + 1)], ilambda, start[(NXS + NXO + 2)])
            names(start) <- c(colnames(XS), colnames(XO), "interceptS", colnames(E),
                              "correlation")
        } else {
            if (length(V) != 1 & length(E) == 1) {
                start <- c(start[(1:NXS)], start[((NXS + 1):(NXS + NXO))], start[(NXS +
                                                                                      NXO + 1)], start[(NXS + NXO + 2)], ikappa)
                names(start) <- c(colnames(XS), colnames(XO), "sigma", "interceptC",
                                  colnames(V))
            } else {
                start <- c(start[(1:NXS)], start[((NXS + 1):(NXS + NXO))], start[(NXS +
                                                                                      NXO + 1)], ilambda, start[(NXS + NXO + 2)], ikappa)
                names(start) <- c(colnames(XS), colnames(XO), "interceptS", colnames(E),
                                  "interceptC", colnames(V))
            }
        }
    }

    ###################### optim function
    theta_HG <- optim(start, loglik_gen, gradlik_gen, method = "BFGS", hessian = T,
                      control = list(fnscale = -1))

    ############################Results#
    if (length(E) == 1 & length(V) == 1) {
        names(theta_HG$par) <- c(colnames(XS), colnames(XO), "Sigma", colnames(E),
                                 "Rho", colnames(V))
    } else {
        if (length(E) != 1 & length(V) == 1) {
            names(theta_HG$par) <- c(colnames(XS), colnames(XO), "interceptS", colnames(E),
                                     "rho")
        } else {
            if (length(E) == 1 & length(V) != 1) {
                names(theta_HG$par) <- c(colnames(XS), colnames(XO), "Sigma", "interceptC",
                                         colnames(V))
            } else {
                names(theta_HG$par) <- c(colnames(XS), colnames(XO), "interceptS",
                                         colnames(E), "interceptC", colnames(V))
            }
        }
    }


    if (length(E) == 1 & length(V) == 1) {
        theta_HG$par <- c(theta_HG$par[1:(NXS+NXO)], exp(theta_HG$par[NXS+NXO+1]),tanh(theta_HG$par[NXS+NXO+2]))
    } else {
        if (length(E) != 1 & length(V) == 1) {
            theta_HG$par <- c(theta_HG$par[1:(NXS+NXO)], theta_HG$par[(NXS+NXO+1):(NXS + NXO + NE+1)],tanh(theta_HG$par[NXS + NXO + NE+2]))
        } else {
            if (length(E) == 1 & length(V) != 1) {
                theta_HG$par <- c(theta_HG$par[1:(NXS+NXO)], exp(theta_HG$par[(NXS+NXO+NE)]), theta_HG$par[(NXS + NXO + NE+1):(NXS + NXO + NE + NV + 1)])
            } else {
                theta_HG$par <- theta_HG$par
            }
        }
    }


    a   <- start
    a1  <- theta_HG$par
    a2  <- theta_HG$value
    a3  <- theta_HG$counts[2]
    a4  <- theta_HG$hessian
    a5  <- solve(-a4)
    a6  <- sqrt(diag(a5))
    a7  <- YSLevels
    a8  <- length(YS)
    a9  <- length(start)
    a10 <- sum(YS == 0)
    a11 <- sum(YS == 1)
    a12 <- ncol(XS)
    a13 <- ncol(XO)
    a14 <- (a8-a9)
    a15 <- -2*a2 + 2*a9
    a16 <- -2*a2 + a9*log(a8)
    if (length(E) == 1) {
        NE <- 1
    } else {
        NE <- dim(model.matrix(~E))[2]
    }
    if (length(V) == 1) {
        NV <- 1
    } else {
        NV <- dim(model.matrix(~V))[2]
    }

    cl <- class(theta_HG)
    result <- list(coefficients=a1,
                   value         =  a2,
                   loglik        = -a2,
                   counts        =  a3,
                   hessian       =  a4,
                   fisher_infoHG =  a5,
                   prop_sigmaHG  =  a6,
                   level         =  a7,
                   nObs          =  a8,
                   nParam        =  a9,
                   N0            = a10,
                   N1            = a11,
                   NXS           = a12,
                   NXO           = a13,
                   df            = a14,
                   aic           = a15,
                   bic           = a16,
                   initial.value = a,
                   NE            = NE,
                   NV            = NV)
    class(result) <- c("HeckmanGe", cl)
    result
}
