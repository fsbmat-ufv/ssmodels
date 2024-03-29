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
    # model.frame requires the parameter to be 'formula'
    mfS <- eval(mfS, parent.frame())
    mtS <- terms(mfS)
    XS <- model.matrix(mtS, mfS)
    NXS <- ncol(XS)
    YS <- model.response(mfS)
    YSLevels <- levels(as.factor(YS))
    ################# Regression Matrix #
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

    ####### Start#
    if (is.null(start))
    start <- c(step2(YS, XS, YO, XO), df)
    #### Optim function#
    theta_tS <- optim(start, loglik_tS, gradlik_tS, method = "BFGS", hessian = T,
        control = list(fnscale = -1))
    ########### Results #
    names(theta_tS$par) <- c(colnames(XS), colnames(XO), "sigma", "rho", "df")
    a   <- start
    a1  <- theta_tS$par
    a2  <- theta_tS$value
    a3  <- theta_tS$counts[2]
    a4  <- theta_tS$hessian
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
    cl <- class(theta_tS)
    result <- list(coefficients=a1,
        value         =  a2,
        loglik        = -a2,
        counts        =  a3,
        hessian       =  a4,
        fisher_infotS =  a5,
        prop_sigmatS  =  a6,
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
        initial.value = a)
    class(result) <- c("HeckmantS", cl)
    result
}
