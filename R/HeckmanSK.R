#' Normal Skew Model fit Function
#'
#' Estimates the parameters of the Normal Skew model
#'
#' @param selection Selection equation.
#' @param outcome Primary Regression Equation.
#' @param lambda Initial start for asymmetry parameter.
#' @param start initial values.
#' @param data Database.
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' HeckmanSK(selectEq, outcomeEq, data = MEPS2001, lambda = 1)
#' @export HeckmanSK
#' @export
HeckmanSK <- function(selection, outcome, data = sys.frame(sys.parent()), lambda, start = NULL) {
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
    ########## Regression Matrix #
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

    ########## Likelihood #
    loglik_SK <- function(start) {
        NXS <- dim(model.matrix(~XS))[2] - 1  #Numero de colunas de XS+1
        NXO <- dim(model.matrix(~XO))[2] - 1  #Numero de colunas de XO+1
        ## parameter indices
        istartS <- 1:NXS
        istartO <- seq(tail(istartS, 1) + 1, length = NXO)
        isigma <- tail(istartO, 1) + 1
        irho <- tail(isigma, 1) + 1
        ilamb1 <- tail(irho, 1) + 1
        g <- start[istartS]
        b <- start[istartO]
        sigma <- start[isigma]
        if (sigma < 0)
            return(NA)
        rho <- start[irho]
        if ((rho < -1) || (rho > 1))
            return(NA)
        lamb1 <- start[ilamb1]
        XS.g <- (XS) %*% g
        XO.b <- (XO) %*% b
        u2 <- YO - XO.b
        z <- u2/sigma
        r <- sqrt(1 - rho^2)
        u <- (1 + (lamb1^2) - ((lamb1 * rho)^2))
        lstar <- -(lamb1 * rho)/sqrt(u)
        # lstar <- (-lamb1*rho)/(sqrt(1+(lamb1^2)-((lamb1*rho)^2))) nuc =
        # function(t){2*dnorm(t)*pnorm(-lstar*t)} h=Vectorize(function(ls)
        # integrate(nuc,-Inf, -ls)$value)
        h <- sn::psn(-XS.g, 0, 1, -lstar)
        ll <- ifelse(YS == 0, log(h), log(2/sigma) + log(dnorm(z)) + log(pnorm(lamb1 *
            z)) + log(pnorm((XS.g + rho * z)/r)))
        return(sum(ll))
    }
    ############## Gradient #
    gradlik_SK <- function(start) {
        NXS <- dim(model.matrix(~XS))[2] - 1  #Numero de colunas de XS+1
        NXO <- dim(model.matrix(~XO))[2] - 1  #Numero de colunas de XO+1
        nObs <- length(YS)
        NO <- length(YS[YS > 0])
        nParam <- NXS + NXO + 3  #Total of parameters

        XS0 <- XS[YS == 0, , drop = FALSE]
        XS1 <- XS[YS == 1, , drop = FALSE]
        YO[is.na(YO)] <- 0
        YO1 <- YO[YS == 1]
        XO1 <- XO[YS == 1, , drop = FALSE]
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
        ilamb1 <- tail(irho, 1) + 1

        g <- start[istartS]
        b <- start[istartO]
        sigma <- start[isigma]
        if (sigma < 0)
            return(matrix(NA, nObs, nParam))
        rho <- start[irho]
        if ((rho < -1) || (rho > 1))
            return(matrix(NA, nObs, nParam))
        lamb1 <- start[ilamb1]
        XS0.g <- as.numeric((XS0) %*% g)
        XS1.g <- as.numeric((XS1) %*% g)
        XO1.b <- as.numeric((XO1) %*% b)
        # u2 <- YO1 - XO1.b
        u2 <- YO1 - XO1.b
        z <- u2/sigma
        r <- (1 - rho^2)
        u <- (1 + (lamb1^2) - ((lamb1 * rho)^2))
        lstar <- -(lamb1 * rho)/sqrt(u)
        dlrho = (-rho * (1 + (lamb1^2) - (lamb1 * rho)^2) + lamb1 * rho * (lamb1 -
            lamb1 * (rho^2)))/((1 + (lamb1^2) - (lamb1 * rho)^2))^(3/2)
        dllam = (-rho/sqrt(u)) + lamb1 * rho * (lamb1 - lamb1 * (rho^2))/(sqrt(u) *
            u)
        omeg <- (XS1.g + rho * z)/sqrt(r)
        K1 <- dnorm(omeg)/pnorm(omeg)
        h <- function(t) {
            sn::psn(-t, 0, 1, -lstar)
        }
        h2 <- function(t) {
            sn::psn(-t, 0, 1, lamb1 * rho/sqrt(u))
        }
        K2 <- dnorm(-XS0.g) * pnorm(lstar * XS0.g)/h(XS0.g)
        eta <- lamb1 * z
        K3 <- dnorm(eta)/pnorm(eta)
        K4 <- (dnorm((XS0.g * (sqrt(1 + lamb1^2)))/(sqrt(u))))/h(XS0.g)
        gradient <- matrix(0, nObs, nParam)
        gradient[YS == 0, istartS] <- w0 * (XS0) * (-2 * K2)
        gradient[YS == 1, istartS] <- w1 * (XS1) * (K1/sqrt(r))
        gradient[YS == 1, istartO] <- w1 * (XO1) * ((z/sigma) - ((lamb1/sigma) *
            K3) - ((rho/(sigma * sqrt(r))) * K1))
        gradient[YS == 1, isigma] <- w1 * ((-1/sigma) + ((z^2)/sigma) - ((lamb1 *
            K3 * z)/sigma) - ((rho * K1 * z)/(sigma * sqrt(r))))
        gradient[YS == 0, irho] <- w0 * ((-(2/sqrt(2 * pi)) * (lamb1 * (1 + lamb1^2)/(sqrt(u) *
            (u + (lamb1 * rho)^2))) * dnorm(sqrt(1 + (lamb1 * rho/sqrt(u))^2) *
            XS0.g))/h2(XS0.g))
        gradient[YS == 1, irho] <- w1 * ((K1 * (rho * XS1.g + z))/((sqrt(r))^3))
        gradient[YS == 0, ilamb1] <- w0 * (dllam * sqrt(2/pi) * (1/(1 + lstar^2)) *
            dnorm(sqrt(1 + lstar^2) * (XS0.g))/h(XS0.g))
        gradient[YS == 1, ilamb1] <- w1 * (K3 * z)
        return(colSums(gradient))
    }

    ####### Start#
    if (is.null(start))
        start <- c(step2(YS, XS, YO, XO), lambda)
    #### Optim function#
    theta_SK <- optim(start,
        loglik_SK,
        gradlik_SK,
        method = "BFGS",
        hessian = T,
        control = list(fnscale = -1))
    ########## Results #
    names(theta_SK$par) <- c(colnames(XS), colnames(XO), "sigma", "rho", "lambda")
    a   <- start
    a1  <- theta_SK$par
    a2  <- theta_SK$value
    a3  <- theta_SK$counts[2]
    a4  <- theta_SK$hessian
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
    cl <- class(theta_SK)
    result <- list(coefficients=a1,
        value         =  a2,
        loglik        = -a2,
        counts        =  a3,
        hessian       =  a4,
        fisher_infoSK =  a5,
        prop_sigmaSK  =  a6,
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
    class(result) <- c("HeckmanSK", cl)
    result
}
