#' Classic Heckman Model fit Function
#'
#' Estimates the parameters of the classic Heckman model
#' via Maximum Likelihood method. The initial start is obtained
#' via the two-step method.
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
    mfS <- eval(mfS, parent.frame())
    mtS <- terms(mfS)
    XS <- model.matrix(mtS, mfS)
    NXS <- ncol(XS)
    YS <- model.response(mfS)
    YSLevels <- levels(as.factor(YS))
    if (length(YSLevels) != 2) {
        stop("the left hand side of the 'selection' formula\n",
            "has to contain", " exactly two levels (e.g. FALSE and TRUE)")
    }
    ##############################################################################
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
    ##############################################################################
    ####### likelihood #
    loglik_HC <- function(start) {
        NXS <- dim(model.matrix(~XS))[2] - 1
        NXO <- dim(model.matrix(~XO))[2] - 1
        ## parameter indices
        istartS <- 1:NXS
        istartO <- seq(tail(istartS, 1) + 1, length = NXO)
        isigma <- tail(istartO, 1) + 1
        irho <- tail(isigma, 1) + 1
        g <- start[istartS]
        b <- start[istartO]
        sigma <- start[isigma]
        if (sigma < 0)
            return(NA)
        rho <- start[irho]
        if ((rho < -1) || (rho > 1))
            return(NA)
        XS.g <- XS %*% g
        XO.b <- XO %*% b
        u2 <- YO - XO.b
        r <- sqrt(1 - rho^2)
        B <- (XS.g + rho/sigma * u2)/r
        ll <- ifelse(YS == 0, (pnorm(-XS.g, log.p = TRUE)), dnorm(u2/sigma, log = TRUE) -
            log(sigma) + (pnorm(B, log.p = TRUE)))
        return(sum(ll))
    }
    ##############################################################################
    ###### Gradient #
    gradlik_HC <- function(start) {
        NXS  <- dim(model.matrix(~XS))[2] - 1
        NXO  <- dim(model.matrix(~XO))[2] - 1
        nObs <- length(YS)
        NO   <- length(YS[YS > 0])
        nParam <- NXS + NXO + 2

        XS0 <- XS[YS == 0, , drop = FALSE]
        XS1 <- XS[YS == 1, , drop = FALSE]
        YO[is.na(YO)] <- 0
        YO1 <- YO[YS == 1]
        XO1 <- XO[YS == 1, , drop = FALSE]
        N0  <- sum(YS == 0)
        N1  <- sum(YS == 1)

        w  <- rep(1, N0 + N1)
        w0 <- rep(1, N0)
        w1 <- rep(1, N1)

        ## parameter indices
        istartS <- 1:NXS
        istartO <- seq(tail(istartS, 1) + 1, length = NXO)
        isigma  <- tail(istartO, 1) + 1
        irho    <- tail(isigma, 1) + 1

        g <- start[istartS]
        b <- start[istartO]
        sigma <- start[isigma]
        if (sigma < 0)
            return(matrix(NA, nObs, nParam))
        rho <- start[irho]
        if ((rho < -1) || (rho > 1))
            return(matrix(NA, nObs, nParam))
        XS0.g <- as.numeric(XS0 %*% g)
        XS1.g <- as.numeric(XS1 %*% g)
        XO1.b <- as.numeric(XO1 %*% b)
        u2 <- YO1 - XO1.b
        r  <- sqrt(1 - rho^2)
        B  <- (XS1.g + rho/sigma * u2)/r
        lambdaB  <- exp(dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))
        gradient <- matrix(0, nObs, nParam)
        gradient[YS == 0, istartS] <- -w0 * XS0 * exp(dnorm(-XS0.g, log = TRUE) -
            pnorm(-XS0.g, log.p = TRUE))
        gradient[YS == 1, istartS] <- w1 * XS1 * lambdaB/r
        gradient[YS == 1, istartO] <- w1 * XO1 * (u2/sigma^2 - lambdaB * rho/sigma/r)
        gradient[YS == 1, isigma]  <- w1 * ((u2^2/sigma^3 - lambdaB * rho * u2/sigma^2/r) -
            1/sigma)
        gradient[YS == 1, irho]    <- w1 * (lambdaB * (u2/sigma + rho * XS1.g))/r^3
        return(colSums(gradient))
    }

    ####### Start#
    if (is.null(start))
    start <- step2(YS, XS, YO, XO)
    #### Optim function#
    theta_HC <- optim(start,
        loglik_HC,
        gradlik_HC,
        method  = "BFGS",
        hessian = T,
        control = list(fnscale = -1))
########################Results
    names(theta_HC$par) <- c(colnames(XS), colnames(XO), "sigma", "rho")
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
    a12 <- ncol(XS)
    a13 <- ncol(XO)
    a14 <- (a8-a9)
    a15 <- -2*a2 + 2*a9
    a16 <- -2*a2 + a9*log(a8)
    cl <- class(theta_HC)
    result <- list(coefficients=a1,
        value         =  a2,
        loglik        = -a2,
        counts        =  a3,
        hessian       =  a4,
        fisher_infoHC =  a5,
        prop_sigmaHC  =  a6,
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
    class(result) <- c("HeckmanCL", cl)
    result
}
