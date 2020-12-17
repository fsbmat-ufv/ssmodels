#' Heckman BS Model fit Function
#'
#' Estimates the parameters of the Heckman-BS model
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
#' @export HeckmanBS
#' @export
HeckmanBS <- function(selection, outcome, data = sys.frame(sys.parent()), start = NULL) {
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
    #### Regression Matrix
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
    ##############################Start#
    if (is.null(start))
    start <- step2(YS, XS, log(YO), XO)
    ##### Likelihood#
    loglik_BS <- function(par) {
        n <- length(YO)
        NXS <- dim(model.matrix(~XS))[2] - 1
        NXO <- dim(model.matrix(~XO))[2] - 1
        ## parameter indices
        igamma <- 1:NXS
        ibeta <- seq(tail(igamma, 1) + 1, length = NXO)
        iphi1 <- tail(ibeta, 1) + 1
        irho <- tail(iphi1, 1) + 1
        gamma <- par[igamma]
        beta <- par[ibeta]
        phi1 <- par[iphi1]
        if (phi1 < 0)
            return(NA)
        rho <- par[irho]
        if ((rho < -1) || (rho > 1))
            return(NA)
        phi2 <- 1
        XS0 <- XS[YS == 0, , drop = FALSE]
        XS1 <- XS[YS == 1, , drop = FALSE]
        YO[is.na(YO)] <- 0
        YO1 <- YO[YS == 1]
        XO1 <- XO[YS == 1, , drop = FALSE]
        N0 <- sum(YS == 0)
        N1 <- sum(YS == 1)

        XS0.g <- exp(as.numeric((XS0) %*% gamma))
        XS1.g <- exp(as.numeric((XS1) %*% gamma))
        XO1.b <- exp(as.numeric((XO1) %*% beta))

        term0 <- ((YO1 * (phi1 + 1)/(phi1 * XO1.b))^(1/2) - ((phi1 * XO1.b)/(YO1 *
            (phi1 + 1)))^(1/2))
        term1 <- exp((-phi1/4) * (term0)^2)
        term2 <- (((phi1 + 1)/(phi1 * XO1.b * YO1))^(1/2) + ((phi1 * XO1.b)/((phi1 +
            1) * (YO1^3)))^(1/2))
        term3 <- (1/(2 * sqrt(2 * pi))) * ((phi1/2)^(1/2))
        term4 <- (((phi2 + 1))/(2 * XS1.g * (1 - rho^2)))^(1/2)
        term5 <- ((phi2 * XS1.g)/(phi2 + 1)) - 1
        term6 <- rho * (phi1/(2 * (1 - rho^2)))^(1/2)
        integrand <- term4 * term5 + term6 * term0
        term7 <- pnorm(integrand, log.p = TRUE)
        term8 <- ((phi2/2)^(1/2)) * (((phi2 + 1)/(phi2 * XS0.g))^(1/2) - ((phi2 *
            XS0.g)/(phi2 + 1))^(1/2))
        FT2 <- pnorm(term8, log.p = TRUE)
        ll <- sum(log(term2) + log(term1) + log(term3) + term7) + sum(FT2)
        return(sum(ll))
    }
    ### Gradient #
    gradlik_BS <- function(par) {
        n <- length(YO)
        NXS <- dim(model.matrix(~XS))[2] - 1
        NXO <- dim(model.matrix(~XO))[2] - 1
        ## parameter indices
        igamma <- 1:NXS
        ibeta <- seq(tail(igamma, 1) + 1, length = NXO)
        iphi1 <- tail(ibeta, 1) + 1
        irho <- tail(iphi1, 1) + 1
        gamma <- par[igamma]
        beta <- par[ibeta]
        phi1 <- par[iphi1]

        nObs <- length(YS)
        NO <- length(YS[YS > 0])
        nParam <- NXS + NXO + 2

        # if(phi1 < 0) return(matrix(NA, nObs, nParam))
        rho <- par[irho]
        # if( ( rho < -1) || ( rho > 1)) return(matrix(NA, nObs, nParam))
        phi2 <- 1
        XS0 <- XS[YS == 0, , drop = FALSE]
        XS1 <- XS[YS == 1, , drop = FALSE]
        YO[is.na(YO)] <- 0
        YO1 <- YO[YS == 1]
        XO1 <- XO[YS == 1, , drop = FALSE]
        N0 <- sum(YS == 0)
        N1 <- sum(YS == 1)

        XS0.g <- exp(as.numeric((XS0) %*% gamma))
        XS1.g <- exp(as.numeric((XS1) %*% gamma))
        XO1.b <- exp(as.numeric((XO1) %*% beta))

        w <- rep(1, N0 + N1)
        w0 <- rep(1, N0)
        w1 <- rep(1, N1)

        mu1 <- exp(as.numeric((XO) %*% beta))
        mu2 <- exp(as.numeric((XS) %*% gamma))

        term0 <- ((YO1 * (phi1 + 1)/(phi1 * XO1.b))^(1/2) - ((phi1 * XO1.b)/(YO1 *
            (phi1 + 1)))^(1/2))
        term1 <- exp((-phi1/4) * (term0)^2)
        term2 <- (((phi1 + 1)/(phi1 * XO1.b * YO1))^(1/2) + ((phi1 * XO1.b)/((phi1 +
            1) * (YO1^3)))^(1/2))
        term3 <- (1/(2 * sqrt(2 * pi))) * ((phi1/2)^(1/2))
        term4 <- (((phi2 + 1))/(2 * XS1.g * (1 - rho^2)))^(1/2)
        term5 <- ((phi2 * XS1.g)/(phi2 + 1)) - 1
        term6 <- rho * (phi1/(2 * (1 - rho^2)))^(1/2)
        integrand <- term4 * term5 + term6 * term0
        term7 <- pnorm(integrand, log.p = TRUE)
        term8 <- ((phi2/2)^(1/2)) * (((phi2 + 1)/(phi2 * XS0.g))^(1/2) - ((phi2 *
            XS0.g)/(phi2 + 1))^(1/2))
        FT2 <- pnorm(term8, log.p = TRUE)
        term9 <- ((-1/2) * (((YO1 * (phi1 + 1))/(phi1 * (XO1.b)))^(1/2) + ((phi1 *
            XO1.b)/(YO1 * (phi1 + 1)))^(1/2))) * (XO1)  #Derivada de term0 em relação a beta
        term10 <- (1/2) * (((XO1.b * phi1)/((YO1^3) * (phi1 + 1)))^(1/2) - ((phi1 +
            1)/(YO1 * phi1 * XO1.b))^(1/2)) * (XO1)  #Derivada de term2 em relação a beta
        term11 <- ((-1/2) * (sqrt(phi2/2)) * (((phi2 + 1)/(phi2 * (XS0.g)))^(1/2) +
            ((phi2 * XS0.g)/((phi2 + 1)))^(1/2))) * (XS0)  #Derivada de term8 em relação a gamma
        term12 <- (-1/2) * ((YO1/((phi1^3) * XO1.b * (phi1 + 1)))^(1/2) + (XO1.b/(YO1 *
            ((phi1 + 1)^3) * phi1))^(1/2))  #Derivada de term0 em relação a phi1
        term13 <- (1/2) * ((XO1.b/(phi1 * ((phi1 + 1)^3) * (YO1^3)))^(1/2) - (1/((phi1^3) *
            (phi1 + 1) * XO1.b * YO1))^(1/2))  #Derivada de term2 em relação a phi1
        term14 <- rho/(2 * ((2 * phi1 * (1 - rho^2))^(1/2)))  #Derivada de term6 em relação a phi1
        term15 <- (rho/(1 - rho^2)) * term4  #Derivada de term4 em relação a rho
        term16 <- (1/(1 - rho^2)) * ((phi1/(2 * (1 - rho^2)))^(1/2))  #Derivada de term6 em relação a rho
        # f1=log(term1), f2=log(term2), f3=log(term3), f4=term7, f5=FT2
        lambda_I <- exp(dnorm(integrand, log = TRUE) - pnorm(integrand, log.p = TRUE))
        lambda_T8 <- exp(dnorm(term8, log = TRUE) - pnorm(term8, log.p = TRUE))
        df1b <- (-phi1/2) * term0 * term9  #Derivada de log(T1) em relação a beta
        df2b <- (term10/term2)  #Derivada de log(T2) em relação a beta0
        df4b <- lambda_I * term9 * term6  #Derivada de log(T7) em relação a beta

        df4g <- lambda_I * ((1/4) * (((XS1.g + 2)/sqrt(XS1.g * (1 - rho^(2)))))) *
            (XS1)  #Derivada de f4 em relação a gamma

        df5g <- lambda_T8 * term11  #Derivada de FT2 em relação a gamma
        df1phi1 <- (-(term0^2)/4) - (phi1/2) * (term0 * term12)  #Derivada do log(T1) em relação a phi1
        df2phi1 <- (1/term2) * term13
        df3phi1 <- 1/(8 * term3 * ((pi * phi1)^(1/2)))
        df4phi1 <- lambda_I * (term0 * term14 + term6 * term12)

        df4rho <- lambda_I * (term5 * term4 * (rho/(1 - rho^2)) + (term0 * term6/(rho *
            (1 - rho^2))))
        gradient <- matrix(0, nObs, nParam)
        gradient[YS == 0, igamma] <- w0 * df5g
        gradient[YS == 1, igamma] <- w1 * (df4g)
        gradient[YS == 1, ibeta] <- w1 * (df1b + df2b + df4b)
        gradient[YS == 1, iphi1] <- w1 * (df1phi1 + df2phi1 + df3phi1 + df4phi1)
        gradient[YS == 1, irho] <- w1 * (df4rho)
        return(colSums(gradient))
    }
    ####### Optim function #
    theta_BS <- optim(start,
        loglik_BS,
        gradlik_BS,
        method = "BFGS",
        hessian = T,
        control = list(fnscale = -1))

    ############# Results #
    names(theta_BS$par) <- c(colnames(XS), colnames(XO), "sigma", "rho")
    a   <- start
    a1  <- theta_BS$par
    a2  <- theta_BS$value
    a3  <- theta_BS$counts[2]
    a4  <- theta_BS$hessian
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

    cl <- class(theta_BS)
    result <- list(coefficients=a1,
        value         =  a2,
        loglik        = -a2,
        counts        =  a3,
        hessian       =  a4,
        fisher_infoBS =  a5,
        prop_sigmaBS  =  a6,
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
    class(result) <- c("HeckmanBS", cl)
    result
}
