## ----setup, include=FALSE-----------------------------------------------------
library(ssmodels)
require(knitr)
require(kfigr)
library(kableExtra)
options(knitr.table.format = "latex")
knitr::opts_chunk$set(echo = TRUE, fig.align = "center", message=FALSE, 
                      warning=FALSE, fig.height=5, fig.width=7.2)

## -----------------------------------------------------------------------------
library("mvtnorm") #Pacote para geração dos dados
set.seed(0) 
#Parâmetros utilizados:
sigma <- 2.9
rho   <- -0.5
n     <- 100
gamma0<- 2.2
gamma1<- -3
beta0 <- -1
beta1 <- 2.7
#Covariáveis
xs    <- runif(n)
xo    <- runif(n)
#Vetor de parâmetros
b     <- rbind(beta0,beta1)
g     <- rbind(gamma0,gamma1)
mu1   <- (as.numeric(model.matrix(~xo) %*% b))
mu2   <- (as.numeric(model.matrix(~xs) %*% g))
#Simulação das variáveis bivariadas com distribuição normal
eps   <- matrix(NA,n,2)
for(k in 1:n){
  eps[k,] <- mvtnorm::rmvnorm(1,
                            c(mu1[k], mu2[k]),
                            matrix(c(sigma, rho*sigma, rho*sigma, 1), 2, 2))
}
#eps <- rmvnorm(1000, c(0,0), matrix(c(sigma, rho*sigma, rho*sigma, 1), 2, 2))
y1 <- cbind(eps[,1])
y2 <- cbind(eps[,2])>0
#Variável de seleção
ys<-1*y2
#Variável de interesse primário
yo<-y1*ys

## ----message=FALSE------------------------------------------------------------
library(ssmodels)
m0 <- HeckmanCL(ys~xs,yo~xo)
summary(m0)
library("sampleSelection")
m1 <- selection(ys~xs, yo ~xo, method = "ml")
summary(m1)

## ----fig1, fig.align= "center", fig.width = 7.2, fig.height= 5, anchor="Figura"----
opar <- par(mar=c(2,2,2,0) + 0.9, mgp=c(2,1,0))
pch <- c(1, 16)
plot(xo, y1, pch=pch[1 + ys], cex=0.5, lwd=0.5, main = "Figura 1: Ajuste dos modelos de Heckman Clássico e modelo linear \n simples a dados simulados com censura.")
# True dependence
abline(a=beta0, b=beta1, lty=1, lwd=2)
# Heckman's model
abline(a=coef(m0)[3], b=coef(m0)[4], lty=2, col="red", lwd=2)
# linear model
cf <- coef(lm(yo ~ xo, subset=ys==1))
abline(a=cf[1], b=cf[2], lty=3, col="blue", lwd=3)
par(opar)

## -----------------------------------------------------------------------------
library("mvtnorm")
set.seed(0)
gamma0 <- 0.5
gamma1 <- -1
beta0 <- -2
beta1 <- 2
lambda0 <- 0.1
lambda1<- 0.5
kappa0 <- 0.3
kappa1 <- -0.5
a <- rbind(lambda0,lambda1)
b <- rbind(beta0,beta1)
g <- rbind(gamma0,gamma1)
d <- rbind(kappa0,kappa1)
n <- 200
xs <- rnorm(n)
xo <- rnorm(n)
mu1 <- (as.numeric(model.matrix(~xo) %*% b))
mu2 <- (as.numeric(model.matrix(~xs) %*% g))
sigma <- exp(as.numeric(model.matrix(~xo) %*% a))
rho   <- tanh(as.numeric(model.matrix(~xo) %*% d))
eps <- matrix(NA,n,2)
for(k in 1:n){
  eps[k,] <- mvtnorm::rmvnorm(1,c(mu1[k],mu2[k]),matrix(c((sigma[k])^2,rho[k]*sigma[k],rho[k]*sigma[k],1),2,2))
}
y1 <- cbind(eps[,1])
y2 <- cbind(eps[,2])>0
ys<-1*y2
yo<-y1*ys

