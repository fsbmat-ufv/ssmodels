% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeckmanGe.R
\name{HeckmanGe}
\alias{HeckmanGe}
\title{Function for fit of the Generalized Heckman Model}
\usage{
HeckmanGe(
  selection,
  outcome,
  outcomeS,
  outcomeC,
  data = sys.frame(sys.parent()),
  start = NULL
)
}
\arguments{
\item{selection}{Selection equation.}

\item{outcome}{Primary Regression Equation.}

\item{outcomeS}{Matrix with Covariables for fit of the Dispersion Parameter.}

\item{outcomeC}{Matrix with Covariates for Adjusting the Correlation Parameter.}

\item{data}{Database.}

\item{start}{initial values.}
}
\description{
Estimates the parameters of the Generalized Heckman model
}
\examples{
data(MEPS2001)
attach(MEPS2001)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
outcomeS <- cbind(age,female,totchr,ins)
outcomeC <- 1
HeckmanGe(selectEq, outcomeEq,outcomeS, outcomeC, data = MEPS2001)
}
