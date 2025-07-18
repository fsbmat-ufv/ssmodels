% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeckmanGe.R
\name{HeckmanGe}
\alias{HeckmanGe}
\title{Generalized Heckman Model Estimation}
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
\item{selection}{A formula specifying the selection equation.}

\item{outcome}{A formula specifying the outcome equation.}

\item{outcomeS}{A formula or matrix specifying covariates for the scale (variance) model.}

\item{outcomeC}{A formula or matrix specifying covariates for the correlation model.}

\item{data}{A data frame containing the variables in the model.}

\item{start}{An optional numeric vector with starting values for the optimization.}
}
\value{
A list containing:
\itemize{
  \item \code{coefficients}: Named vector of estimated model parameters.
  \item \code{value}: Negative of the maximum log-likelihood.
  \item \code{loglik}: Maximum log-likelihood.
  \item \code{counts}: Number of gradient evaluations performed.
  \item \code{hessian}: Hessian matrix at the optimum.
  \item \code{fisher_infoHG}: Approximate Fisher information matrix.
  \item \code{prop_sigmaHG}: Standard errors for the parameter estimates.
  \item \code{level}: Levels of the selection variable.
  \item \code{nObs}: Number of observations in the dataset.
  \item \code{nParam}: Number of estimated parameters.
  \item \code{N0}: Number of censored (unobserved) observations.
  \item \code{N1}: Number of uncensored (observed) observations.
  \item \code{NXS}: Number of covariates in the selection equation.
  \item \code{NXO}: Number of covariates in the outcome equation.
  \item \code{df}: Degrees of freedom (observations minus parameters).
  \item \code{aic}: Akaike Information Criterion.
  \item \code{bic}: Bayesian Information Criterion.
  \item \code{initial.value}: Starting values used for optimization.
  \item \code{NE}: Number of parameters in the scale model.
  \item \code{NV}: Number of parameters in the correlation model.
}
}
\description{
Fits a generalized Heckman sample selection model that allows for heteroskedasticity
in the outcome equation and correlation of the error terms depending on covariates.
The estimation is performed via Maximum Likelihood using the BFGS algorithm.
}
\details{
This function extends the classical Heckman selection model by incorporating models
for the error term's variance (scale) and the correlation between the selection and outcome equations.
The scale model (\code{outcomeS}) allows the error variance of the outcome equation
to depend on covariates, while the correlation model (\code{outcomeC}) allows
the error correlation to vary with covariates.

The optimization is initialized with default or user-supplied starting values,
and the results include robust standard errors derived from the inverse of the observed
Fisher information matrix.
}
\examples{
\dontrun{
data(MEPS2001)
attach(MEPS2001)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
outcomeS <- ~ educ + income
outcomeC <- ~ blhisp + female
HeckmanGe(selectEq, outcomeEq, outcomeS = outcomeS, outcomeC = outcomeC, data = MEPS2001)
}

}
\references{
\insertRef{bastosBarreto}{ssmodels}
}
