% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MEPS2001}
\alias{MEPS2001}
\title{Medical Expenditure Panel Survey}
\format{
An object of class \code{data.frame} with 3328 rows and 22 columns.
}
\source{
2001 Medical Expenditure Panel Survey by the Agency
for Healthcare Research and Quality.
}
\usage{
MEPS2001
}
\description{
The MEPS is a set of large-scale surveys of families, individuals
and their medical providers (doctors, hospitals, pharmacies, etc.)
in the United States. It has data on the health services Americans use,
how often they use them, the cost of these services and how they are paid,
as well as data on the cost and reach of health insurance available
to American workers. The sample is restricted to persons aged between 21
and 64 years and contains a variable response with 3328 observations of
outpatient costs, of which 526 (15.8\%) correspond to unobserved expenditure
values and identified as zero expenditure for adjustment of the models.
It also includes the following explanatory variables:
\itemize{
  \item{educ: education status}
  \item{age: Age}
  \item{income: income}
  \item{female: gender}
  \item{vgood: a numeric vector}
  \item{good: a numeric vector}
  \item{hospexp: a numeric vector}
  \item{totchr: number of chronic diseases}
  \item{ffs: a numeric vector}
  \item{dhospexp: a numeric vector}
  \item{age2: a numeric vector}
  \item{agefem: a numeric vector}
  \item{fairpoor: a numeric vector}
  \item{year01: a numeric vector}
  \item{instype: a numeric vector}
  \item{ambexp: a numeric vector}
  \item{lambexp: log ambulatory expenditures}
  \item{blhisp: ethnicity}
  \item{instype_s1: a numeric vector}
  \item{dambexp: dummy variable, ambulatory expenditures}
  \item{lnambx: a numeric vector}
  \item{ins: insurance status}
}
}
\examples{
data(MEPS2001)
attach(MEPS2001)
hist(lnambx)
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
HeckmanCL(selectEq, outcomeEq, data = MEPS2001)
}
\references{
{
  \insertRef{colin2009microeconometrics}{ssmodels}

  \insertRef{ssmrob}{ssmodels}

  \insertRef{sampleSelection}{ssmodels}
}
}
\keyword{datasets}
