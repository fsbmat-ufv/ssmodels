#' ssmodels: A package for fit the sample selection models.
#'
#' Package that provides models to fit data with sample selection bias problems. Includes:
#' \describe{
#' \item{HeckmanCL(selectEq, outcomeEq, data = data, start)}{Heckman's classic model fit function. Sample selection
#' usually arises in practice as a result of partial observability of the
#' result of interest in a study. In the presence of sample selection, the
#' observed data do not represent a random sample of the population, even
#' after controlling for explanatory variables. That is, data is missing
#' randomly. Thus, standard analysis using only complete cases will lead to
#' biased results. Heckman introduced a sample selection model to analyze
#' this data and proposed a complete likelihood estimation method under the
#' assumption of normality. Such model was called Heckman model or Tobit 2
#' model.}
#' \item{HeckmantS(selectEq, outcomeEq, data = data, df, start)}{Heckman-t model adjustment function. The Heckman-t model
#' maintains the original parametric structure of the Classic Heckman model,
#' but considers a bivariate Student's t distribution as the underlying joint
#' distribution of the selection and primary regression variable and estimates
#' the parameters by maximum likelihood.}
#' \item{HeckmanSK(selectEq, outcomeEq, data = data, lambda, start)}{Heckman-SK model adjustment function. The Heckman-sk
#' model maintains the original parametric structure of the Classic Heckman
#' model, but considers a bivariate Skew-Normal distribution as the underlying
#' joint distribution of the selection and primary regression variable and
#' estimates the parameters by maximum likelihood.}
#' \item{HeckmanBS(selectEq, outcomeBS, data = data, start)}{Heckman-BS model adjustment function. The Heckman-BS model
#' maintains the original parametric structure of the Classic Heckman model,
#' but considers a bivariate Birnbaum-Saunders distribution as the underlying
#' joint distribution of the selection and primary regression variable and
#' estimates the parameters by maximum likelihood.}
#' \item{HeckmanGe(selectEq, outcomeEq,outcomeS, outcomeC, data = data)}{Function for adjustment of Generalized Heckman model. The
#' Generalized Heckman Model generalize the Classic Heckman model by adding
#' covariables to the dispersion and correlation parameters, which allows to
#' identify the covariates responsible for the presence of selection bias and
#' the presence of heteroscedasticity.}
#' }
#'
#' @importFrom stats binomial coef dnorm dt glm lm model.matrix model.response na.pass optim pnorm printCoefmat pt qnorm terms
#' @importFrom utils tail
#' @importFrom Rdpack reprompt
#'
#' @param selection Selection equation.
#' @param outcome Primary Regression Equation.
#' @param outcomeS Matrix with Covariables for fit of the Dispersion Parameter.
#' @param outcomeC Matrix with Covariates for Adjusting the Correlation Parameter.
#' @param df Initial value to the degree of freedom of Heckman-t model.
#' @param lambda Initial value for asymmetry parameter.
#' @param start initial values.
#' @param data Database.
#'
#' @return Applying any package function returns a list of results
#' that include estimates of the fit model parameters, hessian matrix,
#' number of observations, and more. If the initial value is not included
#' in the function argument, an initial value is estimated from the
#' Heckman two-step method setting.
#'
#' @seealso \code{\link{HeckmanCL}}
#' @seealso \code{\link{HeckmantS}}
#' @seealso \code{\link{HeckmanSK}}
#' @seealso \code{\link{HeckmanBS}}
#' @seealso \code{\link{HeckmanGe}}
#'
#' @author Fernando de Souza Bastos, Wagner Barreto de Souza
#'
#' @keywords Heckman
#'
#' @docType package
#' @name ssmodels
NULL

#' Medical Expenditure Panel Survey
#'
#'The MEPS is a set of large-scale surveys of families, individuals
#'and their medical providers (doctors, hospitals, pharmacies, etc.)
#'in the United States. It has data on the health services Americans use,
#'how often they use them, the cost of these services and how they are paid,
#'as well as data on the cost and reach of health insurance available
#'to American workers. The sample is restricted to persons aged between 21
#'and 64 years and contains a variable response with 3328 observations of
#'outpatient costs, of which 526 (15.8\%) correspond to unobserved expenditure
#'values and identified as zero expenditure for adjustment of the models.
#'It also includes the following explanatory variables:
#' \itemize{
#'   \item{educ: education status}
#'   \item{age: Age}
#'   \item{income: income}
#'   \item{female: gender}
#'   \item{vgood: a numeric vector}
#'   \item{good: a numeric vector}
#'   \item{hospexp: a numeric vector}
#'   \item{totchr: number of chronic diseases}
#'   \item{ffs: a numeric vector}
#'   \item{dhospexp: a numeric vector}
#'   \item{age2: a numeric vector}
#'   \item{agefem: a numeric vector}
#'   \item{fairpoor: a numeric vector}
#'   \item{year01: a numeric vector}
#'   \item{instype: a numeric vector}
#'   \item{ambexp: a numeric vector}
#'   \item{lambexp: log ambulatory expenditures}
#'   \item{blhisp: ethnicity}
#'   \item{instype_s1: a numeric vector}
#'   \item{dambexp: dummy variable, ambulatory expenditures}
#'   \item{lnambx: a numeric vector}
#'   \item{ins: insurance status}
#' }
#'
#' @source 2001 Medical Expenditure Panel Survey by the Agency
#' for Healthcare Research and Quality.
#'
#' @references{
#'   \insertRef{colin2009microeconometrics}{ssmodels}
#'
#'   \insertRef{ssmrob}{ssmodels}
#'
#'   \insertRef{sampleSelection}{ssmodels}
#' }
#'
#' @examples
#' data(MEPS2001)
#' attach(MEPS2001)
#' hist(lnambx)
#' selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
#' outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
#' HeckmanCL(selectEq, outcomeEq, data = MEPS2001)
"MEPS2001"

#' U.S. Women's Labor Force Participation
#'
#' The Mroz87 data frame contains data about 753
#' married women. These data are collected within the
#' "Panel Study of Income Dynamics" (PSID). Of the 753
#' observations, the first 428 are for women with positive
#' hours worked in 1975, while the remaining 325
#' observations are for women who did not work for
#' pay in 1975. A more complete discussion of the data
#' is found in \insertCite{mroz1987;textual}{ssmodels}. It also
#' includes the following explanatory variables:
#' \itemize{
#'   \item{lfp: Dummy variable for labor-force participation.}
#'   \item{hours: Wife's hours of work in 1975. }
#'   \item{kids5: Number of children 5 years old or younger.}
#'   \item{kids618: Number of children 6 to 18 years old.}
#'   \item{Age: Wife's age.}
#'   \item{Educ: Wife's educational attainment, in years.}
#'   \item{wage: Wife's average hourly earnings, in 1975 dollars.}
#'   \item{repwage: Wife's wage reported at the time of the 1976 interview.}
#'   \item{hushrs: Husband's hours worked in 1975.}
#'   \item{husage: Husband's age.}
#'   \item{huseduc: Husband's educational attainment, in years.}
#'   \item{huswage: Husband's wage, in 1975 dollars.}
#'   \item{faminc: Family income, in 1975 dollars.}
#'   \item{mtr: Marginal tax rate facing the wife.}
#'   \item{motheduc: Wife's mother's educational attainment, in years.}
#'   \item{fatheduc: Wife's father's educational attainment, in years.}
#'   \item{unem: Unemployment rate in county of residence, in percentage points.}
#'   \item{city: Dummy variable = 1 if live in large city, else 0. }
#'   \item{exper: Actual years of wife's previous labor market experience.}
#'   \item{nwifeinc: Non-wife income.}
#'   \item{wifecoll: Dummy variable for wife's college attendance.}
#'   \item{huscoll: Dummy variable for husband's college attendance.}
#' }
#' @source PSID Staff, The Panel Study of Income Dynamics,
#' Institute for Social ResearchPanel Study of Income
#' Dynamics, University of Michigan, \url{https://psidonline.isr.umich.edu/}
#'
#' @references{
#'   \insertRef{mroz1987}{ssmodels}
#'
#'   \insertRef{ssmrob}{ssmodels}
#'
#'   \insertRef{sampleSelection}{ssmodels}
#'
#'   \insertRef{wooldridge2016}{ssmodels}
#' }
#'
#' @examples
#' # Wooldridge(2016): page 247
#' data(Mroz87)
#' attach(Mroz87)
#' selectEq  <- lfp ~ nwifeinc + educ + exper + I(exper^2) + age + kids5 + kids618
#' outcomeEq <- log(wage) ~ educ + exper + I(exper^2)
#' outcomeS  <- cbind(educ, exper)
#' outcomeC  <- cbind(educ, exper)
#' outcomeBS <- wage ~ educ + exper + I(exper^2)
#' outcomeBS <- wage ~ educ + exper + I(exper^2)
#' HeckmanCL(selectEq, outcomeEq, data = Mroz87)
#' HeckmanBS(selectEq, outcomeBS, data = Mroz87)
#' HeckmanSK(selectEq, outcomeEq, data = Mroz87, lambda = 1)
#' HeckmantS(selectEq, outcomeEq, data = Mroz87, df=5)
#' HeckmanGe(selectEq, outcomeEq, outcomeS, outcomeC, data = Mroz87)
#'
"Mroz87"

#' RAND Health Insurance Experiment
#'
#' 'The RAND Health Insurance Experiment (RAND HIE) was a comprehensive
#' study of health care cost, utilization and outcome in the United States.
#' It is the only randomized study of health insurance, and the only study
#' which can give definitive evidence as to the causal effects of different
#' health insurance plans. For more information about the database visit:
#' \url{https://en.wikipedia.org/w/index.php?title=RAND_Health_Insurance_Experiment&oldid=110166949}
#' accessed september 09, 2019). This data frame contains the following columns:
#' \itemize{
#' \item{plan: HIE plan number.}
#'  \item{site: Participant's  place of residence
#'  when the participant was initially enrolled.}
#'  \item{coins: Coinsurance rate.}
#'  \item{tookphys: Took baseline physical.}
#'  \item{year: Study year.}
#'  \item{zper: Person identifier.}
#'  \item{black: 1 if race of household head is black.}
#'  \item{income: Family income.}
#'  \item{xage: Age in years.}
#'  \item{female: 1 if person is female.}
#'  \item{educdec: Education of household head in years.}
#'  \item{time: Time eligible during the year.}
#'  \item{outpdol: Outpatient expenses:
#'  all covered outpatient medical services
#'  excluding dental care, outpatient psychotherapy,
#'  outpatient drugs or supplies.}
#'  \item{drugdol: Drug expenses:
#'  all covered outpatient and dental drugs.}
#'  \item{suppdol: Supply expenses:
#'  all covered outpatient supplies including dental.}
#'  \item{mentdol: Psychotherapy expenses:
#'  all covered outpatient psychotherapy services including injections
#'  excluding charges for visits in excess of 52 per year,
#'  prescription drugs, and inpatient care.}
#'  \item{inpdol: Inpatient expenses:
#'  all covered inpatient expenses in a hospital, mental hospital,
#'  or nursing home,
#'  excluding outpatient care and renal dialysis.}
#'  \item{meddol: Medical expenses:
#'  all covered inpatient and outpatient services,
#'  including drugs, supplies, and inpatient costs of newborns
#'  excluding dental care and outpatient psychotherapy.}
#'  \item{totadm: Hospital admissions:
#'  annual number of covered hospitalizations.}
#'  \item{inpmis: Incomplete Hospital Records:
#'  missing inpatient records.}
#'  \item{mentvis: Psychotherapy visits:
#'  indicates the annual number of outpatient visits for psychotherapy.
#'  It includes billed visits only.
#'  The limit was 52 covered visits per person per year.
#'  The count includes an initial visit to a psychiatrist or psychologist.}
#'  \item{mdvis: Face-to-Face visits to physicians:
#'  annual covered outpatient visits with physician providers
#'  (excludes dental, psychotherapy, and
#'  radiology/anesthesiology/pathology-only visits).}
#'  \item{notmdvis: Face-to-Face visits to nonphysicians:
#'  annual covered outpatient visits with nonphysician providers
#'  such as speech and physical therapists, chiropractors,
#'  podiatrists, acupuncturists, Christian Science etc.
#'  (excludes dental, healers, psychotherapy,
#'  and radiology/anesthesiology/pathology-only visits).}
#'  \item{num: Family size.}
#'  \item{mhi: Mental health index.}
#'  \item{disea: Number of chronic diseases.}
#'  \item{physlm: Physical limitations.}
#'  \item{ghindx: General health index.}
#'  \item{mdeoff: Maximum expenditure offer.}
#'  \item{pioff: Participation incentive payment.}
#'  \item{child: 1 if age is less than 18 years.}
#'  \item{fchild: \code{female * child}.}
#'  \item{lfam: log of \code{num} (family size).}
#'  \item{lpi: log of \code{pioff} (participation incentive payment).}
#'  \item{idp: 1 if individual deductible plan.}
#'  \item{logc: \code{log(coins+1)}.}
#'  \item{fmde: 0 if \code{idp=1},
#'  \code{ln(max(1,mdeoff/(0.01*coins)))} otherwise.}
#'  \item{hlthg: 1 if self-rated health is good
#'  -- baseline is excellent self-rated health.}
#'  \item{hlthf: 1 if self-rated health is fair
#'  -- baseline is excellent self-rated health.}
#'  \item{hlthp: 1 if self-rated health is poor
#'  -- baseline is excellent self-rated health.}
#'  \item{xghindx: \code{ghindx} (general healt index)
#'  with imputations of missing values.}
#'  \item{linc: log of \code{income} (family income).}
#'  \item{lnum: log of \code{num} (family size).}
#'  \item{lnmeddol: log of \code{meddol} (medical expenses).}
#'  \item{binexp: 1 if \code{meddol} > 0.}
#'  }
#'
#' @source \url{http://cameron.econ.ucdavis.edu/mmabook/mmadata.html}
#'
#' @references{
#'   \insertRef{cameron2005}{ssmodels}
#'
#'   \insertRef{ssmrob}{ssmodels}
#'
#'   \insertRef{sampleSelection}{ssmodels}
#'
#'   \insertRef{wikiRand}{ssmodels}
#' }
#'
#' @keywords RandHIE
#'
#' @examples
#' ##Cameron and Trivedi (2005): Section 16.6
#' data(RandHIE)
#' subsample <- RandHIE$year == 2 & !is.na( RandHIE$educdec )
#' selectEq <- binexp ~ logc + idp + lpi + fmde + physlm + disea +
#'   hlthg + hlthf + hlthp + linc + lfam + educdec + xage + female +
#'   child + fchild + black
#'   outcomeEq <- lnmeddol ~ logc + idp + lpi + fmde + physlm + disea +
#'   hlthg + hlthf + hlthp + linc + lfam + educdec + xage + female +
#'   child + fchild + black
#'   cameron <- HeckmanCL(selectEq, outcomeEq, data = RandHIE[subsample, ])
#'   summary(cameron)
#'
"RandHIE"

#' Panel Study of Income Dynamics
#'
#' The data come from the Panel Study of Income
#' Dynamics, years 1981 to 1992 (also contains earnings
#' data from 1980). The sample consists of 579 white
#' females, who were followed over the considered period.
#' In total, there are 6,948 observations over the 12-year
#' period (1981-1992). This data frame contains the following
#' columns:
#'\itemize{
#'   \item{id: Individual identifier}
#'   \item{year: Survey year}
#'   \item{age: Calculated age in years (based on year and month of birth)}
#'   \item{educ: Years of schooling}
#'   \item{children: Total number of children in family unit, ages 0-17}
#'   \item{s: Participation dummy, =1 if worked (hours>0)}
#'   \item{lnw: Log of real average hourly earnings}
#'   \item{lnw80: Log earnings in 1980}
#'   \item{agesq: Age squared}
#'   \item{children_lag1: Number of children in t-1}
#'   \item{children_lag2: Number of children in t-2}
#'   \item{lnw2: Log of real average hourly earnings}
#'   \item{Lnw: Log of real average hourly earnings}
#' }
#' @source \url{http://simba.isr.umich.edu/}
#'
#' @references{
#'   \insertRef{semykina2013estimation}{ssmodels}
#'
#'   \insertRef{ssmrob}{ssmodels}
#'
#'   \insertRef{sampleSelection}{ssmodels}
#'}
#'
#' @examples
#' data(PSID2)
#' attach(PSID2)
#' hist(Lnw)
#' selectEq <- s ~ educ+ age+ children+ year
#' outcomeEq <- Lnw ~ educ+ age+ children
#' HCinitial(selectEq,outcomeEq, data = PSID2)
#' #Note that the estimated value of rho by the two-step
#' #method is greater than 1
#' summary(HeckmanGe(selectEq,outcomeEq, 1, 1, data = PSID2))
#'
"PSID2"

#' US National Health and Nutrition Examination Study
#'
#' The US National Health and Nutrition Examination Study (NHANES) is a
#' survey data collected by the US National Center for Health Statistics.
#' The survey data dates back to 1999, where individuals of all ages are
#' interviewed in their home annually and complete the health examination
#' component of the survey. The study variables include demographic
#' variables (e.g. age and annual household income), physical measurements
#' (e.g. BMI – body mass index), health variables (e.g. diabetes status),
#' and lifestyle variables (e.g. smoking status). This data frame contains
#' the following columns:
#'\itemize{
#'   \item{id: Individual identifier}
#'   \item{age: Age}
#'   \item{gender: Sex 1=male, 0=female}
#'   \item{educ: Education is dichotomized into high school
#'   and above versus less than high school}
#'   \item{race: categorical variable with five levels}
#'   \item{income: Household income ($1000 per year) was reported
#'   as a range of values in dollar (e.g. 0–4999, 5000–9999, etc.)
#'   and had 10 interval categories. }
#'   \item{Income: Household income ($1000 per year) was reported
#'   as a range of values in dollar (e.g. 0–4999, 5000–9999, etc.)
#'   and had 10 interval categories. }
#'   \item{bmi: body mass index}
#'   \item{sbp: systolic blood pressure}
#' }
#'
#' @source \url{https://wwwn.cdc.gov/nchs/nhanes/ContinuousNhanes/Default.aspx?BeginYear=2003}
#'
#' @references{
#'   \insertRef{ogundimu2019robust}{ssmodels}
#'
#'   \insertRef{little2011subsample}{ssmodels}
#'
#'   \insertRef{ssmrob}{ssmodels}
#'
#'   \insertRef{sampleSelection}{ssmodels}
#'}
#'
#' @examples
#' data("nhanes")
#' attach(nhanes)
#' hist(Income, prob= TRUE, breaks = seq(1, 99, 0.5), xlim = c(1,10),
#' ylim = c(0,0.35), main = "Histogram of Income", xlab = "Category")
#' data2 <- subset(nhanes, !is.na(sbp))
#' data3 <- subset(data2, !is.na(bmi))
#' attach(data3)
#' data <- data3
#' data$YS <- ifelse(is.na(data$Income),0,1)
#' data$educ <- ifelse(data$educ<=2,0,1)
#' attach(data)
#' selectionEq <- YS~age+gender+educ+race
#' outcomeEq   <- sbp~age+gender+educ+bmi
#'
"nhanes"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("If you have questions, suggestions,
  or comments regarding the 'ssmodels' package, please contact: fernando.bastos@ufv.br")
}

.onLoad <- function(lib, pkg){
   Rdpack::Rdpack_bibstyles(package = pkg, authors = "LongNames")
   invisible(NULL)
}
