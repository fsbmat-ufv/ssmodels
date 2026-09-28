# RAND Health Insurance Experiment

'The RAND Health Insurance Experiment (RAND HIE) was a comprehensive
study of health care cost, utilization and outcome in the United States.
It is the only randomized study of health insurance, and the only study
which can give definitive evidence as to the causal effects of different
health insurance plans. For more information about the database visit:
<https://en.wikipedia.org/w/index.php?title=RAND_Health_Insurance_Experiment&oldid=110166949>
accessed september 09, 2019). This data frame contains the following
columns:

- plan: HIE plan number.

- site: Participant's place of residence when the participant was
  initially enrolled.

- coins: Coinsurance rate.

- tookphys: Took baseline physical.

- year: Study year.

- zper: Person identifier.

- black: 1 if race of household head is black.

- income: Family income.

- xage: Age in years.

- female: 1 if person is female.

- educdec: Education of household head in years.

- time: Time eligible during the year.

- outpdol: Outpatient expenses: all covered outpatient medical services
  excluding dental care, outpatient psychotherapy, outpatient drugs or
  supplies.

- drugdol: Drug expenses: all covered outpatient and dental drugs.

- suppdol: Supply expenses: all covered outpatient supplies including
  dental.

- mentdol: Psychotherapy expenses: all covered outpatient psychotherapy
  services including injections excluding charges for visits in excess
  of 52 per year, prescription drugs, and inpatient care.

- inpdol: Inpatient expenses: all covered inpatient expenses in a
  hospital, mental hospital, or nursing home, excluding outpatient care
  and renal dialysis.

- meddol: Medical expenses: all covered inpatient and outpatient
  services, including drugs, supplies, and inpatient costs of newborns
  excluding dental care and outpatient psychotherapy.

- totadm: Hospital admissions: annual number of covered
  hospitalizations.

- inpmis: Incomplete Hospital Records: missing inpatient records.

- mentvis: Psychotherapy visits: indicates the annual number of
  outpatient visits for psychotherapy. It includes billed visits only.
  The limit was 52 covered visits per person per year. The count
  includes an initial visit to a psychiatrist or psychologist.

- mdvis: Face-to-Face visits to physicians: annual covered outpatient
  visits with physician providers (excludes dental, psychotherapy, and
  radiology/anesthesiology/pathology-only visits).

- notmdvis: Face-to-Face visits to nonphysicians: annual covered
  outpatient visits with nonphysician providers such as speech and
  physical therapists, chiropractors, podiatrists, acupuncturists,
  Christian Science etc. (excludes dental, healers, psychotherapy, and
  radiology/anesthesiology/pathology-only visits).

- num: Family size.

- mhi: Mental health index.

- disea: Number of chronic diseases.

- physlm: Physical limitations.

- ghindx: General health index.

- mdeoff: Maximum expenditure offer.

- pioff: Participation incentive payment.

- child: 1 if age is less than 18 years.

- fchild: `female * child`.

- lfam: log of `num` (family size).

- lpi: log of `pioff` (participation incentive payment).

- idp: 1 if individual deductible plan.

- logc: `log(coins+1)`.

- fmde: 0 if `idp=1`, `ln(max(1,mdeoff/(0.01*coins)))` otherwise.

- hlthg: 1 if self-rated health is good – baseline is excellent
  self-rated health.

- hlthf: 1 if self-rated health is fair – baseline is excellent
  self-rated health.

- hlthp: 1 if self-rated health is poor – baseline is excellent
  self-rated health.

- xghindx: `ghindx` (general healt index) with imputations of missing
  values.

- linc: log of `income` (family income).

- lnum: log of `num` (family size).

- lnmeddol: log of `meddol` (medical expenses).

- binexp: 1 if `meddol` \> 0.

## Usage

``` r
RandHIE
```

## Format

An object of class `data.frame` with 20190 rows and 45 columns.

## Source

<https://cameron.econ.ucdavis.edu/mmabook/mmadata.html>

## References

A Colin Cameron, Pravin K Trivedi (2005). *Microeconometrics: methods
and applications*. Cambridge university press. Mikhail Zhelonkin, Marc
G. Genton, Elvezio Ronchetti (2019). *ssmrob: Robust Estimation and
Inference in Sample Selection Models*. R package version 0.7,
<https://CRAN.R-project.org/package=ssmrob>. Ott Toomet, Arne Henningsen
(2008). “Sample Selection Models in R: Package sampleSelection.”
*Journal of Statistical Software*, **27**(7).
<https://www.jstatsoft.org/article/view/v027i07>. Wikipedia contributors
(2019). “RAND Health Insurance Experiment — Wikipedia, The Free
Encyclopedia.”
<https://en.wikipedia.org/w/index.php?title=RAND_Health_Insurance_Experiment&oldid=909771077>.
\[Online; accessed 9-September-2019\].

## Examples

``` r
##Cameron and Trivedi (2005): Section 16.6
data(RandHIE)
subsample <- RandHIE$year == 2 & !is.na( RandHIE$educdec )
selectEq <- binexp ~ logc + idp + lpi + fmde + physlm + disea +
  hlthg + hlthf + hlthp + linc + lfam + educdec + xage + female +
  child + fchild + black
  outcomeEq <- lnmeddol ~ logc + idp + lpi + fmde + physlm + disea +
  hlthg + hlthf + hlthp + linc + lfam + educdec + xage + female +
  child + fchild + black
  cameron <- HeckmanCL(selectEq, outcomeEq, data = RandHIE[subsample, ])
#> Start not provided using default start values.
  summary(cameron)
#> 
#> --------------------------------------------------------------
#>           Classic Heckman Model (Package: ssmodels)           
#> --------------------------------------------------------------
#> --------------------------------------------------------------
#> Maximum Likelihood estimation 
#> optim function with method BFGS - iterations number: 100 
#> Log-Likelihood: -10177.96 
#> AIC: 20431.91 BIC: 20683.7 
#> Number of observations: ( 1293 censored and 4281 observed ) 
#> 38 free parameters ( df = 5536 ) 
#> --------------------------------------------------------------
#> Probit selection equation:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -0.358349   0.185832  -1.928 0.053863 .  
#> logc        -0.100715   0.026690  -3.773 0.000163 ***
#> idp         -0.117050   0.051454  -2.275 0.022954 *  
#> lpi          0.023570   0.008709   2.706 0.006823 ** 
#> fmde         0.001125   0.016017   0.070 0.943990    
#> physlm       0.286320   0.073040   3.920 8.96e-05 ***
#> disea        0.020799   0.003532   5.888 4.14e-09 ***
#> hlthg        0.051141   0.043179   1.184 0.236303    
#> hlthf        0.196376   0.081825   2.400 0.016430 *  
#> hlthp        0.803745   0.210368   3.821 0.000135 ***
#> linc         0.056038   0.016671   3.361 0.000781 ***
#> lfam        -0.039668   0.040779  -0.973 0.330721    
#> educdec      0.039659   0.007604   5.216 1.90e-07 ***
#> xage         0.001170   0.002130   0.549 0.582934    
#> female       0.417754   0.053738   7.774 9.01e-15 ***
#> child        0.109820   0.079482   1.382 0.167124    
#> fchild      -0.420083   0.079016  -5.316 1.10e-07 ***
#> black       -0.556754   0.052312 -10.643  < 2e-16 ***
#> --------------------------------------------------------------
#> Outcome equation:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  2.348848   0.244461   9.608  < 2e-16 ***
#> logc        -0.062249   0.032335  -1.925  0.05426 .  
#> idp         -0.122898   0.063282  -1.942  0.05218 .  
#> lpi          0.016640   0.010002   1.664  0.09624 .  
#> fmde        -0.028474   0.018522  -1.537  0.12427    
#> physlm       0.341393   0.071793   4.755 2.03e-06 ***
#> disea        0.026711   0.003647   7.324 2.76e-13 ***
#> hlthg        0.137045   0.049604   2.763  0.00575 ** 
#> hlthf        0.441995   0.090812   4.867 1.16e-06 ***
#> hlthp        0.935635   0.178087   5.254 1.55e-07 ***
#> linc         0.107146   0.022165   4.834 1.37e-06 ***
#> lfam        -0.155876   0.047356  -3.292  0.00100 ** 
#> educdec      0.014025   0.008712   1.610  0.10751    
#> xage         0.006427   0.002320   2.770  0.00562 ** 
#> female       0.511342   0.061874   8.264  < 2e-16 ***
#> child       -0.178092   0.093013  -1.915  0.05558 .  
#> fchild      -0.533367   0.094176  -5.663 1.56e-08 ***
#> black       -0.475875   0.074851  -6.358 2.21e-10 ***
#> --------------------------------------------------------------
#> Error terms:
#>       Estimate Std. Error t value Pr(>|t|)    
#> sigma  1.47497    0.02969   49.68   <2e-16 ***
#> rho    0.62935    0.05926   10.62   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> --------------------------------------------------------------
```
