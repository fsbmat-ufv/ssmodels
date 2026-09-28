# US National Health and Nutrition Examination Study

The US National Health and Nutrition Examination Study (NHANES) is a
survey data collected by the US National Center for Health Statistics.
The survey data dates back to 1999, where individuals of all ages are
interviewed in their home annually and complete the health examination
component of the survey. The study variables include demographic
variables (e.g. age and annual household income), physical measurements
(e.g. BMI – body mass index), health variables (e.g. diabetes status),
and lifestyle variables (e.g. smoking status). This data frame contains
the following columns:

- id: Individual identifier

- age: Age

- gender: Sex 1=male, 0=female

- educ: Education is dichotomized into high school and above versus less
  than high school

- race: categorical variable with five levels

- income: Household income (\$1000 per year) was reported as a range of
  values in dollar (e.g. 0–4999, 5000–9999, etc.) and had 10 interval
  categories.

- Income: Household income (\$1000 per year) was reported as a range of
  values in dollar (e.g. 0–4999, 5000–9999, etc.) and had 10 interval
  categories.

- bmi: body mass index

- sbp: systolic blood pressure

## Usage

``` r
nhanes
```

## Format

An object of class `data.frame` with 9643 rows and 9 columns.

## Source

<https://wwwn.cdc.gov/nchs/nhanes/ContinuousNhanes/Default.aspx?BeginYear=2003>

## References

Emmanuel O Ogundimu, Gary S Collins (2019). “A robust imputation method
for missing responses and covariates in sample selection models.”
*Statistical methods in medical research*, **28**(1), 102–116. Roderick
J Little, Nanhua Zhang (2011). “Subsample ignorable likelihood for
regression analysis with missing data.” *Journal of the Royal
Statistical Society: Series C (Applied Statistics)*, **60**(4), 591–605.
Mikhail Zhelonkin, Marc G. Genton, Elvezio Ronchetti (2019). *ssmrob:
Robust Estimation and Inference in Sample Selection Models*. R package
version 0.7, <https://CRAN.R-project.org/package=ssmrob>. Ott Toomet,
Arne Henningsen (2008). “Sample Selection Models in R: Package
sampleSelection.” *Journal of Statistical Software*, **27**(7).
<https://www.jstatsoft.org/article/view/v027i07>.

## Examples

``` r
data("nhanes")
attach(nhanes)
#> The following objects are masked from PSID2:
#> 
#>     age, educ, id
#> The following objects are masked from Mroz87 (pos = 4):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 5):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 6):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 7):
#> 
#>     age, educ, income
#> The following objects are masked from Mroz87 (pos = 8):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 9):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 10):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 11):
#> 
#>     age, educ, income
hist(Income, prob= TRUE, breaks = seq(1, 99, 0.5), xlim = c(1,10),
ylim = c(0,0.35), main = "Histogram of Income", xlab = "Category")

data2 <- subset(nhanes, !is.na(sbp))
data3 <- subset(data2, !is.na(bmi))
attach(data3)
#> The following objects are masked from nhanes:
#> 
#>     Income, age, bmi, educ, gender, id, income, race, sbp
#> The following objects are masked from PSID2:
#> 
#>     age, educ, id
#> The following objects are masked from Mroz87 (pos = 5):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 6):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 7):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 8):
#> 
#>     age, educ, income
#> The following objects are masked from Mroz87 (pos = 9):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 10):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 11):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 12):
#> 
#>     age, educ, income
data <- data3
data$YS <- ifelse(is.na(data$Income),0,1)
data$educ <- ifelse(data$educ<=2,0,1)
attach(data)
#> The following objects are masked from data3:
#> 
#>     Income, age, bmi, educ, gender, id, income, race, sbp
#> The following objects are masked from nhanes:
#> 
#>     Income, age, bmi, educ, gender, id, income, race, sbp
#> The following objects are masked from PSID2:
#> 
#>     age, educ, id
#> The following objects are masked from Mroz87 (pos = 6):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 7):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 8):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 9):
#> 
#>     age, educ, income
#> The following objects are masked from Mroz87 (pos = 10):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 11):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 12):
#> 
#>     age, educ, income
#> The following objects are masked from MEPS2001 (pos = 13):
#> 
#>     age, educ, income
selectionEq <- YS~age+gender+educ+race
outcomeEq   <- sbp~age+gender+educ+bmi
```
