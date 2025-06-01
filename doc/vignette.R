## ----setup, include=FALSE-----------------------------------------------------
library(ssmodels)
library(ggplot2)
require(knitr)
require(kfigr)
library(kableExtra)
options(knitr.table.format = "latex")
knitr::opts_chunk$set(echo = TRUE, fig.align = "center", message=FALSE, 
                      warning=FALSE, fig.height=5, fig.width=7.2)

## ----fig4, fig.align= "center", fig.width = 7.2, fig.height= 5, anchor="Figura", fig.alt= "Log of Expenditures Medical"----
library(ssmodels)
#Leitura do dados MEPS2001
data(MEPS2001)
#tornando visiveis as colunas do data-frame
attach(MEPS2001)
barfill <- "grey"
barlines <- "black"
p1 <- ggplot(MEPS2001,aes(ambexp))+geom_histogram(colour = barlines, fill = barfill)+
  scale_x_continuous(name = "(a) Expenditures Medical",
                     breaks = seq(0, 15000, 2500),
                     limits=c(0, 15000))+
  scale_y_continuous(name = "Count",
                     breaks = seq(0, 800, 100),
                     limits=c(0, 800))
p2 <- ggplot(MEPS2001,aes(lambexp))+geom_histogram(colour = barlines, fill = barfill)+
  scale_x_continuous(name = "(b) Log of Expenditures Medical",
                     breaks = seq(0, 11, 1),
                     limits=c(0, 11))+
  scale_y_continuous(name = "Count",
                     breaks = seq(0, 300, 100),
                     limits=c(0, 300))
gridExtra::grid.arrange(p1, p2, ncol=2)

## ----warning=FALSE------------------------------------------------------------
selectEq <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins 
outcomeS <- cbind(age,female,totchr,ins)
outcomeC <- 1
outcomeBS <- ambexp ~ age + female + educ + blhisp + totchr + ins 
mCL <- HeckmanCL(selectEq, outcomeEq, data = MEPS2001)
mBS <- HeckmanBS(selectEq, outcomeBS, data = MEPS2001)
mSK <- HeckmanSK(selectEq, outcomeEq, data = MEPS2001,lambda = 1)
mtS <- HeckmantS(selectEq, outcomeEq, data = MEPS2001,df=12)
mGe <- HeckmanGe(selectEq, outcomeEq,outcomeS, outcomeC, data = MEPS2001)
Parameters <- c("Intercept", "age", "female", "educ", "blhisp", "totchr", "ins", "income",
   "Intercept", "age", "female", "educ", "blhisp", "totchr", "ins", "sigma", "age", "female", "totchr", "ins", "rho", "nu", "lambda")
HBS <- round(mBS$coefficients, digits = 3)
HCL <- round(mCL$coefficients, digits = 3)
HSK <- round(mSK$coefficients, digits = 3)
HtS <- round(mtS$coefficients, digits = 3)
HGe <- round(mGe$coefficients, digits = 3)
Results <- data.frame("Parameters"= Parameters,
                      "HeckmanGe" = c(HGe[1:21], "NA", "NA"), 
                      "HeckmanCL" = c(HCL[1:16], "NA", "NA", "NA", "NA", HCL[17], "NA", "NA"),
                      "HeckmanBS" = c(HBS[1:16], "NA", "NA", "NA", "NA", HBS[17], "NA", "NA"), 
                      "HeckmantS" = c(HtS[1:16], "NA", "NA", "NA", "NA", HtS[17:18], "NA" ),
                      "HeckmanSK" = c(HSK[1:16], "NA", "NA", "NA", "NA", HSK[17], "NA", HSK[18]))
kable(Results, format = "html", align = c("c", "c", "c", "c", "c"))

## ----warning=FALSE------------------------------------------------------------
summary(mCL)

## ----warning=FALSE------------------------------------------------------------
summary(mGe)

## ----warning=FALSE------------------------------------------------------------
summary(mBS)

## ----warning=FALSE------------------------------------------------------------
summary(mtS)

## ----warning=FALSE------------------------------------------------------------
summary(mSK)

## ----warning=FALSE------------------------------------------------------------
library(ssmodels)
data(nhanes)
attach(nhanes)
perc <- function(x,data){
    nna <- ifelse(sum(is.na(x))!=0,summary(x)[[7]],0)
    perc <- ifelse(sum(is.na(x))!=0,(nna/length(data$id))*100,0)
   return(perc)
}
Variables <- c("SBP (mm Hg)", "Age (year)", "Gender", "BMI (Kg/$m^{2}$)", "Education (years)", "Race", "Income ($\\$1000$ per year)", "Numbers Obs.")
perc1 <- round(perc(sbp,nhanes), digits = 2)
perc2 <- round(perc(age,nhanes), digits = 2)
perc3 <- round(perc(gender,nhanes), digits = 2)
perc4 <- round(perc(bmi,nhanes), digits = 2)
perc5 <- round(perc(educ,nhanes), digits = 2)
perc6 <- round(perc(race,nhanes), digits = 2)
perc7 <- round(perc(Income,nhanes), digits = 2)
nObs <- length(Income)
Percentage <- c(perc1, perc2, perc3, perc4, perc5, perc6, perc7, nObs)
df <- subset(nhanes, !is.na(sbp))
df <- subset(df, !is.na(bmi))
attach(df)
perc11 <- round(perc(sbp,df), digits = 2)
perc12 <- round(perc(age,df), digits = 2)
perc13 <- round(perc(gender,df), digits = 2)
perc14 <- round(perc(bmi,df), digits = 2)
perc15 <- round(perc(educ,df), digits = 2)
perc16 <- round(perc(race,df), digits = 2)
perc17 <- round(perc(Income,df), digits = 2)
nObs1 <- length(Income)
Percentage1 <- c(perc11, perc12, perc13, perc14, perc15, perc16, perc17, nObs1)
table <- data.frame("Variables" = Variables, "Percentage of Missing"= Percentage, "Without Missing"= Percentage1)
kable(table, format = "html", align = c("c", "c", "c"))

## ----warning=FALSE, fig.alt="Log Systolic blood pressure"---------------------
library("gridExtra")
barfill <- "grey"
barlines <- "black"
p1 <- ggplot(df, aes(Income)) + geom_histogram( breaks = seq(0, 10, 0.5), aes(y = ..density..), colour = barlines, fill = barfill)+
  scale_x_continuous(name = "Income",
                   breaks = seq(0, 10, 2),
                   limits=c(0, 10))

p2 <- ggplot(df, aes(log(sbp))) + geom_histogram( colour = barlines, fill = barfill) +
  scale_x_continuous(name = "Log Systolic blood pressure")
grid.arrange(p1, p2, ncol=2)

## ----warning = FALSE----------------------------------------------------------
df$YS <- ifelse(is.na(df$Income),0,1)
df$educ <- ifelse(df$educ<=2,0,1)
df$Income <- ifelse(is.na(df$Income),0,df$Income)
attach(df)
selectionEq <- YS~age+gender+educ+race
outcomeEq   <- log(sbp)~age+gender+educ+bmi+Income
outcomeBS   <- sbp~age+gender+educ+bmi+Income
mCL <- HeckmanCL(selectionEq, outcomeEq, data = df)
mBS <- HeckmanBS(selectionEq, outcomeBS, data = df)
mSK <- HeckmanSK(selectionEq, outcomeEq, data = df, lambda = 0)
mtS <- HeckmantS(selectionEq, outcomeEq, data = df, df = 15)

    Parameters <- c("Intercept", "age", "gender", "educ", "race", "Intercept", "age", "gender", "educ", "bmi", "income", "sigma", "rho", "nu", "lambda")
HBS <- round(mBS$coefficients, digits = 5)
HCL <- round(mCL$coefficients, digits = 5)
HSK <- round(mSK$coefficients, digits = 5)
HtS <- round(mtS$coefficients, digits = 5)
Results <- data.frame("Parameters"= Parameters, 
                      "HeckmanCL" = c(HCL[1:13], "NA", "NA"),
                      "HeckmanBS" = c(HBS[1:13], "NA", "NA"), 
                      "HeckmantS" = c(HtS[1:13], HtS[14], "NA"),
                      "HeckmanSK" = c(HSK[1:13], "NA", HSK[14]))
kable(Results, format = "html", align = c("c", "c", "c", "c", "c"))

## ----warning = FALSE----------------------------------------------------------
summary(mCL)

## ----warning = FALSE----------------------------------------------------------
summary(mtS)

## ----warning = FALSE----------------------------------------------------------
summary(mBS)

## ----warning = FALSE----------------------------------------------------------
summary(mSK)

