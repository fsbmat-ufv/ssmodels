# Skew-Normal Sample Selection Model Fit Function

Fits a sample selection model based on the Skew-Normal distribution
using Maximum Likelihood Estimation (MLE). This model allows for
asymmetry in the distribution of the outcome variable's error term,
addressing potential skewness.

## Usage

``` r
HeckmanSK(
  selection,
  outcome,
  data = sys.frame(sys.parent()),
  lambda,
  start = NULL
)
```

## Arguments

- selection:

  A formula specifying the selection equation.

- outcome:

  A formula specifying the outcome equation.

- data:

  A data frame containing the variables.

- lambda:

  Initial start value for the skewness parameter (`lambda`).

- start:

  Optional numeric vector of initial parameter values.

## Value

A list containing:

- `coefficients`: Named vector of estimated model parameters.

- `value`: The (negative) log-likelihood at convergence.

- `loglik`: The maximum log-likelihood.

- `counts`: Number of gradient evaluations.

- `hessian`: Hessian matrix at the optimum.

- `fisher_infoSK`: Approximate Fisher information matrix.

- `prop_sigmaSK`: Standard errors for the estimates.

- `level`: Levels of the selection variable.

- `nObs`: Number of observations.

- `nParam`: Number of model parameters.

- `N0`: Number of censored (unobserved) observations.

- `N1`: Number of observed (uncensored) observations.

- `NXS`: Number of covariates in the selection equation.

- `NXO`: Number of covariates in the outcome equation.

- `df`: Degrees of freedom (observations minus parameters).

- `aic`: Akaike Information Criterion.

- `bic`: Bayesian Information Criterion.

- `initial.value`: Initial parameter values used.

## Details

The function implements MLE for a sample selection model where the
outcome equation's errors follow a Skew-Normal distribution, as proposed
in Ogundimu and Hutton (2016) . The optimization is performed via the
BFGS algorithm.

The results include estimates for:

- Selection equation coefficients.

- Outcome equation coefficients.

- Standard deviation of the error term (`sigma`).

- Correlation between the selection and outcome errors (`rho`).

- Skewness parameter (`lambda`).

- Robust standard errors from the Fisher information matrix.

## References

Emmanuel O Ogundimu, Jane L Hutton (2016). “A Sample Selection Model
with Skew-normal Distribution.” *Scandinavian Journal of Statistics*,
**43**(1), 172–190.

## Examples

``` r
data("Mroz87")
attach(Mroz87)
#> The following objects are masked from MEPS2001 (pos = 3):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 4):
#> 
#>     age, educ
#> The following objects are masked from MEPS2001 (pos = 5):
#> 
#>     age, educ
selectEq <- lfp ~ huswage + kids5 + mtr + fatheduc + educ + city
outcomeEq <- log(wage) ~ educ + city
HeckmanSK(selectEq, outcomeEq, data = Mroz87, lambda = -1.5)
#> Start not provided using default start values.
#> $coefficients
#> (Intercept)     huswage       kids5         mtr    fatheduc        educ 
#>  2.91662397 -0.09441350 -0.38395436 -5.13417100 -0.01560489  0.09497449 
#>        city (Intercept)        educ        city       sigma         rho 
#> -0.04295576  1.10262877  0.07462305  0.11951616  0.96402482 -0.79175493 
#>      lambda 
#> -1.58728652 
#> 
#> $value
#> [1] -874.6076
#> 
#> $loglik
#> [1] -874.6076
#> 
#> $counts
#> gradient 
#>       47 
#> 
#> $hessian
#>             [,1]        [,2]        [,3]        [,4]        [,5]         [,6]
#>  [1,]  -761.7616  -5697.7338  -174.90498  -518.25799  -6631.3241   -9260.0084
#>  [2,] -5697.7338 -54421.2864 -1333.95430 -3704.44565 -51487.4056  -71513.9636
#>  [3,]  -174.9050  -1333.9543  -233.90244  -120.46746  -1668.4387   -2284.7332
#>  [4,]  -518.2580  -3704.4457  -120.46746  -357.14648  -4464.8673   -6250.0363
#>  [5,] -6631.3241 -51487.4056 -1668.43866 -4464.86726 -67197.2071  -83188.1363
#>  [6,] -9260.0084 -71513.9636 -2284.73322 -6250.03628 -83188.1363 -116455.8353
#>  [7,]  -486.7713  -4141.7774  -107.87398  -323.36488  -4444.5593   -6065.8973
#>  [8,]  -382.2667  -2819.7009   -72.26997  -258.59364  -3357.1481   -4694.1550
#>  [9,] -4694.1532 -35735.3049  -969.70451 -3150.92536 -42431.0849  -59706.1584
#> [10,]  -242.6706  -2039.7765   -46.11459  -159.97497  -2255.0940   -3058.9824
#> [11,]   277.3505   1966.9866    58.74978   192.71995   2449.5816    3308.3822
#> [12,]   292.7946   1745.9665   -31.66055   196.14270   2832.1352    3794.3058
#> [13,]    51.7407    391.4114    10.60834    34.91793    447.0649     634.4167
#>              [,7]         [,8]         [,9]       [,10]       [,11]       [,12]
#>  [1,]  -486.77126   -382.26666   -4694.1532  -242.67057   277.35053   292.79462
#>  [2,] -4141.77744  -2819.70086  -35735.3049 -2039.77649  1966.98656  1745.96651
#>  [3,]  -107.87398    -72.26997    -969.7045   -46.11459    58.74978   -31.66055
#>  [4,]  -323.36488   -258.59364   -3150.9254  -159.97497   192.71995   196.14270
#>  [5,] -4444.55932  -3357.14810  -42431.0849 -2255.09397  2449.58163  2832.13518
#>  [6,] -6065.89729  -4694.15496  -59706.1584 -3058.98241  3308.38220  3794.30583
#>  [7,]  -486.77126   -242.67057   -3058.9814  -242.67057   174.47566   176.34119
#>  [8,]  -242.67057  -1088.85633  -13745.6741  -694.14240   831.19533  -380.90117
#>  [9,] -3058.98136 -13745.67406 -179513.5557 -8967.45103 10269.04808 -4458.84169
#> [10,]  -242.67057   -694.14240   -8967.4510  -694.14240   526.04740  -253.74430
#> [11,]   174.47566    831.19533   10269.0481   526.04740 -1554.48296  -136.09969
#> [12,]   176.34119   -380.90117   -4458.8417  -253.74430  -136.09969 -1085.77413
#> [13,]    33.11873    -14.90048    -235.9434   -15.25629   -66.88540   -86.95696
#>            [,13]
#>  [1,]   51.74070
#>  [2,]  391.41142
#>  [3,]   10.60834
#>  [4,]   34.91793
#>  [5,]  447.06493
#>  [6,]  634.41667
#>  [7,]   33.11873
#>  [8,]  -14.90048
#>  [9,] -235.94335
#> [10,]  -15.25629
#> [11,]  -66.88540
#> [12,]  -86.95696
#> [13,]  -25.84295
#> 
#> $fisher_infoSK
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,]  0.5061487522 -8.327244e-03 -4.013654e-03 -0.5371120660 -3.428459e-04
#>  [2,] -0.0083272436  2.709240e-04  2.744186e-04  0.0097748697 -1.270216e-06
#>  [3,] -0.0040136538  2.744186e-04  7.226074e-03  0.0060250463 -8.051846e-05
#>  [4,] -0.5371120660  9.774870e-03  6.025046e-03  0.6309429580  3.486240e-04
#>  [5,] -0.0003428459 -1.270216e-06 -8.051846e-05  0.0003486240  1.324573e-04
#>  [6,] -0.0046989830 -5.475479e-05 -4.440041e-04  0.0010794757 -7.721135e-05
#>  [7,] -0.0008806198 -2.612742e-04  2.981766e-04 -0.0004362881 -7.500774e-05
#>  [8,] -0.0590758978  1.185436e-03  5.864600e-03  0.0482390026 -1.611029e-04
#>  [9,]  0.0029630290 -4.715956e-05 -2.337884e-04 -0.0019006472  1.038548e-05
#> [10,] -0.0020417580  6.424387e-05  1.887737e-04  0.0030487025 -1.476144e-05
#> [11,] -0.0126689095  3.348703e-04  1.542570e-03  0.0147815255 -7.961139e-06
#> [12,]  0.0178885738 -5.470207e-04 -2.748362e-03 -0.0206870693  6.553555e-05
#> [13,]  0.0182748114 -2.666511e-04 -7.861374e-04 -0.0212752137 -1.609516e-04
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,] -4.698983e-03 -8.806198e-04 -0.0590758978  2.963029e-03 -2.041758e-03
#>  [2,] -5.475479e-05 -2.612742e-04  0.0011854355 -4.715956e-05  6.424387e-05
#>  [3,] -4.440041e-04  2.981766e-04  0.0058645997 -2.337884e-04  1.887737e-04
#>  [4,]  1.079476e-03 -4.362881e-04  0.0482390026 -1.900647e-03  3.048702e-03
#>  [5,] -7.721135e-05 -7.500774e-05 -0.0001611029  1.038548e-05 -1.476144e-05
#>  [6,]  4.486229e-04 -1.006737e-04  0.0007568160 -9.218447e-05  6.042436e-05
#>  [7,] -1.006737e-04  7.830085e-03  0.0006403754  7.180598e-05 -2.509722e-03
#>  [8,]  7.568160e-04  6.403754e-04  0.0542936791 -3.354060e-03 -2.454182e-04
#>  [9,] -9.218447e-05  7.180598e-05 -0.0033540597  2.448635e-04 -1.889947e-04
#> [10,]  6.042436e-05 -2.509722e-03 -0.0002454182 -1.889947e-04  5.034535e-03
#> [11,] -1.997339e-04 -4.472223e-05  0.0062881283 -2.045755e-04  3.554003e-04
#> [12,]  2.856504e-04 -5.578415e-07 -0.0089234359  3.722263e-04 -3.574583e-04
#> [13,]  5.129161e-04  6.525616e-04 -0.0029103955 -3.504477e-04 -1.728526e-03
#>               [,11]         [,12]         [,13]
#>  [1,] -1.266891e-02  1.788857e-02  0.0182748114
#>  [2,]  3.348703e-04 -5.470207e-04 -0.0002666511
#>  [3,]  1.542570e-03 -2.748362e-03 -0.0007861374
#>  [4,]  1.478153e-02 -2.068707e-02 -0.0212752137
#>  [5,] -7.961139e-06  6.553555e-05 -0.0001609516
#>  [6,] -1.997339e-04  2.856504e-04  0.0005129161
#>  [7,] -4.472223e-05 -5.578415e-07  0.0006525616
#>  [8,]  6.288128e-03 -8.923436e-03 -0.0029103955
#>  [9,] -2.045755e-04  3.722263e-04 -0.0003504477
#> [10,]  3.554003e-04 -3.574583e-04 -0.0017285264
#> [11,]  2.866391e-03 -2.261126e-03 -0.0065637206
#> [12,] -2.261126e-03  4.370226e-03 -0.0002994570
#> [13,] -6.563721e-03 -2.994570e-04  0.0767131540
#> 
#> $prop_sigmaSK
#>  [1] 0.71144132 0.01645977 0.08500632 0.79431918 0.01150901 0.02118072
#>  [7] 0.08848777 0.23301004 0.01564811 0.07095446 0.05353869 0.06610768
#> [13] 0.27697140
#> 
#> $level
#> [1] "0" "1"
#> 
#> $nObs
#> [1] 753
#> 
#> $nParam
#> [1] 13
#> 
#> $N0
#> [1] 325
#> 
#> $N1
#> [1] 428
#> 
#> $NXS
#> [1] 7
#> 
#> $NXO
#> [1] 3
#> 
#> $df
#> [1] 740
#> 
#> $aic
#> [1] 1775.215
#> 
#> $bic
#> [1] 1835.328
#> 
#> $initial.value
#>  [1]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0 -1.5
#> 
#> attr(,"class")
#> [1] "HeckmanSK" "list"     
```
