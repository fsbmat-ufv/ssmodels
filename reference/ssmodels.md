# ssmodels: Sample Selection Models in R

The `ssmodels` package provides functions to fit data affected by sample
selection bias. It includes several extensions of the classical Heckman
selection model, allowing for different assumptions about the joint
distribution of the selection and outcome equations.

## Details

The following models are implemented:

- `HeckmanCL`:

  Classic Heckman model (Tobit-2).

- `HeckmantS`:

  Heckman model with Student's t-distribution.

- `HeckmanSK`:

  Heckman model with Skew-Normal distribution.

- `HeckmanBS`:

  Heckman model with Birnbaum-Saunders distribution.

- `HeckmanGe`:

  Generalized Heckman model with covariates in the dispersion and
  correlation structures.

The package also includes helper functions for computing Inverse Mills
Ratios (IMR), post-processing parameter vectors, and two-step initial
value estimation.

## References

James J Heckman (1976). “The common structure of statistical models of
truncation, sample selection and limited dependent variables and a
simple estimator for such models.” In *Annals of Economic and Social
Measurement, Volume 5, number 4*, 475–492. NBER.

James J Heckman (1979). “Sample selection bias as a specification
error.” *Econometrica: Journal of the econometric society*, 153–161.

Thomas A Mroz (1987). “The sensitivity of an empirical model of married
women's hours of work to economic and statistical assumptions.”
*Econometrica: Journal of the Econometric Society*, 765–799.

Ott Toomet, Arne Henningsen (2008). “Sample Selection Models in R:
Package sampleSelection.” *Journal of Statistical Software*, **27**(7).
<https://www.jstatsoft.org/article/view/v027i07>.

Yulia V Marchenko, Marc G Genton (2012). “A Heckman selection-t model.”
*Journal of the American Statistical Association*, **107**(497),
304–317.

Emmanuel O Ogundimu, Jane L Hutton (2016). “A Sample Selection Model
with Skew-normal Distribution.” *Scandinavian Journal of Statistics*,
**43**(1), 172–190.

Mikhail Zhelonkin, Marc G Genton, Elvezio Ronchetti (2016). “Robust
inference in sample selection models.” *Journal of the Royal Statistical
Society: Series B (Statistical Methodology)*, **78**(4), 805–827.

Mikhail Zhelonkin, Marc G. Genton, Elvezio Ronchetti (2019). *ssmrob:
Robust Estimation and Inference in Sample Selection Models*. R package
version 0.7, <https://CRAN.R-project.org/package=ssmrob>.

Emmanuel O Ogundimu, Gary S Collins (2019). “A robust imputation method
for missing responses and covariates in sample selection models.”
*Statistical methods in medical research*, **28**(1), 102–116.

Fernando de Souza Bastos, Wagner Barreto-Souza (2020).
“Birnbaum–Saunders sample selection model.” *Journal of Applied
Statistics*.

Fernando de Souza Bastos, Wagner Barreto-Souza, Marc G Genton (2022). “A
Generalized Heckman Model With Varying Sample Selection Bias and
Dispersion Parameters.” *Statistica Sinica*.

## See also

[`HeckmanCL`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmanCL.md),
[`HeckmantS`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmantS.md),
[`HeckmanSK`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmanSK.md),
[`HeckmanBS`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmanBS.md),
[`HeckmanGe`](https://fsbmat-ufv.github.io/ssmodels/reference/HeckmanGe.md)

## Author

Fernando de Souza Bastos, Wagner Barreto de Souza
