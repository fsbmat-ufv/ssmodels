
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ssmodels <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

<!-- badges: end -->

<p style="text-align: justify;">

Com o objetivo de facilitar o ajuste dos modelos de seleção amostral
existentes na literatura, criamos o pacote ssmodels. Nosso pacote
permite o ajuste do modelo de Heckman clássico, e a estimação dos
parâmetros desse modelo via método de máxima verossimilhança e método
de dois passos, além do ajuste dos modelos Heckman-t, introduzido na
literatura em 2012 por @marchenko2012heckman e do modelo Heckman-Skew
introduzido na literatura por @ogundimu2016sample em 2016.
Implementamos, ainda, funções para o ajuste da versão generalizada do
modelo de Heckman que permite a inclusão de covariáveis aos parâmetros
de dispersão e correlação e uma função para o ajuste do modelo
Heckman-BS que utiliza a distribuição Birnbaum-Saunders como
distribuição conjunta das variáveis de seleção e regressão primária.

</p>

## Installation

You can install `ssmodels` from Github using the devtools package by
running the following code:

``` r
#install.package("devtools")
devtools::install_github("fsbmat-ufv/ssmodels")
```

If you want to build the package vignettes (recommended), you’ll want to
run the following code instead

``` r
#install.package("devtools")
devtools::install_github("fsbmat-ufv/ssmodels", build = TRUE, 
                         build_opts = c("--no-resave-data", "--no-manual"))
```
