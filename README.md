# ssmodels

You can install `ssmodels` from Github using the devtools package by running the following code:

```{r}
#install.package("devtools")
devtools::install_github("fsbmat/ssmodels")
```

If you want to build the package vignettes (recommended), you'll want to run the following code instead

```{r}
#install.package("devtools")
devtools::install_github("fsbmat/ssmodels", build = TRUE, 
                         build_opts = c("--no-resave-data", "--no-manual"))
```
