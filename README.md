
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ssmodels <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

<!-- badges: end -->

The goal of ssmodels is to …

## Installation

You can install `ssmodels` from Github using the devtools package by
running the following code:

``` r
#install.package("devtools")
devtools::install_github("fsbmat-ufv/ssmodels")
#> Downloading GitHub repo fsbmat-ufv/ssmodels@HEAD
#> quantreg (5.67  -> 5.75) [CRAN]
#> Rdpack   (1.0.0 -> 2.1 ) [CRAN]
#> Installing 2 packages: quantreg, Rdpack
#> Installing packages into 'C:/Users/Dell/AppData/Local/Temp/Rtmp0OZ47w/temp_libpath34e023ef316f'
#> (as 'lib' is unspecified)
#> package 'quantreg' successfully unpacked and MD5 sums checked
#> package 'Rdpack' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\Dell\AppData\Local\Temp\RtmpQZCo2q\downloaded_packages
#>          checking for file 'C:\Users\Dell\AppData\Local\Temp\RtmpQZCo2q\remotes2ca493a71e0\fsbmat-ufv-ssmodels-f63608a/DESCRIPTION' ...  v  checking for file 'C:\Users\Dell\AppData\Local\Temp\RtmpQZCo2q\remotes2ca493a71e0\fsbmat-ufv-ssmodels-f63608a/DESCRIPTION' (395ms)
#>       -  preparing 'ssmodels':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   v  checking DESCRIPTION meta-information
#>       -  installing the package to process help pages
#>      Loading required namespace: ssmodels
#>       -  saving partial Rd database
#>       -  checking for LF line-endings in source and make files and shell scripts
#>       -  checking for empty or unneeded directories
#>       -  building 'ssmodels_0.0.0.9000.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/Dell/AppData/Local/Temp/Rtmp0OZ47w/temp_libpath34e023ef316f'
#> (as 'lib' is unspecified)
```

If you want to build the package vignettes (recommended), you’ll want to
run the following code instead

``` r
#install.package("devtools")
devtools::install_github("fsbmat-ufv/ssmodels", build = TRUE, 
                         build_opts = c("--no-resave-data", "--no-manual"))
#> Skipping install of 'ssmodels' from a github remote, the SHA1 (f63608a3) has not changed since last install.
#>   Use `force = TRUE` to force installation
```
