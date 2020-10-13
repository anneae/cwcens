
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cwcens

<!-- badges: start -->

<!-- badges: end -->

This pacakage implements non-parametric estimation via a kernel
estimator for the probability in state and restricted mean time in state
in an illness-death model under component-wise censoring. Component-wise
censoring arises when illness can only be measured at a finite set of
times, while death is right censored and thus observed continuously up
to the right censoring time. Component-wise censored composite endpoints
arise often in biostatistical practice. For example, in many oncology
studies, progression-free survival is component-wise censored.

## Installation

You can install the development version of cwcens from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("anneae/cwcens")
```

## Example

First, we will simulate data from an illness-death model. Our dataset
will contain 150 patients. Since `scale21` is greater than zero, the
true model is reversible, meaning individuals can transition back to the
illness-free state (state 1) from the alive with illness state (state
2).

``` r
library(cwcens)
mydat <- simdat(150, scale12=1/.0008, scale13=1/.0002, scale23=1/.0016,
       scale21=1/.0008,
       vital.lfu=c(30.4*36, 30.4*48),
       visit.schedule = 30.4*c(6, 12, 18, 24, 30, 36, 42, 48),
       scatter.sd=10, seed = 678)

head(mydat)
#>       dtime dstatus state2obs laststate1       t1       t2       t3       t4
#> 1 1220.9207       0       Inf  1103.4099 175.8895 356.7024 547.5894 735.8449
#> 2  518.1648       1       Inf   366.3973 187.1082 366.3973       NA       NA
#> 3 1211.5844       1  359.2808   186.9515 186.9515 359.2808 527.3185 728.9607
#> 4 1428.0358       0  912.7325   715.8489 176.6697 361.7226 567.1672 715.8489
#> 5  661.0129       1  538.7025   373.8682 182.7077 373.8682 538.7025       NA
#> 6  778.3871       1  365.3571   191.6876 191.6876 365.3571 542.7675 743.0496
#>         t5       t6      t7 x1 x2 x3 x4 x5 x6 x7 nvisits
#> 1 916.4126 1103.410      NA  1  1  1  1  1  1 NA       7
#> 2       NA       NA      NA  1  1 NA NA NA NA NA       7
#> 3 912.3205 1107.106      NA  1  2  1  1  1  1 NA       7
#> 4 912.7325 1097.262 1276.44  1  1  1  1  2  2  2       7
#> 5       NA       NA      NA  1  1  2 NA NA NA NA       7
#> 6       NA       NA      NA  1  2  2  2 NA NA NA       7
```

The standard approach would …

The kernel approach does not …

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
