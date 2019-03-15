# msir

[![CRAN\_Status](http://www.r-pkg.org/badges/version/msir)](https://cran.r-project.org/package=msir)
[![CRAN\_MonthlyDownloads](http://cranlogs.r-pkg.org/badges/msir)](https://cran.r-project.org/package=msir)
[![CRAN\_dependencies](https://tinyverse.netlify.com/badge/msir)](https://cran.R-project.org/package=mclust)

An [R](https://www.r-project.org/) package implementing *model-based sliced inverse regression* as described in Scrucca (2011).

Model-based Sliced Inverse Regression (MSIR) is a dimension reduction method based on Gaussian finite mixture models which provides an extension to sliced inverse regression (SIR). 

The basis of the MSIR subspace is estimated by modeling the inverse distribution within slice using Gaussian finite mixtures with number of components and covariance matrix parameterization selected by BIC or defined by the user.

## Installation

You can install the released version of **msir** from CRAN using:

```{r}
install.packages("msir")
```

or the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("luca-scr/msir", build_vignettes = TRUE)
```

## Usage

Usage of the main function and some examples are included in vignette **A quick tour of msir**, which is available as

```{r}
vignette("msir")
```

The vignette is also available in the *Vignette* section on the navigation bar on top of the package's web page.

## References

Scrucca, L. (2011) Model-based SIR for dimension reduction. *Computational Statistics & Data Analysis*, 55(11), 3010-3026.
