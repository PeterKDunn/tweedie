
# tweedie <a href="https://CRAN.R-project.org/package=tweedie/"><img src="man/figures/tweedie_hex.png" align="right" height="138" /></a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- ![R-CMD-check](https://github.com/rstudio/rmarkdown/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rstudio/rmarkdown/actions/workflows/R-CMD-check.yaml)
[![CRAN release](https://www.r-pkg.org/badges/version/rmarkdown)](https://cran.r-project.org/package=rmarkdown)
[![Codecov test coverage](https://codecov.io/gh/rstudio/rmarkdown/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rstudio/rmarkdown?branch=main) -->

<!-- badges: end -->

The **tweedie** package allows likelihood computations for Tweedie
distributions.

Apart from special cases (the normal, Poisson, gamma, inverse Gaussian
distributions), Tweedie distributions do not have closed-form density
functions or distribution functions. This package uses fast numerical
algorithms (infinite oscillating integrals; infinite series) to evaluate
the Tweedie density functions and distribution functions.

## Installation

You can install the development version of **tweedie** from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("PeterKDunn/tweedie")
#> ℹ Loading metadata database
#> ✔ Loading metadata database ... done
#> 
#> 
#> → Package library at '/tmp/RtmpXD4av3/temp_libpath17e217b4a6586'.
#> → Will install 4 packages.
#> → Will update 1 package.
#> → Will download 4 CRAN packages (1.63 MB).
#> → Will download 1 package with unknown size.
#> + cli               3.6.6 [bld][cmp][dl] (644.13 kB)
#> + lifecycle         1.0.5 [bld][dl] (107.14 kB)
#> + rlang             1.3.0 [bld][cmp][dl] (784.08 kB)
#> + statmod           1.5.2 [bld][cmp][dl] (93.23 kB)
#> + tweedie   3.2.0 → 3.2.0 [bld][cmp][dl] (GitHub: 944a2b1)
#> ℹ Getting 4 pkgs (1.63 MB) and 1 pkg with unknown size
#> ✔ Got lifecycle 1.0.5 (source) (107.14 kB)
#> ✔ Got rlang 1.3.0 (source) (784.08 kB)
#> ✔ Got statmod 1.5.2 (source) (93.23 kB)
#> ✔ Got tweedie 3.2.0 (source) (358.74 kB)
#> ✔ Got cli 3.6.6 (source) (644.13 kB)
#> ℹ Building cli 3.6.6
#> ℹ Building rlang 1.3.0
#> ℹ Building statmod 1.5.2
#> ✔ Built statmod 1.5.2 (1.9s)
#> ✔ Installed statmod 1.5.2  (20ms)
#> ✔ Built cli 3.6.6 (6.7s)
#> ✔ Installed cli 3.6.6  (1s)
#> ✔ Built rlang 1.3.0 (7.9s)
#> ✔ Installed rlang 1.3.0  (32ms)
#> ℹ Building lifecycle 1.0.5
#> ✔ Built lifecycle 1.0.5 (904ms)
#> ✔ Installed lifecycle 1.0.5  (1s)
#> ℹ Packaging tweedie 3.2.0
#> ✔ Packaged tweedie 3.2.0 (493ms)
#> ℹ Building tweedie 3.2.0
#> ✔ Built tweedie 3.2.0 (3.6s)
#> ✔ Installed tweedie 3.2.0 (github::PeterKDunn/tweedie@944a2b1) (1s)
#> ✔ 1 pkg + 4 deps: upd 1, added 4, dld 5 (NA B) [19.6s]
```

To also build and install the vignette (this takes a little longer, and
installs **knitr** and **rmarkdown** if needed):

``` r
options(pkg.build_vignettes = TRUE)
pak::pak("PeterKDunn/tweedie", dependencies = TRUE)
#> 
#> → Package library at '/tmp/RtmpXD4av3/temp_libpath17e217b4a6586'.
#> ℹ No downloads are needed
#> ✔ 1 pkg + 4 deps: kept 5 [972ms]
vignette("tweedie", package = "tweedie")
#> Warning: vignette 'tweedie' not found
```

Alternatively, using **remotes**:

``` r
remotes::install_github("PeterKDunn/tweedie", build_vignettes = TRUE, dependencies = TRUE)
#> Using GitHub PAT from the git credential store.
#> Skipping install of 'tweedie' from a github remote, the SHA1 (944a2b16) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

`install.packages("tweedie")` includes the vignette automatically.

## Tweedie distributions

Tweedie distributions are exponential dispersion models, with a mean
$\mu$ and a variance $\phi \mu^\xi$, for some dispersion parameter
$\phi > 0$ and a power index $\xi$ (sometimes called $p$) that uniquely
defines the distribution within the Tweedie family (for all real values
of $\xi$ **not** between 0 and 1).

Special cases of the Tweedie distributions are:

- the *normal* distribution, with $\xi = 0$ (i.e., the variance is
  $\phi$ and not related to the mean);
- the *Poisson* distribution, with $\xi = 1$ and $\phi = 1$ (i.e., the
  variance is the same as the mean);
- the *gamma* distribution, with $\xi = 2$; and
- the *inverse Gaussian* distribution, with $\xi = 3$.

For all other values of $\xi$, the probability functions and
distribution functions have no closed forms.

For $\xi < 1$, applications are limited (non-existent so far?), but have
support on the entire real line and $\mu > 0$.

For $1 < \xi < 2$, Tweedie distributions can be represented as a Poisson
sum of gamma distributions. These distributions are continuous for
$Y > 0$ but have a discrete mass at $Y = 0$.

For $\xi \ge 2$, the distributions have support on the positive reals.

The vignette contains examples.

<!-- REMEMBER: commit and push the resulting figure files, so they display on GitHub and CRAN. -->
