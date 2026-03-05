
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saasi

<!-- badges: start --> <!-- badges: end -->

Saasi is an ancestral state reconstruction method that accounts for
variation in sampling rates among locations or traits.

## Installation

You can install the development version of saasi from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")

remotes::install_github("MAGPIE-SFU/saasi")
```

## Building documentation

This package uses roxygen2 and pkgdown.

To generate RD files:

``` r
devtools::document()
```

To build the pkgdown website:

``` r
pkgdown::build_site()
```

To rebuild the README.md:

``` r
devtools::build_readme()
```

Do **not edit** the `README.md` file, instead work on the `README.rmd`
and generate a new `README.md`.
