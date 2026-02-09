
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saasi

<!-- badges: start -->
<!-- badges: end -->

This vignette demonstrates how to use the `saasi` package for ancestral
state reconstruction.

## Installation

You can install the development version of saasi from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("MAGPIE-SFU/saasi")
```

# Overview

The `saasi` function requires three main inputs:

1.  **A phylogenetic tree** (class `phylo`) that is:

    - Rooted and binary
    - Has branch lengths in units of time (all positive)
    - Contains tip states (`tree$tip.state`) with no missing values

    Use `check_tree_compatibility()` to verify compatibility and
    `prepare_tree_for_saasi()` to prepare your tree.

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
