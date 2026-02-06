
<!-- README.md is generated from README.Rmd. Please edit that file -->

``` r
library(saasi)
library(diversitree)
#> Loading required package: ape
library(tidytree)
#> If you use the ggtree package suite in published research, please cite the appropriate paper(s):
#> 
#> LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR Jones, T Bradley, H Zhu, Y
#> Guan, Y Jiang, G Yu. treeio: an R package for phylogenetic tree input and output with richly
#> annotated and associated data. Molecular Biology and Evolution. 2020, 37(2):599-603. doi:
#> 10.1093/molbev/msz240
#> 
#> Shuangbin Xu, Lin Li, Xiao Luo, Meijun Chen, Wenli Tang, Li Zhan, Zehan Dai, Tommy T. Lam, Yi
#> Guan, Guangchuang Yu. Ggtree: A serialized data object for visualization of a phylogenetic tree
#> and annotation data. iMeta 2022, 1(4):e56. doi:10.1002/imt2.56
#> 
#> Attaching package: 'tidytree'
#> The following objects are masked from 'package:ape':
#> 
#>     drop.tip, keep.tip
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(ape)
library(phytools)
#> Loading required package: maps
library(readr)
```

# Introduction

<!-- badges: start -->
<!-- badges: end -->

This vignette demonstrates how to use the `saasi` package for ancestral
state reconstruction.

Saasi (Sampling-Aware Ancestral State Inference) is an ancestral state
reconstruction method that accounts for variation in sampling rates
among locations or traits. Unlike traditional methods that assume
uniform sampling, saasi explicitly models heterogeneous sampling rates,
leading to more accurate ancestral state estimates in unevenly sampled
phylogenies.

NOTE: THIS DOCUMENT IS SUBJECT TO CHANGE (FUNCTION NAMES, ARGUMENTS
ETC).

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

2.  **Birth-death-sampling parameters** (class `data.frame`):

    - Speciation rate
    - Extinction/removal rate
    - Sampling rate

    Use `estimate_bds_parameters()` to estimate speciation and sampling
    rates. Users should specify the extinction rate based on domain
    knowledge.

3.  **Transition rate matrix Q** (class `matrix`):

    - Rates of transition between states

    Use `estimate_transition_rates()` to estimate transition rates.

The output is a data frame containing the probability of each state for
each internal node of the phylogenetic tree.

## Example: Ebola 2013-2016 West African Ebola Epidemic

### Read tree and load metadata

For this example, we use data from Nextstrain:
<https://nextstrain.org/ebola/ebov-2013?c=country>

``` r
tree <- ape::read.tree(system.file("extdata", "nextstrain_ebola_ebov-2013_timetree.nwk", package = "saasi"))
metadata <- readr::read_tsv(system.file("extdata", "nextstrain_ebola_ebov-2013_metadata.tsv", package = "saasi"))
#> Rows: 1493 Columns: 8
#> ── Column specification ─────────────────────────────────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr  (7): strain, country, division, author, author__url, accession, accession__url
#> date (1): date
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

First, create a data frame containing strains and states. The first
column of the data frame should match `tree$tip.label`, and the second
column should contain the state of interest. Users can skip this step if
`tree$tip.state` already exists.

``` r
tip_data <- data.frame(
  tip_label = metadata$strain,
  state = metadata$country
)
```

### Modify the tree to be compatible with saasi

Users can apply the function `prepare_tree_for_saasi()` to make the tree
compatible with saasi. To check if the tree is compatible, use the
function `check_tree_compatibility()`.

``` r
ebola_tree <- prepare_tree_for_saasi(tree, tip_data)
```

### Estimating transition rates

Users can estimate the transition rate matrix Q using the function
`estimate_transition_rates()`. Users can specify the model for the
transition rate matrix: 1. Equal rate `ER`, 2. Symmetric rate `SYM`, 3.
All rates different `ARD`, and Custom `custom_q`.

``` r
Q <- estimate_transition_rates(ebola_tree, method = 'simmap', matrix_structure = 'SYM')
#> make.simmap is sampling character histories conditioned on
#> the transition matrix
#> 
#> Q =
#>                  Guinea    Liberia Sierra Leone
#> Guinea       -0.4716929  0.2602826    0.2114103
#> Liberia       0.2602826 -0.3695978    0.1093152
#> Sierra Leone  0.2114103  0.1093152   -0.3207255
#> (estimated using likelihood);
#> and (mean) root node prior probabilities
#> pi =
#>       Guinea      Liberia Sierra Leone 
#>    0.3333333    0.3333333    0.3333333
#> Done.
```

### Estimating speciation and sampling rates

Users can estimate speciation and sampling rates using the function
`estimate_bds_parameters()`. Users should have some knowledge about the
expected removal rate (1/mu) for the disease of interest, as saasi
requires the extinction rate mu to be known. To obtain better estimates
for the overall speciation and sampling rates, users should have
knowledge about the maximum and minimum R0, as well as the upper and
lower bounds for the total removal rate (1/(mu+psi)).

For Ebola, we assume the total infectious period (1/(mu+psi)) is between
20 and 40 days. Converting to years, the total removal rate (mu + psi)
is between 9.125 (365/40) and 18.25 (365/20). If we assume mu = 5, this
gives us bounds for the overall sampling rate psi between 4.125 and
13.25. In this example, we assume R0 is between 1.5 and 3.

``` r
rates <- estimate_bds_parameters(
    ebola_tree,
    mu = 5, 
    r0_max = 3, 
    r0_min = 1.5,
    psi_max = 15,
    infectious_period_min = 20/365, # convert days to years
    infectious_period_max = 40/365, # convert days to years
    n_starts = 100)
```

### Run saasi analysis

Users can run the saasi analysis by specifying the following parameters:

- **phy** - A `phylo` object containing the phylogenetic tree
- **q_matrix** - An n × n transition rate matrix (class `matrix`), where
  n is the number of states
- **lambda** - Speciation rate(s) for each state. If
  `length(lambda) = 1`, all states are assumed to have the same
  speciation rate; otherwise, provide a vector of length n
- **mu** - Extinction rate(s) for each state. If `length(mu) = 1`, all
  states are assumed to have the same extinction rate; otherwise,
  provide a vector of length n
- **psi** - Sampling rate(s) for each state. If `length(psi) = 1`, all
  states are assumed to have the same sampling rate; otherwise, provide
  a vector of length n
- **prior** - Prior probabilities for the root node. If `NULL`
  (default), uses equal probabilities across all states

``` r
saasi_ebola <- saasi(phy = ebola_tree,            # phylogenetic tree
                     q_matrix = Q,                # transition rate matrix
                     lambda = rates$lambda,       # speciation rate
                     mu = rates$mu,               # extinction rate
                     psi = rates$psi,             # sampling rate
                     prior = NULL)                # root prior (equal probabilities)
#> Tree is compatible with SAASI
```

### Plot and save the result

Users can use the built-in function `plot_saasi()` to visualize results.
Set `save_file = NULL` to display without saving.

``` r
p1 <- plot_saasi(ebola_tree, saasi_ebola, save_file = NULL)
```

<img src="man/figures/README-Plot the result-1.png" width="100%" />

### Run saasi with a different set of parameters

Users can specify different parameters for different locations. For
example, if we assume that Liberia samples at half the rate of other
countries:

``` r
pars2 <- create_params_template(colnames(Q), lambda = rates$lambda, mu = rates$mu, psi = c(rates$psi, rates$psi/2, rates$psi))

saasi_ebola2 <- saasi(phy = ebola_tree,                           # phylogenetic tree
                     q_matrix = Q,                                # transition rate matrix
                     lambda = rates$lambda,                       # speciation rate
                     mu = rates$mu,                               # extinction rate
                     psi = c(rates$psi, rates$psi/2, rates$psi),  # sampling rate
                     prior = NULL)  
#> Tree is compatible with SAASI

p2 <- plot_saasi(ebola_tree, saasi_ebola2, save_file = NULL)
```

<img src="man/figures/README-Run saasi with a different set of parameters-1.png" width="100%" />

### Interpreting results

The `saasi()` function returns a data frame with probability
distributions over states for each internal node. Higher probabilities
indicate greater confidence in that ancestral state. When comparing
equal versus unequal sampling rates, you may notice differences in
ancestral state probabilities, particularly for nodes ancestral to
under-sampled locations like Liberia.
