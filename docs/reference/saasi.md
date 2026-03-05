# Sampling Aware Ancestral State Inference

Reconstructs ancestral states for internal nodes of a phylogenetic tree
while accounting for heterogeneous sampling rates across states or
locations.

## Usage

``` r
saasi(phy, Q, pars, sensitivity_test = FALSE)
```

## Arguments

- phy:

  A `phylo` object containing the phylogenetic tree. The tree must be
  rooted, binary, and contain `tip.state`.

- Q:

  A numeric transition rate matrix \\\exp{(n \times n)}\\ where n is the
  number of states. Row and column names represent states. Off-diagonal
  elements are transition rates; diagonal elements are set such that
  rows sum to zero.

- pars:

  A data frame with the other parameters used in the ancestral state
  reconstruction algorithm. Must have the following column names: state,
  prior, lambda, mu, and psi. The prior values refer to the baseline
  probabilities of the states (used at the root of the tree).

- sensitivity_test:

  A boolean indicating if a sensitivity test should be conducted. The
  default value is `FALSE`.

## Value

A data frame with state probabilities for each internal node in `phy`.

## Examples

``` r
head(demo_pars) #contains state, lambda, mu, psi for each state
#>   state rootprior lambda  mu psi
#> 1     1 0.3333333    1.5 0.3 0.1
#> 2     2 0.3333333    1.5 0.3 0.5
#> 3     3 0.3333333    1.5 0.3 0.5
result <- saasi(phy = demo_tree_prepared, 
                Q = demo_Q,
                pars = demo_pars) 
#> Tree is compatible with SAASI
head(result) 
#>              1         2            3
#> 14 0.693604070 0.2318290 0.0745668922
#> 15 0.263875644 0.6429590 0.0931653929
#> 16 0.001350377 0.9972886 0.0013609914
#> 17 0.000821105 0.9982972 0.0008816567
#> 18 0.004837719 0.9862182 0.0089440393
#> 19 0.004747732 0.9433805 0.0518717196
```
