# Sampling-Aware Ancestral State Inference

Reconstructs ancestral states for internal nodes of a phylogenetic tree
while accounting for heterogeneous sampling rates across states or
locations.

## Usage

``` r
saasi(phy, Q, pars, sensitivity_test = FALSE)
```

## Arguments

- phy:

  An object of class `phylo` specifying the phylogenetic tree. The tree
  must be rooted, binary, and contain the `tip.state` attribute. Use the
  [`check_tree_compatibility()`](https://magpie-sfu.github.io/saasi/reference/check_tree_compatibility.md)
  function to verify that `phy` is satisfies these conditions.

- Q:

  An \\n\times n\\ transition rate matrix, where \\n\\ is the number of
  states. Row and column names must correspond to the states.
  Off-diagonal elements are the transition rates between states;
  diagonal elements are set such that rows sum to zero.

- pars:

  A `data.frame` with the other parameters used in the ancestral state
  reconstruction algorithm. Must have the following columns:

  - `state`: a vector of state names. These should match the column and
    row names of `Q`.

  - `lambda`: a vector of birth rates for the birth-death-sampling
    process.

  - `mu`: a vector of death rates for the birth-death-sampling process.

  - `psi`: a vector of sampling rates for the birth-death-sampling
    process.

  It can also have `root_prior` that specifies the *a priori*
  distribution of the state of the root; if one is not provided a
  uniform distribution over the root state will be used.

- sensitivity_test:

  A boolean indicating if a sensitivity test should be conducted (see
  Details). The default value is `FALSE`.

## Value

A `data.frame` with a state probabilities for each internal node in
`phy`.

## Details

If `sensitivity_test=TRUE`, then a test of the sensitivity of `saasi` to
the birth-death-sampling rate parameters is conducted. The values of
\\\lambda\\, \\\mu\\, and \\\psi\\ are perturbed by sampling values at
random from +/- 10% of the input value. The ancestral state
reconstruction is then run on the perturbed parameters; this is repeated
10 times. To summarise sensitivity, the modal state for each internal
node is first computed for each perturbed reconstruction. Then the total
number of differences in modal states are computed between the perturbed
reconstructions and the original reconstruction. The average number of
these total differences across all 10 perturbed reconstructions is
reported.

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
