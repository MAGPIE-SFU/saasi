# Estimate the transition rate matrix of the discrete trait process

`estimate_transition_rates()` computes a maximum likelihood estimate of
the transition rate matrix \\Q\\ governing the discrete trait process.
The discrete traits evolve according to a Markov process with transition
rate matrix \\Q\\ such that element \\q\_{ij}\\ is the instantaneous
rate of transitioning from state \\i\\ to state \\j\\.

## Usage

``` r
estimate_transition_rates(tree, matrix_structure = "ER", method = "ace")
```

## Arguments

- tree:

  A phylo object with a tip.state attribute assigning traits to all
  tips.

- matrix_structure:

  The form of the transition rate matrix. Can either be a numeric matrix
  or one of the following strings: `"ER"`, `"SYM"`, and `"ARD"` (see
  details). The default value is `"ER"`.

- method:

  The method used to perform the estimation. Possible values are `"ace"`
  or `"fitMk"`.

## Value

A transition rate matrix that has been fit to the observed phylogeny.
This matrix will be compatible with other `saasi` functions.

## Details

The `model` parameters controls the structure of \\Q\\. `"ER"` specifies
equal rates for all transitions, \\q\_{ij}=q\\ for all traits \\i,j\\.
`"SYM"` specifies a symmetric transition rate matrix,
\\q\_{ij}=q\_{ji}\\ for all traits \\i,j\\. `"ARD"` places no structural
constraints and allows all traits to be different, \\q\_{ij}\ne
q\_{kl}\\ for all traits \\i,j,k,l\\. If another grouping of transition
rates is desired this can be input as a numeric matrix such that each
entry has a matrix with a positive integer and all entries with the same
value will share a rate.

The `method` parameter allows the user to select a preferred method
between `ace` and `fitMk`, however if the preferred method fails the
other will be used instead. If `ace` is selected, the maximum likelihood
estimator is computed using the
[`ape::ace()`](https://rdrr.io/pkg/ape/man/ace.html) function. If
`fitMk` is selected, the
[`phytools::fitMk()`](https://rdrr.io/pkg/phytools/man/fitMk.html)
function is used to fit \\Q\\ to the provided phylogenetic tree.

## Examples

``` r
# Load a timed phylogenetic tree for Ebola
data(ebola_tree)

# Use the simmap function to estimate the rate transition matrix
# - Impose a symmetric form on the matrix
Q <- estimate_transition_rates(ebola_tree, matrix_structure = "SYM")

# Use ace to estimate the rate transition matrix
# - Impose a structure such that transitions to or from Guinea happen at a different rate than Liberia or Sierra Leone
struct <- matrix(c(0, 1, 1, 1, 0, 2, 1, 2, 0), nrow=3, ncol=3)
Q <- estimate_transition_rates(ebola_tree, matrix_structure = struct)
```
