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
# Use the demo tree
demo_tree
#> 
#> Phylogenetic tree with 13 tips and 12 internal nodes.
#> 
#> Tip labels:
#>   sp2, sp4, sp8, sp9, sp10, sp26, ...
#> Node labels:
#>   nd1, nd2, nd6, nd22, nd27, nd74, ...
#> 
#> Rooted; includes branch length(s).

# Use the simmap function to estimate the rate transition matrix
# - Impose equal rates for the matrix structure
estimate_transition_rates(demo_tree_prepared, matrix_structure = "ER")
#>            1          2          3
#> 1 -1.7548671  0.8774336  0.8774336
#> 2  0.8774336 -1.7548671  0.8774336
#> 3  0.8774336  0.8774336 -1.7548671

if (FALSE) { # \dontrun{
# Use ace to estimate the rate transition matrix
# - Impose a structure such that transitions to or from the first state
#   happen at a different rate than second and third states
struct <- matrix(c(0, 1, 2,
                   1, 0, 2,
                   2, 2, 0), nrow=3, ncol=3)
estimate_transition_rates(demo_tree_prepared, matrix_structure = struct)
} # }
```
