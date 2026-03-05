# Simulate a phylogenetic tree from a birth-death-sampling process

`sim_bds_tree()` simulates a phylogenetic tree with state annotations
for all internal nodes and tips. A birth-death-sampling process is used
to generate the tree topology and branch lengths; this will include
speciation, sampling, and mutation events. A continuous-time Markov
process is used to assign states to all nodes. The resulting tree will
be compatible with
[`saasi()`](https://magpie-sfu.github.io/saasi/reference/saasi.md). The
tree is post-processed to remove tips at the present and ensure a
minimum number of tips.

## Usage

``` r
sim_bds_tree(
  params_df,
  q_matrix,
  x0,
  max_taxa = 100,
  max_t = 100,
  include_extinct = FALSE,
  min_tip = 1
)
```

## Arguments

- params_df:

  A `data.frame` specifying the parameters from which the tree will be
  simulated. Must contain

  - `lambda`: vector of speciation rates

  - `mu`: vector of extinction rates

  - `psi`: vector of sampling rates

  The number of rows should be equal to the number of states.

- q_matrix:

  State transition rate matrix.

- x0:

  The root state.

- max_taxa:

  Maximum number of nodes allowed in the initial simulated tree. Default
  value is 100.

- max_t:

  Maximum depth allowed in the returned tree. Default value is 100.

- include_extinct:

  Boolean. If `TRUE`, extinct taxa are included in the returned tree.
  Default value is `FALSE`.

- min_tip:

  Minimum number of tips required in the post-processed tree. If the
  processed tree has fewer tips, the simulation is repeated until this
  condition is satisfied. Default value is 1.

## Value

An object of class `phylo` that is compatible with
[`saasi()`](https://magpie-sfu.github.io/saasi/reference/saasi.md).

## Examples

``` r
# Simulate a tree from the demo model with root state set to 1
demo_pars
#>   state rootprior lambda  mu psi
#> 1     1 0.3333333    1.5 0.3 0.1
#> 2     2 0.3333333    1.5 0.3 0.5
#> 3     3 0.3333333    1.5 0.3 0.5
demo_Q
#>      1    2    3
#> 1 -0.6  0.3  0.3
#> 2  0.3 -0.6  0.3
#> 3  0.3  0.3 -0.6
sim_bds_tree(demo_pars, demo_Q, 1)
#> 
#> Phylogenetic tree with 37 tips and 36 internal nodes.
#> 
#> Tip labels:
#>   sp1, sp2, sp3, sp4, sp5, sp6, ...
#> Node labels:
#>   nd1, nd2, nd9, nd14, nd19, nd20, ...
#> 
#> Rooted; includes branch length(s).
```
