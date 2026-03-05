# Simulate a birth/death/sampling tree with post-processing

This tree can be passed to
[saasi](https://magpie-sfu.github.io/saasi/reference/saasi.md), and will
include speciation, extinction, sampling, and mutation events. The tree
is post-processed to remove tips at the present and ensure a minimum
number of tips.

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

  Parameters that are used in the simulation. Should contain state
  names, speciation, extinction, and sampling rates,

- x0:

  Natural number used as root state in the returned tree. Must be a
  state declared in `params_df`.

- max_taxa:

  Maximum number of nodes allowed in the initial simulated tree.

- max_t:

  Maximum depth allowed in the returned tree.

- include_extinct:

  Boolean declaring whether extinct taxa are included in the returned
  tree.

- min_tip:

  Minimum number of tips required in the post-processed tree. If the
  processed tree has fewer tips, the simulation is repeated until this
  condition is satisfied. Default is 1.

- q_martix:

  Transition rate matrix.

## Value

A `phylo` phylogenetic tree (`ape` format) with post-processing applied.
The tree includes `tip.state` attribute with character state labels.
