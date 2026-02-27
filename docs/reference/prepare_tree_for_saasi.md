# Prepare phylogenetic tree for saasi

Prepare phylogenetic tree for saasi

## Usage

``` r
prepare_tree_for_saasi(
  tree,
  tip_data = NULL,
  resolve_polytomies = TRUE,
  fix_branches = TRUE,
  min_branch_length = 1e-05,
  drop_states = character(0)
)
```

## Arguments

- tree:

  A phylo object

- tip_data:

  A named vector where names = tip labels, values = states, or a df with
  2 columns: column 1 = tip labels, column 2 = states. Optional if tree
  already has tip states

- resolve_polytomies:

  Resolve polytomies with multi2di

- fix_branches:

  Fix zero/negative branch lengths

- min_branch_length:

  Replace value for zero or negative branch lengths (default: 1e-5)

- drop_states:

  States to be removed from the tree. NA is always included. Add
  additional values to drop states. Example:c("Not Collected").

## Value

A prepared phylo object compatible to SAASI.
