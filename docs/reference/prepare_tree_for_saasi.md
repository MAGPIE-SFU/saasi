# Reformat a phylogenetic tree to be compatible with `saasi`

`prepare_tree_for_saasi()` will take in a `phylo` object and reformat it
so that the output is compabtible with the `saasi` function. This
reformatting includes:

- Adding tip annotations

- Resolving polytomies

- Imposing a minimum branch length

- Removing undesired states

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

  An object of class `phylo`.

- tip_data:

  Either a named vector or a `data.frame`. If a named vector, the
  elements of the vector must be the state of each tip and the names
  must be the labels of those tips. If a `data.frame`, there must be 2
  columns: the first with the tip labels and the second with the tip
  states. This can be omitted if `tree` already has a `tip.states`
  attribute.

- resolve_polytomies:

  Boolean. Indicates whether polytomies should be resolved. Default
  value is `TRUE`. Setting to `FALSE` will produce a tree that is
  incompatible with `saasi`.

- fix_branches:

  Boolean. Indicates whether zero and negative branch lengths should be
  fixed. Default value is `TRUE`. Setting to `FALSE` will produce a tree
  that is incompatible with `saasi`.

- min_branch_length:

  Numeric. The minimum branch length for the resultant tree. Any short
  branch lengths will be increased accordingly. Default value is `1e-5`.

- drop_states:

  An optional character vector of states to remove from the tree.

## Value

An object of class `phylo` that is compatible with the
[`saasi()`](https://magpie-sfu.github.io/saasi/reference/saasi.md)
function.

## Details

The reformatting is conducted as follows:

- Polytomies are resolved using
  [`ape::multi2di()`](https://rdrr.io/pkg/ape/man/multi2di.html).

- The length of any branches with zero or negative length is set to
  `min_branch_length`.

- All tips with an `NA` value for their state annotation are removed in
  addition to any states specified in `drop_states`.

## See also

[`add_tip_states()`](https://magpie-sfu.github.io/saasi/reference/add_tip_states.md)
to manually add state annotations to tip states or
[`drop_tips_by_state()`](https://magpie-sfu.github.io/saasi/reference/drop_tips_by_state.md)
to manually remove tips from the tree based on their state.

## Examples

``` r
# Check if the demo tree is compatible
check_tree_compatibility(demo_tree)
#> Tree is not compatible with SAASI
#> NA tip states present. Remove with drop_tips_by_state(tree, NA)
#> No tip states. Attach with attach_tip_states() or prepare_tree_for_saasi()
#> [1] FALSE

# Since tip annotations are missing, we need state data for the tips
head(demo_metadata)
#>   node state
#> 1  sp2     1
#> 2  sp4     2
#> 3  sp8     2
#> 4  sp9     1
#> 5 sp10     3
#> 6 sp26     2

# Reformat the demo tree with the tip states
tree_prepared <- prepare_tree_for_saasi(demo_tree, demo_metadata)

# Check the compability of the reformatted tree
check_tree_compatibility(tree_prepared)
#> Tree is compatible with SAASI
#> [1] TRUE
```
