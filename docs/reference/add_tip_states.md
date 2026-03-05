# Add tip states to a phylogenetic tree

`add_tip_states()` adds tip states to a `phylo` object. This function is
used within the
[`prepare_tree_for_saasi()`](https://magpie-sfu.github.io/saasi/reference/prepare_tree_for_saasi.md)
function.

## Usage

``` r
add_tip_states(tree, tip_data)
```

## Arguments

- tree:

  An object of class `phylo`.

- tip_data:

  Either a named vector or a `data.frame`. If a named vector, the
  elements of the vector must be the state of each tip and the names
  must be the labels of those tips. If a `data.frame`, there must be 2
  columns: the first with the tip labels and the second with the tip
  states.

## Value

An object of class `phylo` that matches the `tree` input with an
additional attribute `tip.state` that contains the tip states.

## See also

[`drop_tips_by_state()`](https://magpie-sfu.github.io/saasi/reference/drop_tips_by_state.md)
to remove specific traits or
[`prepare_tree_for_saasi()`](https://magpie-sfu.github.io/saasi/reference/prepare_tree_for_saasi.md)
for reformatting `phylo` objects to be `saasi`-compatible.
