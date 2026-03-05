# Drop tips from a phylogenetic tree by their state values

`drop_tips_by_state()` will remove from a `phylo` object all tips in a
specified state(s).

## Usage

``` r
drop_tips_by_state(tree, drop_values)
```

## Arguments

- tree:

  An object of class `phylo` with a `tip.state` attribute. See
  [`add_tip_states()`](https://magpie-sfu.github.io/saasi/reference/add_tip_states.md)
  for adding tip annotations to a `phylo` object.

- drop_values:

  Character vector of states to remove.

## Value

An object of class `phylo` with the specified states removed.

## See also

[`add_tip_states()`](https://magpie-sfu.github.io/saasi/reference/add_tip_states.md)
to add tip annotations

## Examples

``` r
# This demo tree has 13 tips in 3 possible states
demo_tree_prepared
#> 
#> Phylogenetic tree with 13 tips and 12 internal nodes.
#> 
#> Tip labels:
#>   sp2, sp4, sp8, sp9, sp10, sp26, ...
#> Node labels:
#>   nd1, nd2, nd6, nd22, nd27, nd74, ...
#> 
#> Rooted; includes branch length(s).
table(demo_tree_prepared$tip.state)
#> 
#> 1 2 3 
#> 2 6 5 

# Remove all tips in state 3
new_tree <- drop_tips_by_state(demo_tree_prepared, "3")

# The resulting tree has 8 tips, having removed all tips in state 3
new_tree
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   sp2, sp4, sp8, sp9, sp26, sp30, ...
#> Node labels:
#>   nd1, nd2, nd6, nd22, nd27, nd32, ...
#> 
#> Rooted; includes branch length(s).
table(new_tree$tip.state)
#> 
#> 1 2 
#> 2 6 
```
