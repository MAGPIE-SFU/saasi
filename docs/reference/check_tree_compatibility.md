# Check if phylogenetic tree is compatible with `saasi`

Verifies the compatibility of a `phylo` object with the
[`saasi()`](https://magpie-sfu.github.io/saasi/reference/saasi.md)
function without modifying the object. For `tree` to be compatible with
[`saasi()`](https://magpie-sfu.github.io/saasi/reference/saasi.md), it
must:

- be rooted,

- have the correct number of internal nodes (the number of tips - 1),

- have no unary nodes or polytomies,

- have positive branch lengths, and

- have valid tip state annotations for all tips.

If the tree is not compatible with
[`saasi()`](https://magpie-sfu.github.io/saasi/reference/saasi.md), the
user can run
[`prepare_tree_for_saasi()`](https://magpie-sfu.github.io/saasi/reference/prepare_tree_for_saasi.md)
to resolve the issues.

## Usage

``` r
check_tree_compatibility(tree)
```

## Arguments

- tree:

  An object of class `phylo`.

## Value

Logical. Returns `TRUE` if the tree is compatible with `saasi` and
`FALSE` otherwise.

## Examples

``` r
# Check if this demo tree is compatible with saasi:
check_tree_compatibility(demo_tree)
#> Tree is not compatible with SAASI
#> NA tip states present. Remove with drop_tips_by_state(tree, NA)
#> No tip states. Attach with attach_tip_states() or prepare_tree_for_saasi()
#> [1] FALSE

# Now check another tree which has been properly prepared:
check_tree_compatibility(demo_tree_prepared)
#> Tree is compatible with SAASI
#> [1] TRUE
```
