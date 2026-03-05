# SAASI demonstration: pre-prepared tree

A simulated birth-death-sampling tree with 13 tips and 3 states. The
tree has pre-prepared to be compatible with `saasi`.

## Usage

``` r
demo_tree_prepared
```

## Format

An object of class `phylo` with `tip.state` attribute

## Examples

``` r
data(demo_tree_prepared)
demo_tree_prepared$tip.state
#>  [1] "1" "2" "2" "1" "3" "2" "3" "2" "3" "2" "3" "2" "3"
```
