# Plot phylogenetic tree with saasi results

Plot phylogenetic tree with saasi results

## Usage

``` r
plot_saasi(
  tree,
  saasi_result,
  colors = NULL,
  tip_cex = 0.5,
  node_cex = 0.2,
  save_file = NULL,
  width = 3000,
  height = 3000,
  res = 300
)
```

## Arguments

- tree:

  A phylo object containing tip.state

- saasi_result:

  Output from saasi()

- colors:

  Character vector of colors. If NULL, colors are automatically
  generated based on the number of states (default: NULL)

- tip_cex:

  Size of tip circles (default: 0.5)

- node_cex:

  Size of node pie charts (default: 0.2)

- save_file:

  Character. File path to save the plot (e.g. "tree.png"). If NULL, plot
  is drawn to current device (default: NULL)

- width:

  Plot width in pixels when saving (default: 3000)

- height:

  Plot height in pixels when saving (default: 3000)

- res:

  Plot resolution in dpi when saving (default: 300)

## Value

Invisibly returns the color vector used

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic plot
plot_saasi(tree, saasi_result)

# Custom colors and save to file
plot_saasi(tree, saasi_result,
           colors = c("red", "blue", "green"),
           node_cex = 0.3,
           save_file = "tree.png")
} # }
```
