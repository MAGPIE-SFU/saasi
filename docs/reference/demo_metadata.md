# Metadata for demonstration tree

A data frame containing the node name and tip states

## Usage

``` r
demo_metadata
```

## Format

A data frame with 13 rows (one per tip) and 2 columns

- node:

  node name for each tip

- state:

  Character state for each tip (1, 2, or 3)

## Examples

``` r
data(demo_metadata)
demo_metadata
#>    node state
#> 1   sp2     1
#> 2   sp4     2
#> 3   sp8     2
#> 4   sp9     1
#> 5  sp10     3
#> 6  sp26     2
#> 7  sp27     3
#> 8  sp30     2
#> 9  sp31     3
#> 10 sp59     2
#> 11 sp64     3
#> 12 sp72     2
#> 13 sp75     3
demo_metadata$states
#> NULL
```
