# Drop tips by state

Removes tips that matches any of the specified values.

## Usage

``` r
drop_tips_by_state(tree, drop_values)
```

## Arguments

- tree:

  A phylo object with tip.state

- drop_values:

  Character vector of states to remove. Include NA to drop tips with NA
  states. Example: c("Not Collected", NA)

## Value

A phylo object with the specified tips removed and tip.state correctly
updated.
