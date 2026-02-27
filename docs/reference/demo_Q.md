# Transition rate matrix for demonstration tree

A 3x3 transition rate matrix (Q matrix) used to simulate the
demonstration tree. This matrix defines the instantaneous rates of
transition between character states.

## Usage

``` r
demo_Q
```

## Format

A 3x3 numeric matrix

## Examples

``` r
data(demo_Q)
demo_Q
#>      1    2    3
#> 1 -0.6  0.3  0.3
#> 2  0.3 -0.6  0.3
#> 3  0.3  0.3 -0.6
# Check that rows sum to 0
rowSums(demo_Q)
#> 1 2 3 
#> 0 0 0 
```
