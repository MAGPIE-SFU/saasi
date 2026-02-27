# Birth-death-sampling parameters for demonstration tree

A data frame containing the birth-death-sampling parameters used to
simulate the demonstration tree. These parameters control
diversification and sampling rates for each character state.

## Usage

``` r
demo_pars
```

## Format

A data frame with 3 rows (one per state) and 5 columns:

- state:

  Character state names (1, 2, or 3)

- rootprior:

  Prior probability of root being in each state (1/3 each)

- lambda:

  Speciation rate (1.5 for all states)

- mu:

  Extinction rate (0.3 for all states)

- psi:

  Sampling rate (0.1 for state 1, 0.5 for states 2 and 3)

## Examples

``` r
data(demo_pars)
demo_pars
#>   state rootprior lambda  mu psi
#> 1     1 0.3333333    1.5 0.3 0.1
#> 2     2 0.3333333    1.5 0.3 0.5
#> 3     3 0.3333333    1.5 0.3 0.5
# Compare sampling rates across states
demo_pars$psi
#> [1] 0.1 0.5 0.5
```
