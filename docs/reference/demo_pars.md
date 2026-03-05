# SAASI demonstration: birth-death-sampling process parameters

A `data.frame` specifying the parameters of the birth-death-sampling
process that is used to simulate the demonstration tree. These
parameters control diversification and sampling rates for each state.

## Usage

``` r
demo_pars
```

## Format

A data frame with 3 rows (one per state) and 5 columns:

- `state`:

  The character label for each state

- `rootprior`:

  Prior probability distribution over states for the root node

- `lambda`:

  Numeric speciation rate for each state

- `mu`:

  Numeric extinction rate for each state

- `psi`:

  Numeric Sampling rate for each state

## Examples

``` r
data(demo_pars)
demo_pars
#>   state rootprior lambda  mu psi
#> 1     1 0.3333333    1.5 0.3 0.1
#> 2     2 0.3333333    1.5 0.3 0.5
#> 3     3 0.3333333    1.5 0.3 0.5
```
