# Estimate rate parameters of a birth-death-sampling process

`estimate_bds_parameters()` estimates the speciation rate \\\lambda\\
and sampling rate \\\psi\\ of a birth-death-sampling process given a
phylogenetic tree and death rate \\\mu\\. The function uses
user-specified epidemiological constraints on the sampling rate, basic
reproduction number, and the duration of the infectious period to
improve estimability.

## Usage

``` r
estimate_bds_parameters(
  phy,
  mu,
  trim = c(0.1, 0.5),
  n_starts = 100,
  force_two_step = FALSE,
  psi_max = 7,
  r0_min = 1,
  r0_max = 5,
  infectious_period_min = NULL,
  infectious_period_max = NULL
)
```

## Arguments

- phy:

  An object of class `phylo` with branch lengths. Must be rooted and
  binary.

- mu:

  The extinction rate. This must be non-negative.

- trim:

  A numeric vector of length 2 specifying the quantiles for removing
  node times in the lineages-through-time regression (see Details).
  Default value is `c(0.10, 0.50)`.

- n_starts:

  The number of starting points for multi-start optimization. Default
  value is `100`.

- force_two_step:

  Logical. If `TRUE`, uses two-step estimates even if they violate
  constraints (see Details). Default value is `FALSE`.

- psi_max:

  Maximum sampling rate per unit time. Must be strictly positive.
  Default value is `7`.

- r0_min:

  Minimum basic reproductive number (R0 = lambda/(mu+psi)). Must be
  strictly positive. Default 1.

- r0_max:

  Maximum basic reproductive number. Must be strictly positive. Default
  5.

- infectious_period_min:

  Minimum infectious period (1/(mu+psi)). Must be strictly positive.
  Default NULL.

- infectious_period_max:

  Maximum infectious period (1/(mu+psi)). Must be strictly positive.
  Default NULL.

## Value

A list with the following attributes:

- `lambda`: Estimated speciation rate

- `psi`: Estimated sampling rate

- `mu`: Input death rate

- `r0`: Estimated basic reproductive number (\\\lambda/(\mu+\psi)\\)

- `infectious_period`: Estimated infectious period (\\1/(\mu+\psi)\\)

- `net_diversification`: Estimated net diversification rate
  (\\\lambda-\mu-\psi\\)

## Details

Estimation is carried out by first estimating the reproduction number
\\\frac{\lambda}{\mu+\psi}\\ and the net diversification rate
\\\lambda-\mu-\psi\\, then solving for the unknown \\\lambda\\ and
\\\psi\\. The net diversification rate is estimated by performing a
lineages-through-time (LTT) regression on node times that are within the
quantiles specified in the `trim` input and the reproduction number is
estimated by first computing maximum likelihood estimates (MLEs) of
\\\lambda\\ and \\\psi\\ by numerically optimizing a likelihood and then
compute the reproduction number. The stability of the MLEs is assessed
by checking if the estimate of \\R_0\\ is within 0.02 of 1.

If the estimates of \\\lambda\\ and \\\psi\\ do not satisfy all of the
constraints specified by the `psi_max`, `r0_min`, `r0_max`,
`infectious_period_min`, and `infectious_period_max` parameters, then
the MLEs of \\\lambda\\ and \\\psi\\ are returned, unless
`force_two_step=TRUE` in which case the original estimates of
\\\lambda\\ and \\\psi\\ are returned regardless of whether or not the
constraints are satisfied. If LTT regression fails then the MLEs of
\\\lambda\\ and \\\psi\\ will be returned.

The MLEs are obtained numerically using the `L-BFGS-S` algorithm (see
optim for implementation details) `n_starts` times from randomly chosen
starting points and selecting the best maximizer of all attempts. Due to
the user-imposed constraints, not all randomly generated starting points
are valid, so `100*n_starts` attempts are made to generate valid
starting points. This may result in fewer than `n_starts` optimizations
actually being performed.

If the two-step method produces estimates that violate constraints or is
numerically unstable, the function falls back to the MLE estimates from
step 1.

## See also

[`estimate_transition_rates()`](https://magpie-sfu.github.io/saasi/reference/estimate_transition_rates.md)

## Examples

``` r
data(ebola_tree)

# Use the following information about the 2013 Ebola outbreak to obtain estimates
# - Average removal rate of 5
# - Infectious periods range from 20 to 40 days
# - An upper bound on the sampling rate of 15
# - Plausible values for the basic reproduction number are between 1.5 and 3
BDS_fit <- estimate_bds_parameters(ebola_tree, mu = 5,
                                               psi_max = 15,
                                               infectious_period_min = 20/365,
                                               infectious_period_max = 40/365,
                                               r0_min = 1.5,
                                               r0_max = 3)
```
