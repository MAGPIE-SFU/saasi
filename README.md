
<!-- README.md is generated from README.Rmd. Please edit that file -->

\# saasi

<!-- badges: start --> <!-- badges: end -->

Saasi is an ancestral state reconstruction method that accounts for
variation in sampling rates among locations or traits.

\## Installation

You can install the development version of saasi from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")

remotes::install_github("MAGPIE-SFU/saasi")
```

For this particular branch, use:

``` r
# install.packages("remotes")

remotes::install_github("MAGPIE-SFU/saasi", ref = "saasi-maintaining")
```

Before running this demo, please make sure the following packages are
installed:

``` r
knitr::opts_chunk$set(echo = TRUE)
library(saasi)
library(diversitree)
#> Loading required package: ape
library(tidytree)
#> If you use the ggtree package suite in published research, please cite
#> the appropriate paper(s):
#> 
#> G Yu. Data Integration, Manipulation and Visualization of Phylogenetic
#> Trees (1st ed.). Chapman and Hall/CRC. 2022. ISBN: 9781032233574
#> 
#> Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
#> ggtree: an R package for visualization and annotation of phylogenetic
#> trees with their covariates and other associated data. Methods in
#> Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628
#> 
#> Attaching package: 'tidytree'
#> The following objects are masked from 'package:ape':
#> 
#>     drop.tip, keep.tip
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(ape)
library(phytools)
#> Loading required package: maps
library(readr)
```

## saasi package

This is a demo showing how to use the saasi package.

saasi requires the following argument: (1) a phylogenetic tree (class
`phylo`). The tree should satisfy the following properties: rooted,
binary, branch length in units of time and positive, has tree\$tip.state
and has no missing states. Check `check_tree_compatibility` function for
details. In addition, `prepare_tree_for_saasi` function helps the user
to prepare a tree that is compatible with saasi. (2) speciation,
extinction and sampling rates (class `data.frame`). These rates can be
estimated using the function `estimate_bds_parameters`. (3) a transition
rate matrix (class `matrix`). The output will be a data frame that
contains the probabilities of each state for each internal node of the
phylogenetic tree. The Q matrix can be estimated using the function
`estimate_transition_rates`.

## Example: Ebola 2013-2016 West African Ebola Epidemic

For the next example, download the data from Nextstrain:
<https://nextstrain.org/ebola/ebov-2013?c=country>

``` r
# Read the tree file and metadata (replace with your own directory)
tree <- read.tree("testing_data_nextstrain/nextstrain_ebola_ebov-2013_timetree.nwk")
metadata <- read_tsv("testing_data_nextstrain/nextstrain_ebola_ebov-2013_metadata.tsv")
#> Rows: 1493 Columns: 8
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr  (7): strain, country, division, author, author__url, accession, accessi...
#> date (1): date
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

# Create a data.frame contains strains and states. 
# The first column should match tree$tip.label, 
# the second column is the region or country that you are interested in.
tip_data <- data.frame(
  tip_label = metadata$strain,
  state = metadata$country
)

# Check compatibility
check_tree_compatibility(tree)
#> Tree is not compatible with SAASI
#> Missing internal nodes
#> Unary nodes present. Suppress with ape::collapse.branches() or re-run multi2di()
#> Polytomies present. Resolve with multi2di(tree, tol = 1e-8)
#> Zero-length branches present. Fix with: tree$edge.length[tree$edge.length <= 0] <- 1e-5
#> Negative branches present. Fix with: tree$edge.length[tree$edge.length <= 0] <- 1e-5
#> NA tip states present. Remove with drop_tips_by_state(tree, NA)
#> No tip states. Attach with attach_tip_states() or prepare_tree_for_saasi()
#> [1] FALSE

# The tree is not compatible. Specifically, the tree is nonbinary, it does
# not contains tip states, contains negative/zero branch length, and contains
# polytomy.
# To fix it, run the following command:

ebola_tree <- prepare_tree_for_saasi(tree,tip_data)

# Check compatibility again
check_tree_compatibility(ebola_tree)
#> Tree is compatible with SAASI
#> [1] TRUE

# The tree is compatible with saasi.

# Estimate the transition rate matrix Q
Q <- estimate_transition_rates(ebola_tree,method = 'simmap',model = 'SYM')
#> make.simmap is sampling character histories conditioned on
#> the transition matrix
#> 
#> Q =
#>                  Guinea    Liberia Sierra Leone
#> Guinea       -0.4718132  0.2604072    0.2114061
#> Liberia       0.2604072 -0.3697032    0.1092960
#> Sierra Leone  0.2114061  0.1092960   -0.3207021
#> (estimated using likelihood);
#> and (mean) root node prior probabilities
#> pi =
#>       Guinea      Liberia Sierra Leone 
#>    0.3333333    0.3333333    0.3333333
#> Done.

# Estimate the rates
# For ebola, we assume the total infectious period 1/(mu+psi) is between 
# 20 - 40 days (20/365 - 40/365 years)
# The bound for mu+psi is [9.125,18.25]
# We further assume psi > mu (user can also assume mu > psi)
# Set mu = 5
rates <- estimate_bds_parameters(
    ebola_tree,
    mu = 5, 
    r0_max = 3, 
    r0_min = 1.5,
    psi_max = 15,
    infectious_period_min = 20/365, # convert days to years
    infectious_period_max = 40/365, # convert days to years
    n_starts = 100)

# Setting up parameters
pars1 <- create_params_template(colnames(Q),lambda = rates$lambda,mu = rates$mu,psi = rates$psi)

# Run saasi analysis
saasi_ebola <- saasi(ebola_tree,pars1,Q)

# Plot and save the result, set save_file = NULL to not save file
p1 <- plot_saasi(ebola_tree,saasi_ebola,save_file = "ebola_equal_psi.png")

# Try a different set up, we assume Liberia samples less than other countries
pars2 <- create_params_template(colnames(Q),lambda = rates$lambda,mu = rates$mu,psi = c(rates$psi,rates$psi/2,rates$psi))

# Rerun analysis and save the result, set save_file = NULL to not save file
saasi_ebola2 <- saasi(ebola_tree,pars2,Q)

p2 <- plot_saasi(ebola_tree,saasi_ebola2,save_file = "ebola_unequal_psi.png")
```

``` r
# Notice that the transition rate matrix Q is relatively small compare to the other rates.
# This is due to lack of transition events between states.
# For sanity check, randomly shuffle the tip.state

testing_ebola <- ebola_tree
set.seed(123) 
testing_ebola$tip.state <- sample(testing_ebola$tip.state)

Q <- estimate_transition_rates(testing_ebola,method = 'simmap',model = 'SYM')
#> make.simmap is sampling character histories conditioned on
#> the transition matrix
#> 
#> Q =
#>                 Guinea   Liberia Sierra Leone
#> Guinea       -342.1827  171.0913     171.0913
#> Liberia       171.0913 -342.1827     171.0913
#> Sierra Leone  171.0913  171.0913    -342.1827
#> (estimated using likelihood);
#> and (mean) root node prior probabilities
#> pi =
#>       Guinea      Liberia Sierra Leone 
#>    0.3333333    0.3333333    0.3333333
#> Done.

# Now Q gets extremely large
```
