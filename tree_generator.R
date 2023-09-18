library("ape")
library("diversitree")

get_rand_phy <- function(seed, max_taxa, λ, μ, q) {
  set.seed(3)
  pars <- c(λ, λ, μ, μ, q[[1, 2]], q[[2, 1]])
  phy <- tree.bisse(pars, max.taxa = max_taxa, x0 = 0)
  return(phy)
}
