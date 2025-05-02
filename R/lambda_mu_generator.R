#' Helper fn for `probability_density`.
#'
#' @noRd
c1 <- function(lambda, mu, psi) {
  return(abs(sqrt((lambda - mu - psi)^2 + (4 * lambda * psi))))
}

#' Helper fn for `probability_density`.
#'
#' @noRd
c2 <- function(lambda, mu, psi) {
  return(-(lambda - mu - psi) / (c1(lambda, mu, psi)))
}

#' Helper fn for `probability_density`.
#'
#' @noRd
q <- function(lambda, mu, psi, t) {
  d1 <- 2 * (1 - (c2(lambda, mu, psi))^2)
  d2 <- exp(-c1(lambda, mu, psi) * t) * (1 - c2(lambda, mu, psi))^2
  d3 <- exp(c1(lambda, mu, psi) * t) * (1 + c2(lambda, mu, psi))^2
  return(d1 + d2 + d3)
}

#' Helper fn for `mle_lm`.
#'
#' @noRd
probability_density <- function(lambda, mu, psi,
                                internal_node_times, leaf_node_times) {
  # Add checks for invalid inputs
  sorted_x <- sort(internal_node_times)
  sorted_y <- sort(leaf_node_times)
  d1 <- length(internal_node_times) * log(lambda)
  d2 <- sum(-log(q(lambda, mu, psi, sorted_x)))
  d3 <- sum(log(psi) + log(q(lambda, mu, psi, sorted_y)))
  ret <- d1 + d2 + d3
  return(ret)
}

#' Estimate speciation/extinction rates for a tree
#'
#' This is done by method of maximum likelihood. Mostly implemented by
#' \href{https://github.com/yexuansong}{@yexuansong}, from methods described in
#' \href{https://doi.org/10.1093/molbev/msr217}{Stadler et al. (2012)}.
#'
#' @param phy A `phylo` phylogenetic tree (`ape` format).
#' @param lambda An initial "guess" for speciation, used in subsequent formulae.
#' @param mu An initial "guess" for extinction, used in subsequent formulae.
#' @param psi Sampling rate.
#' @param method See `method` parameter in \link{mle}.
#' @param lower A two-element vector containing the lower bound of the speciation 
#' and extinction values (geater than 0), the default is c(0.001,0.001).
#' @param upper A two-element vector containing the upper bound of the speciation 
#' and extinction values, should be similar to the sampling rate, still need to find
#' a better justification for the upcoming version.
#' @return A two-element vector containing speciation and extinction values
#' estimated from `phy`.
#'@export
mle_lm <- function(phy, lambda, mu, psi, method = "L-BFGS-B", lower = c(0.001,0.001), upper) {
  node_depths <- ape::node.depth.edgelength(phy)
  node_times <- max(node_depths) - node_depths
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- phy[["Nnode"]] * 2 + 1
  nleaf_node <- nnode - phy[["Nnode"]]
  leaf_node_times <- node_times[1:nleaf_node]
  internal_node_times <- node_times[(nleaf_node + 1):nnode]

  negative_log_likelihood <- function(lambda, mu) {
    positive_ret <- probability_density(lambda, mu, psi,
                                        internal_node_times, leaf_node_times)
    return(-positive_ret)
  }

  fit <- stats4::mle(negative_log_likelihood,
                     start = list(lambda = lambda, mu = mu),
                     method = method,
                     lower = lower,
                     upper = upper)

  ret <- unname(c(stats4::coef(fit)[1], stats4::coef(fit)[2]))
  return(ret)
}
