# c1
c1 <- function(lambda, mu, psi){
  return(abs(sqrt((lambda - mu - psi)^2 + (4 * lambda * psi))))
}

# c2
c2 <- function(lambda, mu, psi){
  return(-(lambda - mu - psi) / (c1(lambda, mu, psi)))
}

q <- function(lambda, mu, psi, t){
  d1 <- 2 * (1 - (c2(lambda, mu, psi))^2)
  d2 <- exp(-c1(lambda, mu, psi) * t) * (1 - c2(lambda, mu, psi))^2
  d3 <- exp(c1(lambda, mu, psi) * t) * (1 + c2(lambda, mu, psi))^2
  return(d1 + d2 + d3)
}

probability_density <- function(lambda, mu, psi,
                                internal_node_times, leaf_node_times) {
  # Add checks for invalid inputs
  sorted_x <- sort(internal_node_times)
  sorted_y <- sort(leaf_node_times)
  d1 <- length(internal_node_times) *log(lambda)
  d2 <- sum(-log(q(lambda, mu, psi, sorted_x)))
  d3 <- sum(log(psi) + log(q(lambda, mu, psi, sorted_y)))
  ret <- d1 + d2 + d3
  return(ret)
}

#' foo
#'@export
sim_mle_lm <- function(phy, lambda = 2, mu = 0.5, psi = 0.5,
                       method = "L-BFGS-B") {
  node_depths <- ape::node.depth.edgelength(phy)
  node_times <- max(node_depths) - node_depths
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- phy[["Nnode"]] * 2 + 1
  nleaf_node <- nnode - phy[["Nnode"]]
  leaf_node_times <- node_times[1:nleaf_node]
  internal_node_times <- node_times[(nleaf_node+1):nnode]

  # Returns negative log likelihood
  nLL <- function(lambda, mu) {
    positive_ret <- probability_density(lambda, mu, psi,
                                        internal_node_times, leaf_node_times)
    return(-positive_ret)
  }

  fit <- stats4::mle(nLL,
                     start = list(lambda = lambda, mu = mu),
                     method = method,
                     lower = c(0.1, 0.05),
                     upper = c(10, 10))
  
  ret <- unname(c(stats4::coef(fit)[1], stats4::coef(fit)[2]))
  return(ret)
}
