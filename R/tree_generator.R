#' Random tree generation
#'
#' Get a random tree that can be passed to [asrproject::asr].
#'
#' @param seed Integer vector used in random number generation.
#' @param max_taxa Maximum number of leaf nodes to generate.
#' @param plot If `TRUE`, a plot of the returned tree is generated.
#' @return A `phylo` phylogenetic tree (ape format).
#' @inheritParams asr
#' @export
get_rand_phy <- function(seed, max_taxa, params_df, q_matrix, plot = FALSE) {
  set.seed(seed)

  q_vector <- as.vector(t(q_matrix))
  q_vector <- q_vector[!is.na(q_vector)]

  pars <- c(params_df$lambda,
            params_df$mu,
            q_vector)

  phy <- diversitree::tree.musse(pars, max.taxa = max_taxa, x0 = 1)

  if (plot) {
    nstate <- nrow(params_df)
    plot(diversitree::history.from.sim.discrete(phy, 1:nstate), phy)
    ape::tiplabels(frame = "circle", cex = 0.5)
    ape::nodelabels(frame = "circle", cex = 0.5)
  }

  return(phy)
}
