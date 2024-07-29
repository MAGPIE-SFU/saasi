#' TODO
#'
#' two_sampled_states_phy <-
#' get_rand_phy(3, 20, "two_sampled_states_params.csv",
#' "two_sampled_states_q.csv", plot = TRUE)
#'
#' @param seed TODO
#' @param max_taxa TODO
#' @param params_df TODO
#' @param q_matrix_file TODO
#' @param plot TODO
#' @return TODO
#' @export
get_rand_phy <- function(seed, max_taxa, params_df, q_matrix_file,
                         plot = FALSE) {
  set.seed(seed)

  q_matrix <- get_q_matrix(q_matrix_file)
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
