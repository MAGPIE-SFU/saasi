library("ape")
library("diversitree")

get_rand_phy <- function(seed, max_taxa, params_df, q_matrix) {
  set.seed(seed)

  q_vector <- as.vector(t(q_matrix))
  q_vector <- q_vector[!is.na(q_vector)]
  pars <- c(params_df$lambda,
            params_df$mu,
            q_vector)

  # TODO bad to assume root is 1?
  phy <- tree.musse(pars, max.taxa = max_taxa, x0=1)
  return(phy)
}
