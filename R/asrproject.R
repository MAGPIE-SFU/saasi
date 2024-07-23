#' Test package.
#'
#' @param phy TODO has tip.state
#' @param params_file TODO
#' @param q_matrix_file TODO
#' @return TODO phylo
#' @export
get_tree <- function(phy, params_file, q_matrix_file) {
  # TODO validate parameters
  params_df <- read.csv(file.path(params_file))
  q_matrix <- get_q_matrix(q_matrix_file)

  nstate <- nrow(params_df)
  nleaf_node <- phy[["Nnode"]]
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- nleaf_node * 2 + 1

  return("Hello world!")
}

#' TODO
#'
#' @param q_matrix_csv TODO
#' @return TODO
#' @noRd
get_q_matrix <- function(q_matrix_file) {
  q_matrix <- as.matrix(read.csv(file.path(q_matrix_file)))
  rownames(q_matrix) <- q_matrix[, 1]
  q_matrix <- q_matrix[, -c(1)]
  diag(q_matrix) <- NA
  class(q_matrix) <- "numeric"
  return(params_df)
}
