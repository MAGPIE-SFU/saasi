#' Construct a rate matrix q using the output of ace 
#' 
#' This function is here because ace returns rates and an index matrix, not a matrix of rates.
#' 
#' @param aceoutput the output of ace 
#' @return A rate matrix q in the form saasi needs 
#' @export
extract_ace_q <- function(aceoutput) { 
    ind = aceoutput$index.matrix
    acerates = aceoutput$rates 
    q=ind
    rownames(q)=colnames(aceoutput$lik.anc)
    colnames(q)=rownames(q)
    for (i in 1:nrow(ind)) {
        for (j in 1:ncol(ind)) {
            q[i,j] <-  acerates[ind[i,j]]
        }
    }
    return(q)
}


#' Transition rate adjustment
#' 
#' This function adjust the estimated transition rates using `ace` based on 
#' different sampling fractions.
#' 
#' @param qij_matrix The matrix representation of the transition rates, can be
#' obtained using `extract_ace_q` function.
#' @param state The state that has higher or lower sampling rate relative to 
#' the other states. If `state` is not a character but an index (numeric value), pick
#' the index relate to the index of the matrix.
#' @param sampling_diff Sampling difference the state of interest and other states.
#' This should be a single positive value indicates the sampling difference (ratio).
#' If 0<sampling_diff<1, this implies that
#' the state of interest samples less than other states. If 
#' sampling_diff > 1, this implies that the state of interest samples more than other
#' states. 
#' @return A new transition rate matrix that is adjusted by the sampling differences.
#' @export
q_adjust <- function(qij_matrix,state,sampling_diff){
  
  if(sampling_diff <= 0){
    stop("The sampling difference should be greater than 0")
  }
  
  if(is.character(state)){
    index <- which(colnames(qij_matrix) == state)
  } else {
    index <- state # this should be a numeric value, which represent the index.
  }
  updated_qij <- qij_matrix
  
  x <- c(1, 2, 5, 10, 50, 100)
  y <- c(1, 0.75, 0.6, 0.5, 0.25, 0.2)
  model <- nls(y ~ a * x^b, start = list(a = 1, b = -0.5))
  fit_func_lower <- function(x_new) {
    predict(model, newdata = data.frame(x = x_new))
  }
  
  x <- c(1, 2, 5, 10, 50, 100)
  y <- c(1, 1.5, 2, 2.5, 4, 5)
  model2 <- nls(y ~ a * x^b, start = list(a = 1, b = -0.5))
  fit_func_higher <- function(x_new2) {
    predict(model2, newdata = data.frame(x = x_new2))
  }
  
  # these values are tested using simulations
  
  if(sampling_diff <= 1){
    updated_qij[,index] = updated_qij[,index]/fit_func_higher(1/sampling_diff)
    updated_qij[index,] = updated_qij[index,]/fit_func_lower(1/sampling_diff)
  } else {
    updated_qij[,index] = updated_qij[,index]/fit_func_lower(sampling_diff)
    updated_qij[index,] = updated_qij[index,]/fit_func_higher(sampling_diff)
  }
  return(updated_qij)
}

# Tree visualization using ggtree 
# 
# This function creates the ggtree object that visualize the tree, which includes
# tip states and inferred node state probabilities. Note, the `ggtree` package 
# should be loaded before using this function.
# 
# @param phy A `phylo` phylogenetic tree.
# @param asi Ancestral state inference using `saasi` or other methods.
# @return A phylogenetic tree that contain tip states and inferred node state probabilities.
# @export
# plot_ggtree <- function(phy,asi){
#   if (!requireNamespace("ggtree", quietly = TRUE)) {
#     stop("Package 'ggtree' is required but not installed. Please install it with BiocManager::install(`ggtree`)", call. = FALSE)
#   }
# 
# }



