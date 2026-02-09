#' Plot phylogenetic tree with saasi results
#'
#' @param tree A phylo object containing tip.state 
#' @param saasi_result Output from saasi()
#' @param colors Character vector of colors. If NULL, colors are automatically
#'   generated based on the number of states (default: NULL)
#' @param tip_cex Size of tip circles (default: 0.5)
#' @param node_cex Size of node pie charts (default: 0.2)
#' @param save_file Character. File path to save the plot (e.g. "tree.png").
#'   If NULL, plot is drawn to current device (default: NULL)
#' @param width Plot width in pixels when saving (default: 3000)
#' @param height Plot height in pixels when saving (default: 3000)
#' @param res Plot resolution in dpi when saving (default: 300)
#' @return Invisibly returns the color vector used
#' @export
#' @examples
#' # Basic plot
#' plot_saasi(tree, saasi_result)
#'
#' # Custom colors and save to file
#' plot_saasi(tree, saasi_result,
#'            colors = c("red", "blue", "green"),
#'            node_cex = 0.3,
#'            save_file = "tree.png")
#' 
plot_saasi <- function(tree,
                       saasi_result,
                       colors = NULL,
                       tip_cex = 0.5,
                       node_cex = 0.2,
                       save_file = NULL,
                       width = 3000, 
                       height = 3000,
                       res = 300) {
  n_states <- length(unique(tree$tip.state))
  
  if(is.null(colors)){
    colors <- grDevices::rainbow(n_states)
  }
  
  if(length(colors) < n_states){
    stop("Number of colors must be >= number of states")
  }
  
  # Map tip states to colors
  state_factor <- factor(tree$tip.state)
  tip_colors <- colors[as.numeric(state_factor)]
  pch_value <- 21
  
  # Open file device if saving
  if(!is.null(save_file)){
	  grDevices::png(save_file, width = width, height = height, res = res)
  }
  
  # Plotting 
  plot_phy <- ape::ladderize(tree, right = FALSE)
  plot(plot_phy, show.tip.label = FALSE)
  
  ape::tiplabels(bg = tip_colors, cex = tip_cex, adj = 0.5, pch = pch_value)
  
  ape::nodelabels(pie = saasi_result,
                  piecol = colors[1:ncol(saasi_result)],
                  cex = node_cex)
  
  graphics::legend("topleft",
         legend = levels(state_factor),
         col = colors[1:n_states],
         pch = 19,
         pt.cex = 2,
         bty = "n",
         cex = 0.8)
  
  if(!is.null(save_file)){
	  grDevices::dev.off()
  }
  invisible(colors)
}



#' Creates a data frame for saasi 
#'
#' @param states Vector of state names or number of states
#' @param lambda Speciation rate per unit time (single value or vector per state)
#' @param mu Extinction rate per unit time (single value or vector per state)
#' @param psi Sampling rate per unit time (single value or vector per state)
#' @param prior Prior probabilities (default: uniform; or vector per state), should sum up to 1.
#'
#' @return A data frame for saasi
#'
#' @export
#' @examples
#' params <- create_params_template(
#'   states = c("A", "B"),
#'   lambda = 0.5,
#'   mu = 0.01,
#'   psi = c(0.1, 0.9)  
#' )
create_params_template <- function(states,
                                   lambda = 0.5,
                                   mu = 0.01,
                                   psi = 0.5,
                                   prior = NULL) {
  
  if(is.numeric(states) && length(states) == 1){
    n_states <- states
    state_names <- seq_len(n_states)
  } 
  else{
    state_names <- states
    n_states <- length(states)
  }
  
  if(length(lambda) == 1){
    lambda <- rep(lambda, n_states)
  } 
  if(length(mu) == 1){
    mu <- rep(mu, n_states)
  } 
  if(length(psi) == 1){
    psi <- rep(psi, n_states)
  } 
  if(is.null(prior)){
    prior <- rep(1/n_states, n_states)
  } 
  
  if(length(lambda) != n_states || length(mu) != n_states || 
      length(psi) != n_states || length(prior) != n_states){
    stop("All parameter vectors must have the same length as number of states")}
  
  params_df <- data.frame(
    state = state_names,
    prior = prior,
    lambda = lambda,
    mu = mu,
    psi = psi
  )
  return(params_df)
}


#' Estimate the transition rate matrix of the discrete trait process
#' 
#' `estimate_transition_rates()` computes a maximum likelihood estimate of the transition rate matrix \eqn{Q} governing the discrete trait process.
#' The discrete traits evolve according to a Markov process with transition rate matrix \eqn{Q} such that element \eqn{q_{ij}} is the instantaneous rate of transitioning from state \eqn{i} to state \eqn{j}.
#' 
#' The `model` parameters controls the structure of \eqn{Q}. 
#' `"ER"` specifies equal rates for all transitions, \eqn{q_{ij}=q} for all traits \eqn{i,j}. 
#' `"SYM"` specifies a symmetric transition rate matrix, \eqn{q_{ij}=q_{ji}} for all traits \eqn{i,j}.
#' `"ARD"` places no structural constraints and allows all traits to be different, \eqn{q_{ij}\ne q_{kl}} for all traits \eqn{i,j,k,l}.
#' If another grouping of transition rates is desired this can be input as a numeric matrix such that each entry has a matrix with a positive integer and all entries with the same value will share a rate.
#' 
#' The `method` parameter allows the user to select a preferred method between `ace` and `simmap`, however if the preferred method fails the other will be used instead.
#' If `ace` is selected, the maximum likelihood estimator is computed using the [ape::ace()] function.
#' If `simmap` is selected, the [phytools::make.simmap()] function is used to fit \eqn{Q} to the provided phylogenetic tree.
#' 
#' @param tree A phylo object with a tip.state attribute assigning traits to all tips.
#' @param matrix_structure The form of the transition rate matrix. Can either be a numeric matrix or one of the following strings: `"ER"`, `"SYM"`, and `"ARD"` (see details). The default value is `"ER"`.
#' @param method The method used to perform the estimation. Possible values are `"ace"` or `"simmap"`.
#' @return A transition rate matrix that has been fit to the observed phylogeny. This matrix will be compatible with other `saasi` functions.
#' @export
#' @examples
#' # Load a timed phylogenetic tree for Ebola
#' data(ebola_tree)
#' 
#' # Use the simmap function to estimate the rate transition matrix
#' # - Impose a symmetric form on the matrix
#' Q <- estimate_transition_rates(tree, matrix_structure = "SYM")
#' 
#' # Use ace to estimate the rate transition matrix
#' # - Impose a structure such that transitions to or from Guinea happen at a different rate than Liberia or Sierra Leone
#' struct <- matrix(c(0, 1, 1, 1, 0, 2, 1, 2, 0), nrow=3, ncol=3)
#' Q <- estimate_transition_rates(tree, matrix_structure = struct)
estimate_transition_rates <- function(tree, 
                                      matrix_structure = "ER",
                                      method = "ace") {

  if(is.character(matrix_structure)) {
    if( !(matrix_structure %in% c("ER", "SYM", "ARD")) )
      stop("matrix_structure must be either a numeric matrix or one of the following strings: 'ER', 'SYM', or 'ARD'")
  } else if(!is.matrix(matrix_structure)) {
    stop("matrix_structure must be either a numeric matrix or one of the following strings: 'ER', 'SYM', or 'ARD'")
  }
  
  # Convert tip.state to factor if needed
  tip_states <- as.factor(tree$tip.state)
  
  q_matrix <- NULL
  primary_method <- method
  if(method == "ace"){
    fallback_method = "simmap"
  }
  else{
    fallback_method = "ace"
  }

  if(primary_method == "ace"){
    q_matrix <- try_ace(tree, tip_states, matrix_structure)
  } 
  else{
    q_matrix <- try_simmap(tree, tip_states, matrix_structure)
  }
  
  # If primary failed, try the other method
  
  if(is.null(q_matrix)){
    warning("Primary method failed")
    if(fallback_method == "ace"){
      q_matrix <- try_ace(tree, tip_states, matrix_structure)
    } 
    else{
      q_matrix <- try_simmap(tree, tip_states, matrix_structure)
    }
  }
  
  # If both failed
  if(is.null(q_matrix)){
    stop("Both ace and simmap failed to estimate transition rates. Check your tip.states, it may contains NA.")
  }

  # Replace NAs on diagonal of the matrix
  diag(q_matrix) <- NA
  diag(q_matrix) <- -rowSums(q_matrix, na.rm=T)
  return(q_matrix)
}


#' Estimate Q using ace
#' @noRd
try_ace <- function(tree, tip_states, model){
  result <- tryCatch(
    withCallingHandlers(
      {ace_result <- ape::ace(tip_states, tree, type = "discrete", model = model)
        if(is.null(ace_result$index.matrix) || is.null(ace_result$rates)){
          stop("ace failed")
        }
        ind <- ace_result$index.matrix
        ace_rates <- ace_result$rates
        q <- ind
        if(!is.null(ace_result$lik.anc)){
          rownames(q) <- colnames(ace_result$lik.anc)
          colnames(q) <- rownames(q)
        }
        for (i in 1:nrow(ind)) {
          for (j in 1:ncol(ind)) {
            q[i, j] <- ace_rates[ind[i, j]]
          }
        }
        q
      },
      warning = function(w) {
        message("ace warning: ", conditionMessage(w))
      }
    ),
    error = function(e) {
      message("ace failed: ", e$message)
      NULL
    }
  )
  return(result)
}

#' Estimate Q using simmap
#' @noRd
try_simmap <- function(tree, tip_states, model) {
  result <- tryCatch(
    withCallingHandlers(
      {states_named <- tree$tip.state
      names(states_named) <- tree$tip.label
      simmap_result <- phytools::make.simmap(tree, states_named, model = model, Q = "empirical")
      simmap_result$Q
      },
      warning = function(w) {
        message("simmap warning: ", conditionMessage(w))
      }
    ),
    error = function(e) {
      message("simmap failed: ", e$message)
      NULL
    }
  )
  return(result)
}
