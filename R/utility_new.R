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


#' Estimates transition rates using ace or simmap.
#' Default using ace, if fails use simmap instead
#'
#' @param tree A phylo object with tip.state attribute
#' @param model Transition model: "ER" (equal rates), "SYM" (symmetric), 
#'   "ARD" (all rates different), or "custom" (default: "ER")
#' @param custom_q Custom Q matrix (only used if model = "custom")
#' @param method Which method to use.
#' @return A transition rate matrix (Q matrix) ready for SAASI
#' @export
#' @examples
#' # Symmetric model
#' data(ebola_tree)
#' q_matrix <- estimate_transition_rates(ebola_tree, model = "SYM")
estimate_transition_rates <- function(tree, 
                                      model = "ER", 
                                      custom_q = NULL,
                                      method = "ace") {
    if(model == "custom"){
    if(is.null(custom_q)){
      stop("Must provide custom_q if model = 'custom'")
    }
  }
  
  if(!model %in% c("ER", "SYM", "ARD")){
    stop("model must be one of the following: 'ER', 'SYM', 'ARD', or 'custom'")
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
    q_matrix <- try_ace(tree, tip_states, model)
  } 
  else{
    q_matrix <- try_simmap(tree, tip_states, model)
  }
  
  # If primary failed, try the other method
  
  if(is.null(q_matrix)){
    warning("Primary method failed")
    if(fallback_method == "ace"){
      q_matrix <- try_ace(tree, tip_states, model)
    } 
    else{
      q_matrix <- try_simmap(tree, tip_states, model)
    }
  }
  
  # If both failed
  if(is.null(q_matrix)){
    stop("Both ace and simmap failed to estimate transition rates. Check your tip.states, it may contains NA.")
  }
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
