#' Plot a phylogenetic tree with SAASI tip annotations
#'
#' `plot_saasi()` plots a phylogenetic tree annoted with [saasi()] output. Tips are annotated with their states.
#' Internal nodes are annotated with a pie chart showing the marginal distribution over states for each node.
#' The [ape::plot.phylo()] function is used for plotting.
#'
#' @param tree An object of class `phylo` with `tip.state` attribute. 
#' @param saasi_result Output of `saasi()` containing internal node annotations.
#' @param colours Optional character vector specifying the node colour for each state. By default, colours are automatically
#'   generated based on the number of states.
#' @param tip_cex Size of the tips. Default value is 0.5.
#' @param node_cex Size of node pie charts. Default value is 0.2.
#' @param save_file Optional string specifying the file path where the plot will be saved.
#'   By default the plot is displayed on the current device.
#' @param width Width of the plot in pixels. Only relevant if saving to file. Default value is 3000 pixels.
#' @param height Height of the plot in pixels. Only relevant if saving to file. Default value is 3000 pixels.
#' @param res Plot resolution in dpi. Only relevant if saving to file. Default value is 300 dpi.
#' @return If a save file is specified, the plot will be saved as directed. 
#' If not, then the plot will be displayed on the users device.
#' The function returns the output of `plot.phylo()`, a list containing plot specifications.
#' @export
#' @examples
#' # Run SAASI
#' saasi_res <- saasi(demo_tree_prepared, demo_Q, demo_pars)
#' 
#' # Plot the results using the default settings
#' plot_saasi(demo_tree_prepared, saasi_res)
#'
#' \dontrun{
#' # Use custom colours and save the result to tree.png
#' plot_saasi(tree, saasi_result,
#'            colours = c("orange", "darkblue", "pink"),
#'            node_cex = 0.3,
#'            save_file = "tree.png")
#' }
#' 
plot_saasi <- function(tree,
                       saasi_result,
                       colours = NULL,
                       tip_cex = 0.5,
                       node_cex = 0.2,
                       save_file = NULL,
                       width = 3000, 
                       height = 3000,
                       res = 300) {
  n_states <- length(unique(tree$tip.state))
  
  if(is.null(colours)){
    colours <- grDevices::rainbow(n_states)
  }
  
  if(length(colours) < n_states){
    stop("You must specify at least as many colours as there are states.")
  }
  
  # Map tip states to colours
  state_factor <- factor(tree$tip.state)
  tip_colours <- colours[as.numeric(state_factor)]
  pch_value <- 21
  
  # Open file device if saving
  if(!is.null(save_file)){
    grDevices::png(save_file, width = width, height = height, res = res)
  }
  
  # Plotting 
  plot_phy <- ape::ladderize(tree, right = FALSE)
  res <- plot(plot_phy, show.tip.label = FALSE)
  
  ape::tiplabels(bg = tip_colours, cex = tip_cex, adj = 0.5, pch = pch_value)
  
  if(!is.null(saasi_result)){
  ape::nodelabels(pie = saasi_result,
                  piecol = colours[1:ncol(saasi_result)],
                  cex = node_cex)
  }
  
  graphics::legend("topleft",
                   legend = levels(state_factor),
                   col = colours[1:n_states],
                   pch = 19,
                   pt.cex = 2,
                   bty = "n",
                   cex = 0.8)
  
  if(!is.null(save_file)){
    grDevices::dev.off()
  }
  return(invisible(res))
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
#' The `method` parameter allows the user to select a preferred method between `ace` and `fitMk`, however if the preferred method fails the other will be used instead.
#' If `ace` is selected, the maximum likelihood estimator is computed using the [ape::ace()] function.
#' If `fitMk` is selected, the [phytools::fitMk()] function is used to fit \eqn{Q} to the provided phylogenetic tree.
#' 
#' @param tree A phylo object with a tip.state attribute assigning traits to all tips.
#' @param matrix_structure The form of the transition rate matrix. Can either be a numeric matrix or one of the following strings: `"ER"`, `"SYM"`, and `"ARD"` (see details). The default value is `"ER"`.
#' @param method The method used to perform the estimation. Possible values are `"ace"` or `"fitMk"`.
#' @return A transition rate matrix that has been fit to the observed phylogeny. This matrix will be compatible with other `saasi` functions.
#' @export
#' @examples
#' # Load a timed phylogenetic tree for Ebola
#' data(ebola_tree)
#' 
#' # Use the simmap function to estimate the rate transition matrix
#' # - Impose a symmetric form on the matrix
#' Q <- estimate_transition_rates(ebola_tree, matrix_structure = "SYM")
#' 
#' # Use ace to estimate the rate transition matrix
#' # - Impose a structure such that transitions to or from Guinea happen at a different rate than Liberia or Sierra Leone
#' struct <- matrix(c(0, 1, 1, 1, 0, 2, 1, 2, 0), nrow=3, ncol=3)
#' Q <- estimate_transition_rates(ebola_tree, matrix_structure = struct)
estimate_transition_rates <- function(tree, 
                                      matrix_structure = "ER",
                                      method = "ace") {
  
  # check that matrix_structure is a valid input
  if(is.character(matrix_structure)) {
    if( !(matrix_structure %in% c("ER", "SYM", "ARD")) )
      stop("matrix_structure must be either a numeric matrix or one of the following strings: 'ER', 'SYM', or 'ARD'")
  } else if(!is.matrix(matrix_structure)) {
    stop("matrix_structure must be either a numeric matrix or one of the following strings: 'ER', 'SYM', or 'ARD'")
  }
  
  # check that tree has class phylo and has valid tip states
  if( class(tree) != "phylo" )
    stop("tree must be of class \"phylo\".")
  if( !("tip.state" %in% names(tree)) )
    stop("tree must contain a tip.state attribute containing the tip states.")
  if( sum(is.na(tree$tip.state)) > 0 )
    stop("Tip states contain at least 1 NA value, please remove NA states using the drop_tips_by_state() function.")
  # Convert tip.state to factor if needed
  tip_states <- as.factor(tree$tip.state)
  
  q_matrix <- NULL
  primary_method <- method
  if(method == "ace"){
    fallback_method = "fitMk"
  }
  else if(method == "fitMk"){
    fallback_method = "ace"
  } else{
    stop("User must select either \"ace\" or \"fitMk\" as the primary estimation method.")
  }
  
  if(primary_method == "ace"){
    q_matrix <- try_ace(tree, tip_states, matrix_structure)
  } 
  else{
    q_matrix <- try_fitMk(tree, tip_states, matrix_structure)
  }
  
  # If primary failed, try the other method
  
  if(is.null(q_matrix)){
    warning(paste0(primary_method, " failed, using ", fallback_method, " instead."))
    if(fallback_method == "ace"){
      q_matrix <- try_ace(tree, tip_states, matrix_structure)
    } 
    else{
      q_matrix <- try_fitMk(tree, tip_states, matrix_structure)
    }
  }
  
  # If both failed
  if(is.null(q_matrix)){
    stop("Both ace and fitMk failed to estimate transition rates. Check your tip.states, it may contains NA.")
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
try_fitMk <- function(tree, tip_states, model) {
  result <- tryCatch(
    withCallingHandlers(
      {states_named <- tree$tip.state
      names(states_named) <- tree$tip.label
      est <- phytools::fitMk(tree, states_named, model=model)
      q_matrix <- matrix(est$rates[est$index.matrix], nrow=nrow(est$index.matrix))
      rownames(q_matrix) <- est$states
      colnames(q_matrix) <- rownames(q_matrix)
      q_matrix 
      },
      warning = function(w) {
        message("fitMk warning: ", conditionMessage(w))
      }
    ),
    error = function(e) {
      message("fitMk failed: ", e$message)
      NULL
    }
  )
  return(result)
}