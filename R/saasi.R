#' Sampling Aware Ancestral State Inference
#'
#' Reconstructs ancestral states for internal nodes of a phylogenetic tree while 
#' accounting for heterogeneous sampling rates across states or locations.
#'
#' If `sensitivity_test=TRUE`, then a test of the sensitivity of `saasi` to the birth-death-sampling rate parameters is conducted.
#' The values of \eqn{\lambda}, \eqn{\mu}, and \eqn{\psi} are perturbed by sampling values at random from +/- 10% of the input value.
#' The ancestral state reconstruction is then run on the perturbed parameters; this is repeated 10 times.
#' To summarise sensitivity, the modal state for each internal node is first computed for each perturbed reconstruction.
#' Then the total number of differences in modal states are computed between the perturbed reconstructions and the original reconstruction.
#' The average number of these total differences across all 10 perturbed reconstructions is reported.
#'
#' @param phy An object of class `phylo` specifying the phylogenetic tree. 
#' The tree must be rooted, binary, and contain the `tip.state` attribute. 
#' Use the [check_tree_compatibility()] function to verify that `phy` is satisfies these conditions.
#' @param Q An \eqn{n\times n} transition rate matrix, where \eqn{n} is the number of states. 
#' Row and column names must correspond to the states.
#' Off-diagonal elements are the transition rates between states; diagonal elements are set such that rows sum to zero.
#' @param pars A `data.frame` with the other parameters used in the ancestral state reconstruction algorithm. 
#' Must have the following columns:
#' - `state`: a vector of state names. These should match the column and row names of `Q`.
#' - `lambda`: a vector of birth rates for the birth-death-sampling process.
#' - `mu`: a vector of death rates for the birth-death-sampling process.
#' - `psi`: a vector of sampling rates for the birth-death-sampling process.
#' 
#' It can also have `root_prior` that specifies the *a priori* distribution of the state of the root; if one is not provided a uniform distribution over the root state will be used.
#' @param sensitivity_test A boolean indicating if a sensitivity test should be conducted (see Details). The default value is `FALSE`.
#'   
#' @return A `data.frame` with a state probabilities for each internal node in `phy`. 
#'   
#' @examples
#' head(demo_pars) #contains state, lambda, mu, psi for each state
#' result <- saasi(phy = demo_tree_prepared, 
#'                 Q = demo_Q,
#'                 pars = demo_pars) 
#' head(result) 
#' @export
saasi <- function(phy, Q, pars, sensitivity_test=FALSE) {
  ## INPUT CHECKS --------------------------------------------------------------
  
  ## phy
  # Check compatibility of phy with SAASI
  if( !check_tree_compatibility(phy) ) {
    stop("`phy` is not compatible with `saasi`. Please use the `prepare_tree_for_saasi` function to reformat your `phylo` object.")
  }
  
  # Check if the tree is likely not a timed tree
  # - If the median branch length is less than 1e-3, flag it as being potentially not a timed tree
  if( quantile(phy$edge.length, probs=0.5) <= 1e-3 )
    warning("The median branch length in your tree is less than 0.001. Please double-check that your tree is a timed tree.")
  
  
  ## Q
  # Check that q_matrix is compatible, it must:
  # - be matrix with no NAs
  # - be a square matrix
  # - have rows summing to 0 (a threshold of 1e-10 is used to account small numerical issues)
  # - both rows and columns must have names
  if( any(is.na(Q)) ) {
    stop("Q contains NAs.")
  }
  if( !is.matrix(Q)) {
    stop("Q must be a matrix. Please consult the documentation for details.")
  }
  if( nrow(Q) != ncol(Q) ) {
    stop(paste0("The transition rate matrix must be square. Current input has ", nrow(Q), " rows and ", ncol(Q), " columns."))
  }
  if( any(abs(rowSums(Q)) > 1e-10) ) {
    stop("The rows of the transition rate matrix must sum to 0.")
  }
  if( is.null(rownames(Q)) || is.null(colnames(Q)) ) {
    stop("Both rows and columns of the transition rate matrix should have names corresponding to the traits.")
  }
  
  #  If Q is entered as NULL, estimate them using the default values
  if( is.null(Q) ){
    warning("No transition rate matrix provided, estimating with ace using an equal rates matrix structure. If the user desires an alternative matrix structure refer to the estimate_transition_rates function.")
    Q <- estimate_transition_rates(phy)
  }
 
  
  ## pars
  # Input should be a data.frame with no NAs
  if( any(is.na(pars)) )
    stop("pars contains NAs.")
  if( !is.data.frame(pars))
    stop("pars must be a data.frame. Please consult the documentation for details.")
  
  # Check that pars has all of the right column names
  par_names <- c("state", "lambda", "mu", "psi")
  missing_pars <- par_names[!(par_names %in% names(pars))] # vector of expected parameters that are missing in the input dataframe
  if( length(missing_pars) > 0 ){
    stop("The following columns are missing from pars: ", paste(missing_pars, sep=","))
  }
  
  # If no prior is provided, add a uniform prior
  if( !("root_prior" %in% names(pars)) ) {
    pars$root_prior <- rep(1, nrow(pars)) / nrow(pars)
  } else if( (abs(sum(pars$root_prior) - 1) > 1e-5) || min(pars$root_prior) < 0 ) { # check that root_prior is a valid probability
    stop("Please ensure that root_prior is a valid probability: all entries must be non-negative and sum to 1.")
  }
  
  # Check birth-death-sampling rates
  # - mu: all must be non-negative
  # - lambda: all must be non-negative and at least one must be positive
  # - psi: all must be non-negative and all observed states must be positive
  with(pars, {
    if( any(mu < 0) ) {
      stop("All death rates (mu) must be non-negative.")
    }
    if( all(!is.null(psi)) ) {
      if( any(psi < 0) ) {
        stop("All sampling rates (psi) must be non-negative.")
      } else if( any(sort(colnames(Q)[sum(psi > 0) > 0]) != sort(unique(phy$tip.state))) ) {
        stop("Sampling rate (psi) must be positive for all observed states.")
      }
    }
    if( all(!is.null(lambda)) ) {
      if( any(lambda < 0) ) {
        stop("All birth rates (lambda) must be non-negative.")
      } else if( sum(lambda > 0) == 0 ) {
        stop("At least one birth rate (lambda) must be positive.")
      }
    }
  })
  
  # Check that pars & Q have the same states
  # - phy can have missing states, but cannot have states not present in Q or pars
  if( any(sort(colnames(Q)) != sort(pars$state)) ){
    stop("The states specified in Q and pars do not match.")
  }
  if( any(!(unique(phy$tip.state) %in% pars$state)) ) {
    stop("There is at least one state in phy that is not in pars or Q.")
  }
  

  ## PRE-PROCESSING ------------------------------------------------------------
  # reformat input and instantiate derived quantities prior to executing main algorithm

  # Convert states to numerical values and preserve state names in separate column/attribute
  if( !is.numeric(phy$tip.state) ) {
    phy$tip.statename <- phy$tip.state
    phy$tip.state <- as.numeric(factor(phy$tip.state))
    pars$statename <- pars$state
    pars$state <- as.numeric(factor(pars$state))
    # Reorder pars and Q so that row/column corresponds to the numeric state value
    pars <- pars[order(pars$state),]
    Q <- Q[order(as.numeric(factor(rownames(Q)))), order(as.numeric(factor(colnames(Q))))]
  }
  
  nstate <- nrow(pars)
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- phy[["Nnode"]] * 2 + 1
  # Root node ID == number of leaf nodes + 1 == number of internal nodes + 2.
  root_node <- phy[["Nnode"]] + 2
  
  node_depths <- ape::node.depth.edgelength(phy)
  max_depth <- max(node_depths)
  post_order_edges <- ape::reorder.phylo(phy, "postorder")[["edge"]]
  
  topology_df <- get_topology_df(nnode,
                                 node_depths,
                                 max_depth,
                                 post_order_edges)
  
  
  ## Run SAASI -----------------------------------------------------------------
  backwards_likelihoods_list <- get_backwards_likelihoods_list(phy,
                                                               pars,
                                                               Q,
                                                               nstate,
                                                               nnode,
                                                               post_order_edges,
                                                               topology_df)
  
  state_probabilities_list <- get_state_probabilities_list(
    phy,
    pars,
    Q,
    nstate,
    nnode,
    root_node,
    post_order_edges,
    topology_df,
    backwards_likelihoods_list
  )

  
  # Perform sensitivity testing, if requested
  if( !sensitivity_test ){
    return(state_probabilities_list)
  } else{
    cat("Initial run complete. Running sensitivity analysis.","\n")
    # number of times to perturb the BDS rates and run saasi again
    n_tests <- 10
    
    # re-run saasi for n_tests perturbed parameter dataframes
    sens_res <- lapply(1:n_tests, function(x) {
      # generate perturbations by randomly generating a uniform number between 10% below and 10% above the provided value
      temp <- pars
      temp$lambda <- runif(nrow(pars), temp$lambda*0.9, temp$lambda*1.1)
      temp$mu <- runif(nrow(pars), temp$mu*0.9, temp$mu*1.1)
      temp$psi <- runif(nrow(pars), temp$psi*0.9, temp$psi*1.1)
      
      # re-compute the backwards likelihoods
      bw_likelihood <- get_backwards_likelihoods_list(phy,
                                                      temp,
                                                      Q,
                                                      nstate,
                                                      nnode,
                                                      post_order_edges,
                                                      topology_df)
      # re-compute and return the state probabilities for the internal nodes
      get_state_probabilities_list(phy,
                                   temp,
                                   Q,
                                   nstate,
                                   nnode,
                                   root_node,
                                   post_order_edges,
                                   topology_df,
                                   bw_likelihood)
    })
    
    # summarise sensitivity by computing the modal states for 
    #  the input parameter values (base) and for the perturbed results
    #  then count the total number of nodes that differ
    base_modes <- modal_states(state_probabilities_list)
    sens_modes <- lapply(sens_res, modal_states)
    sens_total_diff <- sapply(sens_modes, function(x) sum(x != base_modes))
    
    cat("Sensitivity analysis: Input parameters were perturbed", n_tests, "times and an average of",
        mean(sens_total_diff), "of all", nnode, "nodes had a different modal state.")
    
    # compute the sensitivity summaries
    return(state_probabilities_list)
  }
}

#' Get data frame representation of tree topology.
#'
#' @return Data frame showing the parents, children, and depth of nodes.
#' @noRd
get_topology_df <- function(nnode, node_depths, max_depth, post_order_edges) {
  topology_df <- data.frame(id = seq_len(nnode),
                            t_root = NA,
                            left = NA,
                            right = NA,
                            parent = NA)
  
  # Populate topology df with time distances from present day
  invisible(lapply(seq_len(length(node_depths)), function(i) {
    topology_df$t_root[topology_df$id == i] <<- max_depth - node_depths[i]
  }))
  
  # Populate topology df with parent/children connections
  invisible(lapply(seq(1, length(post_order_edges[, 1]), 2), function(i) {
    node <- post_order_edges[i, 1]
    left <- post_order_edges[i, 2]
    right <- post_order_edges[i + 1, 2]
    
    topology_df$parent[topology_df$id == left] <<- node
    topology_df$parent[topology_df$id == right] <<- node
    topology_df$left[topology_df$id == node] <<- left
    topology_df$right[topology_df$id == node] <<- right
  }))
  
  return(topology_df)
}

#' Get list of likelihoods calculated during backwards traversal of ancestral
#' state reconstruction algorithm.
#'
#' @return List of state likelihoods used in backwards time equations.
#' `list[[x]][[y]]` is the likelihood for state y in node x.
#' @noRd
get_backwards_likelihoods_list <- function(phy,
                                           params_df,
                                           Q,
                                           nstate,
                                           nnode,
                                           post_order_edges,
                                           topology_df) {
  backwards_likelihoods_list <- rep(list(rep(0, nstate)), nnode)
  
  # Populate leaf node likelihoods for backwards time equations
  invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
    state <- phy[["tip.state"]][[i]]
    state_freq <- params_df$psi[params_df$state == state]
    backwards_likelihoods_list[[i]][[state]] <<- state_freq
  }))
  
  # Populate internal node likelihoods for backwards time equations
  invisible(lapply(seq(1, length(post_order_edges[, 1]), 2), function(i) {
    node <- post_order_edges[[i, 1]]
    left <- post_order_edges[[i, 2]]
    right <- post_order_edges[[i + 1, 2]]
    
    tf <- topology_df$t_root[topology_df$id == node]
    left_t0 <- topology_df$t_root[topology_df$id == left]
    right_t0 <- topology_df$t_root[topology_df$id == right]
    
    left_likelihoods <- backwards_likelihoods_list[[left]]
    right_likelihoods <- backwards_likelihoods_list[[right]]
    likelihoods <- get_backwards_likelihoods(abs(left_likelihoods),
                                             abs(right_likelihoods),
                                             left_t0, right_t0, tf,
                                             params_df, Q)
    #abs_likelihoods <- abs(likelihoods)
    #norm_likelihoods <- abs_likelihoods / sum(abs_likelihoods)
    backwards_likelihoods_list[[node]] <<- likelihoods # note changed by CC 
  }))
  
  return(backwards_likelihoods_list)
}

#' Get list of probabilities ultimately generated by ancestral state
#' reconstruction algorithm.
#'
#' @return List of ancestral state probabilities.
#' `list[[x]][[y]]` is the probability of state y in node x.
#' @noRd
get_state_probabilities_list <- function(phy,
                                         params_df,
                                         Q,
                                         nstate,
                                         nnode,
                                         root_node,
                                         post_order_edges,
                                         topology_df,
                                         backwards_likelihoods_list) {
  state_probabilities_list <- rep(list(rep(0, nstate)), nnode)
  
  # Populate leaf node state probabilities
  invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
    state <- phy[["tip.state"]][[i]]
    state_probabilities_list[[i]][[state]] <<- 1
  }))
  
  state_probabilities_list[[root_node]] <- (
    backwards_likelihoods_list[[root_node]] * params_df$root_prior
    / (sum(backwards_likelihoods_list[[root_node]] * params_df$root_prior))
  )
  
  # Populate internal node state probabilities
  invisible(lapply(seq(length(post_order_edges[, 1]), 1, -2), function(i) {
    node <- post_order_edges[[i, 1]]
    if (node == root_node) {
      return()
    }
    
    parent <- topology_df$parent[topology_df$id == node]
    parent_state_probabilities <- (
      state_probabilities_list[[parent]] * backwards_likelihoods_list[[node]]
    )
    abs_parent_state_probabilities <- abs(parent_state_probabilities)
    norm_probabilities <- (
      params_df$lambda*abs_parent_state_probabilities / sum(params_df$lambda*abs_parent_state_probabilities)
    )
    
    t0 <- topology_df$t_root[topology_df$id == parent]
    tf <- topology_df$t_root[topology_df$id == node]

    likelihoods <- get_forwards_likelihoods(norm_probabilities,
                                            t0, tf,
                                            params_df, Q)
    #abs_likelihoods <- abs(likelihoods)
    state_probabilities_list[[node]] <<- (
      backwards_likelihoods_list[[node]] * likelihoods
      / sum(backwards_likelihoods_list[[node]] * likelihoods)
    )
  }))
  internal_node_ids <- root_node:nnode
  
  node_prob <- matrix(unlist(state_probabilities_list[internal_node_ids]),ncol=nstate,byrow = TRUE)
  node_lik <- as.data.frame(node_prob,row.names = internal_node_ids)
  colnames(node_lik) = c(1:nstate)
  if(!is.null(phy$tip.statename)){
      colnames(node_lik) = levels(factor(phy$tip.statename))
  }
  return(node_lik)
}

#' Generates state likelihoods for a node from the state likelihoods of the node
#' children.
#'
#' @return Vector of state likelihoods.
#' @noRd
get_backwards_likelihoods <- function(left_likelihoods, right_likelihoods,
                                      left_t0, right_t0, tf,
                                      params_df, q_matrix) {
  left_sol <- backwards_likelihoods_helper(left_likelihoods,
                                           left_t0, tf,
                                           params_df, q_matrix)
  right_sol <- backwards_likelihoods_helper(right_likelihoods,
                                            right_t0, tf,
                                            params_df, q_matrix)
  likelihoods <- 2*params_df$lambda * left_sol * right_sol
  return(likelihoods / sum(likelihoods))
}

#' State likelihoods from one child node only.
#'
#' @return Vector of state likelihoods.
#' @noRd
backwards_likelihoods_helper <- function(child_likelihoods,
                                         t0,
                                         tf,
                                         params_df,
                                         q_matrix) {
  # Number of states == n
  nstate <- nrow(params_df)
  
  func <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      # States 1...n
      dD_equations_list <- lapply(seq_len(nstate), function(i) {
        # With i =/= j:
        # * psi[-i] is psi[j]
        # * q[i,][-i] is q[i, j]
        # * y[i + nstate] is Ei
        # * y[1:nstate][-i] is Dj
        # See derivation sxn in https://www.biorxiv.org/content/10.1101/2025.05.20.655151v1 for details
        return(
          -(lambda[i] + mu[i] + psi[i] + sum(q[i, ][-i])) * y[i]
          + 2 * lambda[i] * y[i + nstate] * y[i]
          + sum(q[i, ][-i] * y[1:nstate][-i])
        )
      })
      
      # States 1...n
      dE_equations_list <- lapply(seq_len(nstate), function(i) {
        # With i =/= j:
        # * psi[-i] is psi[j]
        # * q[i,][-i] is q[i, j]
        # * y[i + nstate] is Ei
        # * y[nstate + 1:nstate][-i] is Ej
        # See derivation sxn in https://www.biorxiv.org/content/10.1101/2025.05.20.655151v1 for details
        return(
          mu[i] - (lambda[i] + mu[i] + psi[i] + sum(q[i, ][-i]))
          * y[i + nstate] + lambda[i] * y[i + nstate]^2
          + sum(q[i, ][-i] * y[nstate + 1:nstate][-i])
        )
      })
      
      return(list(c(dD_equations_list, dE_equations_list)))
    })
  }
  
  
  # D1...Dn are child likelihoods, and E1...En are 1
  y <- c(child_likelihoods, rep(1, nstate))
  # Need to explicitly name index or events_df does not work
  names(y) <- seq_len(nstate * 2)
  
  times <- seq(0, tf, by = tf / 100)
  # brj: resolve warnings about t0 not in times
  times <- sort(unique(c(times, t0)))
  
  parms <- list(lambda = params_df$lambda,
                mu = params_df$mu,
                psi = params_df$psi,
                q = q_matrix,
                nstate = nstate)
  
  # Force D1...Dn at t0 to be same as children
  events_df <- data.frame(var = seq_len(nstate),
                          time = rep(t0),
                          value = child_likelihoods,
                          method = rep("replace", nstate))
  
  sol <- deSolve::ode(y, times, func, parms, events = list(data = events_df))
  
  ret <- utils::tail(sol, n = 1)[1 + 1:nstate]
  if (any(is.nan(ret))) {
    # if the value is nan, assign the likelihood to the previous likelihood,
    # because the value is too close so that the ode cannot tell the difference.
    return(child_likelihoods)
  }
  return(ret)
}


#' Generates state likelihoods for a node from the state likelihoods of the node
#' parent.
#'
#' @return Vector of state likelihoods.
#' @noRd
get_forwards_likelihoods <- function(parent_state_probabilities, t0, tf,
                                     params_df, q_matrix) {
  # Number of states == n
  nstate <- nrow(params_df)
  
  func <- function(x, y, parms) {
    with(as.list(c(y, parms)), {
      # States 1...n
      dD_equations_list <- lapply(seq_len(nstate), function(i) {
        # NOTE: for the anagenetic speciation, the only way that a
        # state could change is due to transistion events
        # see saasi.clad for cladogenetic change equation
        return(
          -(sum(q[i, ][-i]) * y[i])
          + sum(t(q)[i, ][-i] * y[1:nstate][-i])
        )
      })
      
      # States 1...n
      dE_equations_list <- lapply(seq_len(nstate), function(i) {
        # With i =/= j:
        # * psi[-i] is psi[j]
        # * q[i,][-i] is q[i, j]
        # * y[i + nstate] is Ei
        # * y[nstate + 1:nstate][-i] is Ej
        # See derivation sxn in doi:10.1111/j.2041-210X.2012.00234.x for details
        return(
          mu[i] - (lambda[i] + mu[i] + psi[i] + sum(q[i, ][-i]))
          * y[i + nstate] + lambda[i] * y[i + nstate]^2
          + sum(q[i][-i] * y[nstate + 1:nstate][-i])
        )
      })
      
      return(list(c(dD_equations_list, dE_equations_list)))
    })
  }
  
  # D1...Dn are parent state probabilities, and E1...En are 1
  y <- c(parent_state_probabilities, rep(1, nstate))
  
  # Increment time in the positive direction because otherwise the ode solver
  # can run into errors with negative numbers being smaller than machine min.
  times <- seq(0, t0, by = t0 / 100)
  parms <- list(lambda = params_df$lambda,
                mu = params_df$mu,
                psi = params_df$psi,
                q = q_matrix,
                nstate = nstate)
  
  # brj: can't reproduce warning
  sol <- deSolve::ode(y, times, func, parms, method = "ode45", rtol = 1e20)
  
  # Closest index to tf
  closest_index <- which.min(abs(sol[, 1] - tf))
  likelihoods <- unname(sol[closest_index, 1 + 1:nstate])
  
  return(likelihoods)
}


#' Compute modal state for all nodes
#' 
#' `modal_states()` is a post-processing function that extracts the
#' node with the greatest probability at each node.
#' 
#' @param node_probs A list of the same format as output by `saasi()`.
#' @return A vector corresponding to the nodes where the elements of the vector
#' are the states with the highest probability at each node.
#' @noRd
modal_states <- function(node_probs) {
  n_nodes <- length(node_probs[[1]])
  temp <- sapply(1:n_nodes, function(i) {
    which.max(sapply(node_probs, function(x) x[i]))
  })
  return(temp)
}