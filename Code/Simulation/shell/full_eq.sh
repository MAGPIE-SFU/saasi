#!/bin/bash

# Simple ASR methods comparison
N_SIMS=${1:-10}
NTAXA=${2:-100}
MAX_TIME=${3:-100}
OUTPUT_DIR=${4:-"simple_asr"}

mkdir -p "$OUTPUT_DIR"

# Create R script
cat > "$OUTPUT_DIR/run_asr.R" << 'EOF'
library(saasi)
library(diversitree)
library(ape)
library(phytools)
library(tidytree)

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
        # See derivation sxn in doi:10.1111/j.2041-210X.2012.00234.x for details
        return(
          -(lambda[i] + mu[i] + psi[i] + sum(q[i, ][-i])) * y[i]
          + 2 * lambda[i] * y[i + nstate] * y[i]
          #+ sum(q[i, ][-i] * y[1:nstate][-i])
          +  sum(t(q)[i, ][-i] * y[1:nstate][-i])
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
          #+ sum(q[i][-i] * y[nstate + 1:nstate][-i])
          + sum(t(q)[i, ][-i] * y[nstate + 1:nstate][-i])
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
  parms <- list(lambda = params_df$lambda,
                mu = params_df$mu,
                psi = params_df$psi,
                q = q_matrix,
                nstate = nstate)

  # Force D1...Dn at t0 to be same as children
  # events_df <- data.frame(var = seq_len(nstate),
  #                         time = rep(t0),
  #                         value = child_likelihoods,
  #                         method = rep("replace", nstate))
  # 
  # # Suppress warnings about t0 not in times
  # suppressWarnings(
  #   sol <- deSolve::ode(y, times, func, parms, events = list(data = events_df))
  # )
  if (t0 > 0) {
    events_df <- data.frame(var = seq_len(nstate),
                            time = rep(t0),
                            value = child_likelihoods,
                            method = rep("replace", nstate))
    suppressWarnings(
      sol <- deSolve::ode(y, times, func, parms, events = list(data = events_df))
    )
  } else {
    # No events needed when t0=0, initial conditions already set
    sol <- deSolve::ode(y, times, func, parms)
  }
  #print(sol)

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

  # Suppress warnings about initial conditions guessed as 0
  suppressWarnings(
    # Run opposite directions because of positively increasing x. Should not
    # affect result.
    sol <- deSolve::ode(y, times, func, parms, method = "ode45", rtol = 1e-8)
  )

  # Closest index to tf
  closest_index <- which.min(abs(sol[, 1] - tf))
  likelihoods <- unname(sol[closest_index, 1 + 1:nstate])
  
  return(likelihoods)
}


#' Sampling Aware Ancestral State Inference
#'
#' Get the internal node state probabilities of a tree with defined leaf states.
#'
#' @param phy A `phylo` phylogenetic tree. The tree needs to be a rooted binary tree. Must contain
#' `tip.state`. `tip.state` can be names of the states, or numeric values (1:n).
#' @param params_df Data frame containing non-q parameters used in ancestral
#' state reconstruction algorithm. Must have the following column names:
#' `state`, `prior`, `lambda`, `mu`, and `psi`. The `prior` values refer to the baseline
#' probabilities of the states (used at the root of the tree). 
#'
#' Example:
#' | **state** | **prior** | **lambda** | **mu** | **psi** |
#' | :--- | :--- | :--- | :--- | :--- |
#' | 1 | 1/3 | 3 | 0.02 | 1 |
#' | 2 | 1/3 | 3 | 0.02 | 1 |
#' | 3 | 1/3 | 3 | 0.02 | 1 |
#' @param q_matrix Numeric q matrix used in ancestral state reconstruction
#' algorithm. Row and column indices or names represent states.
#' @return A data frame listing the state probabilities of every node in `phy`. The row names correspond
#' to the node IDs. 
#' @export
saasi <- function(phy, params_df, q_matrix) {
  
  # Adding warning messages 
  
  # Checking if the tip states are presented as numeric values
  if(!is.numeric(phy$tip.state)){
    # stop("Convert tip state to numeric values.")
      if( !all(sort(rownames(q_matrix))==sort(params_df$state))) {
          stop("State names and q names need to agree and be present only once") 
      } else {
          phy$tip.statename = phy$tip.state
          phy$tip.state = as.numeric(factor(phy$tip.state)) # makes it numeric 
      params_df$statename = params_df$state
      params_df$state = as.numeric(factor(params_df$statename)) # numeric state
      params_df = params_df[order(params_df$state), ]
      # reorder q if necessary, so its order corresponds to the numeric state 
      q_matrix = q_matrix[ order(as.numeric(factor(rownames(q_matrix)))),
                           order(as.numeric(factor(colnames(q_matrix))))]
      }
  }
  
  # Checking if the tree is binary 
  if(!ape::is.binary.phylo(phy)){
    stop("The phylogenetic tree needs to be a binary tree.")
  }
  
  # Checking if the tree is a rooted tree
  if(!ape::is.rooted.phylo(phy)){
    stop("The phylogenetic tree needs to be a rooted tree.")
  }
  
  
  required_cols <- c("state", "prior", "lambda", "mu", "psi")
  
  # Check if input is a data frame
  if (!is.data.frame(params_df)) {
    stop("Input must be a data frame.")
  }
  
  # Check for required column names
  if (!all(required_cols %in% colnames(params_df))) {
    stop(sprintf("Missing required columns: %s", 
                 paste(setdiff(required_cols, colnames(params_df)), collapse = ", ")))
  }
  
  # Check that 'state' is a sequence of 1-based natural numbers
  if (!is.numeric(params_df$state)) {
    stop("'state' column must be a sequence of numerical values (1, 2, ..., n).")
  }
  
  # Checking the number of unique states in tip.state matches the number of states
  # in params_df$state
  
  if(!length(unique(phy$tip.state)) == length(params_df$state)){
    stop("The number of states does not match.")
  }
  
  nstate <- nrow(params_df)
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
  
  backwards_likelihoods_list <- get_backwards_likelihoods_list(phy,
                                                               params_df,
                                                               q_matrix,
                                                               nstate,
                                                               nnode,
                                                               post_order_edges,
                                                               topology_df)
  #print(backwards_likelihoods_list)
  
  state_probabilities_list <- get_state_probabilities_list(
    phy,
    params_df,
    q_matrix,
    nstate,
    nnode,
    root_node,
    post_order_edges,
    topology_df,
    backwards_likelihoods_list
  )
  
  # if (plot) {
  #   # https://stackoverflow.com/a/17735894
  #   highest_likelihoods <-
  #     apply(state_probabilities_df, 1, function(e) which(e == max(e)))
  #   
  #   # https://colorbrewer2.org/?type=qualitative&scheme=Set3&n=12
  #   # Will serve up to n states == len of vector; afterwards recycles colors
  #   state_colors <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
  #                     "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
  #                     "#ccebc5", "#ffed6f")
  #   # Map likelihoods to colors
  #   highest_likelihood_colors <-
  #     sapply(highest_likelihoods,
  #            function(e) state_colors[[e %% (length(state_colors) + 1)]])
  #   
  #   ape::plot.phylo(phy, label.offset = 0.05)
  #   ape::tiplabels(highest_likelihoods[1:(phy[["Nnode"]] + 1)],
  #                  frame = "circle",
  #                  cex = cex,
  #                  bg = highest_likelihood_colors)
  #   ape::nodelabels(highest_likelihoods[-(1:(phy[["Nnode"]] + 1))],
  #                   frame = "circle",
  #                   cex = cex,
  #                   bg = highest_likelihood_colors)
  # }
  # 
  return(state_probabilities_list)
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
#' list[[x]][[y]] is the likelihood for state y in node x.
#' @noRd
get_backwards_likelihoods_list <- function(phy,
                                           params_df,
                                           q_matrix,
                                           nstate,
                                           nnode,
                                           post_order_edges,
                                           topology_df) {
  backwards_likelihoods_list <- lapply(1:nnode, function(i) rep(0, nstate))
  
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
                                             params_df, q_matrix)
    abs_likelihoods <- abs(likelihoods)
    norm_likelihoods <- abs_likelihoods / sum(abs_likelihoods)
    backwards_likelihoods_list[[node]] <<- norm_likelihoods
    
    #backwards_likelihoods_list[[node]] <<- likelihoods # note changed by CC 
  }))
  
  return(backwards_likelihoods_list)
}

#' Get list of probabilities ultimately generated by ancestral state
#' reconstruction algorithm.
#'
#' @return List of ancestral state probabilities.
#' list[[x]][[y]] is the probability of state y in node x.
#' @noRd
get_state_probabilities_list <- function(phy,
                                         params_df,
                                         q_matrix,
                                         nstate,
                                         nnode,
                                         root_node,
                                         post_order_edges,
                                         topology_df,
                                         backwards_likelihoods_list) {
  #state_probabilities_list <- rep(list(rep(0, nstate)), nnode)
  #state_probabilities_list <- rep(list(rep(0, nstate)), nnode)
  state_probabilities_list <- lapply(1:nnode, function(i) rep(0, nstate))
  # Populate leaf node state probabilities
  invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
    state <- phy[["tip.state"]][[i]]
    state_probabilities_list[[i]][[state]] <<- 1
  }))
  
  state_probabilities_list[[root_node]] <- (
    backwards_likelihoods_list[[root_node]] * params_df$prior
    / (sum(backwards_likelihoods_list[[root_node]] * params_df$prior))
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
                                            params_df, q_matrix)
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

#' Convert state probabilities list to data frame.
#'  List of state probabilities.
#' Data frame of state probabilities in each node in phy.
#' 
# get_state_probabilities_df <- function(phy, nstate, state_probabilities_list) {
#   state_probabilities_df <-
#     as.data.frame(do.call(rbind, state_probabilities_list))
#   row.names(state_probabilities_df) <-
#     c(phy[["tip.label"]], phy[["node.label"]])
#   names(state_probabilities_df) <- seq_len(nstate)
#   return(state_probabilities_df)
# }


#' Helper fn for `probability_density`.
#'
#' @noRd
c1 <- function(lambda, mu, psi) {
  return(abs(sqrt((lambda - mu - psi)^2 + (4 * lambda * psi))))
}

#' Helper fn for `probability_density`.
#'
#' @noRd
c2 <- function(lambda, mu, psi) {
  return(-(lambda - mu - psi) / (c1(lambda, mu, psi)))
}

#' Helper fn for `probability_density`.
#'
#' @noRd
q <- function(lambda, mu, psi, t) {
  d1 <- 2 * (1 - (c2(lambda, mu, psi))^2)
  d2 <- exp(-c1(lambda, mu, psi) * t) * (1 - c2(lambda, mu, psi))^2
  d3 <- exp(c1(lambda, mu, psi) * t) * (1 + c2(lambda, mu, psi))^2
  return(d1 + d2 + d3)
}

#' Helper fn for `mle_lm`.
#'
#' @noRd
probability_density <- function(lambda, mu, psi,
                                internal_node_times, leaf_node_times) {
  # Add checks for invalid inputs
  sorted_x <- sort(internal_node_times)
  sorted_y <- sort(leaf_node_times)
  d1 <- length(internal_node_times) * log(lambda)
  d2 <- sum(-log(q(lambda, mu, psi, sorted_x)))
  d3 <- sum(log(psi) + log(q(lambda, mu, psi, sorted_y)))
  ret <- d1 + d2 + d3
  return(ret)
}

#' Estimate speciation/extinction rates for a tree
#'
#' This is done with a maximum likelihood method,  implemented mainly by
#' \href{https://github.com/yexuansong}{@yexuansong}, from methods described in
#' \href{https://doi.org/10.1093/molbev/msr217}{Stadler et al. (2012)}.
#'
#' @param phy A `phylo` phylogenetic tree (`ape` format).
#' @param lambda An initial "guess" for speciation, used in subsequent formulae.
#' @param mu An initial "guess" for extinction, used in subsequent formulae.
#' @param psi Sampling rate.
#' @param method See `method` parameter in \link{mle}.
#' @param lower A two-element vector containing the lower bound of the speciation 
#' and extinction values (geater than 0), the default is c(0.001,0.001).
#' @param upper A two-element vector containing the upper bound of the speciation 
#' and extinction values, should be similar to the sampling rate, still need to find
#' a better justification for the upcoming version.
#' @return A two-element vector containing speciation and extinction values
#' estimated from `phy`.
#'@export
mle_lm <- function(phy, lambda, mu, psi, method = "L-BFGS-B", lower = c(0.001,0.001), upper) {
  node_depths <- ape::node.depth.edgelength(phy)
  node_times <- max(node_depths) - node_depths
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- phy[["Nnode"]] * 2 + 1
  nleaf_node <- nnode - phy[["Nnode"]]
  leaf_node_times <- node_times[1:nleaf_node]
  internal_node_times <- node_times[(nleaf_node + 1):nnode]
  
  negative_log_likelihood <- function(lambda, mu) {
    positive_ret <- probability_density(lambda, mu, psi,
                                        internal_node_times, leaf_node_times)
    return(-positive_ret)
  }
  
  fit <- stats4::mle(negative_log_likelihood,
                     start = list(lambda = lambda, mu = mu),
                     method = method,
                     lower = lower,
                     upper = upper)
  
  ret <- unname(c(stats4::coef(fit)[1], stats4::coef(fit)[2]))
  names(ret) = c("speciation","extinction")
  return(ret)
}


args <- commandArgs(trailingOnly = TRUE)
sim_num <- as.numeric(args[1])
ntaxa <- as.numeric(args[2])
max_time <- as.numeric(args[3])
output_dir <- args[4]

set.seed(sim_num)

# Generate tree
pars <- data.frame(state=c(1,2), prior=c(0.5,0.5), lambda=c(3,1.5), mu=c(0.1,0.05), psi=c(0.5,0.5))
q_matrix <- matrix(0.3, 2, 2)
diag(q_matrix) <- NA


phy <- sim_bds_tree(pars, q_matrix, x0=1, max_taxa=ntaxa, max_t=max_time, include_extinct=FALSE)
h <- history.from.sim.discrete(phy, 1:2)

# Process tree
node_depths <- node.depth.edgelength(phy)
tmrca <- max(node_depths)
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
if(length(tips_to_drop) > 0) {
    phy <- drop.tip(phy, tips_to_drop)
    phy$tip.state <- phy$tip.state[setdiff(names(phy$tip.state), tips_to_drop)]
}

phy_info <- as_tibble(phy)%>% mutate(State = c(factor(h$tip.state),factor(h$node.state))[label])
write.csv(phy_info, file=paste0(output_dir,"/true_tree_",sim_num,".csv"))

# Export tree and metadata
write.tree(phy, file=paste0(output_dir, "/tree_", sim_num, ".nwk"))
metadata_df <- data.frame(
    node = phy$tip.label,
    state = paste0("State", phy$tip.state),
    numeric_state = phy$tip.state
)
write.csv(metadata_df, file=paste0(output_dir, "/metadata_", sim_num, ".csv"), row.names=FALSE)
write.table(metadata_df, file=paste0(output_dir, "/metadata_", sim_num, ".txt"), sep=",", row.names=FALSE)

timing_results <- data.frame(
    simulation = integer(0),
    method = character(0),
    time_seconds = numeric(0),
    success = logical(0)
)

# ACE
ace_start <- Sys.time()
ace_success <- FALSE
ace_result <- ace(phy$tip.state,phy,type="discrete",model="ER")
write.csv(ace_result$lik.anc, file=paste0(output_dir,"/ace_",sim_num,".csv"))
ace_success <- TRUE
ace_time <- as.numeric(Sys.time() - ace_start)
timing_results <- rbind(timing_results, data.frame(simulation=sim_num, method="ace", time_seconds=ace_time, success=ace_success))

# SAASI
saasi_start <- Sys.time()
saasi_success <- FALSE
tryCatch({
    saasi_result <- saasi(phy, pars, q_matrix)
    write.csv(saasi_result, file=paste0(output_dir, "/saasi_", sim_num, ".csv"))
    saasi_success <- TRUE
}, error = function(e) {})
saasi_time <- as.numeric(Sys.time() - saasi_start)
timing_results <- rbind(timing_results, data.frame(simulation=sim_num, method="saasi", time_seconds=saasi_time, success=saasi_success))

# SAASI with pars estimations
saasi_par_start <- Sys.time()
saasi_par_success <- FALSE
tryCatch({
	estimates <- mle_lm(phy,lambda = 2, mu = 0.1, psi = 1, lower = c(0.001,0.001), upper = c(5,5))
	q_matrix_est <- extract_ace_q(ace_result)
	pars_est <- data.frame(state=c(1,2),
                       prior=c(0.5,0.5),
                       lambda=rep(estimates[1],2),
                       mu=rep(estimates[2],2),
                       psi=c(.5,.5))
	saasi_est_result <- saasi(phy,pars_est, q_matrix_est)
	write.csv(saasi_est_result, file=paste0(output_dir,"/saasi_pars_est_",sim_num,".csv"))
	saasi_par_success <- TRUE
}, error = function(e){})
saasi_par_est_time <- as.numeric(Sys.time() - saasi_par_start)
timing_results <- rbind(timing_results,data.frame(simulation=sim_num,method = "saasi_par_est",time_seconds=saasi_par_est_time,success=saasi_par_success))

# SIMMAP
simmap_start <- Sys.time()
simmap_success <- FALSE
tryCatch({
    # Calculate custom pi corrected for sampling bias
    rho_1 <- 0.25  # sampling rate for state 1
    rho_2 <- 1.0   # sampling rate for state 2
    
    # Observed counts in tree
    obs_counts <- table(phy$tip.state)
    n_1 <- as.numeric(obs_counts["1"])
    n_2 <- as.numeric(obs_counts["2"])
    
    # Correct for sampling bias
    true_n_1 <- n_1 / rho_1
    true_n_2 <- n_2 / rho_2
    
    # Calculate corrected pi
    pi_corrected <- c(true_n_1, true_n_2) / sum(true_n_1, true_n_2)
    
    # Use corrected pi in fitMk
    smap_model <- fitMk(phy, phy$tip.state, model="ER", pi=pi_corrected)
    phy_smap <- simmap(smap_model)
    sm <- summary(phy_smap)
    simmap_lik <- as.data.frame(sm$ace)
    write.csv(simmap_lik, file=paste0(output_dir, "/simmap_", sim_num, ".csv"))
    simmap_success <- TRUE
}, error = function(e) {})
simmap_time <- as.numeric(Sys.time() - simmap_start)
timing_results <- rbind(timing_results, data.frame(simulation=sim_num, method="simmap", time_seconds=simmap_time, success=simmap_success))

# Save timing
write.csv(timing_results, file=paste0(output_dir, "/timing_", sim_num, ".csv"), row.names=FALSE)
EOF

# Run simulations
for i in $(seq 1 $N_SIMS); do
    # Create subfolder for this simulation
    SIM_DIR="$OUTPUT_DIR/sim_$i"
    mkdir -p "$SIM_DIR"
    
    # Run R methods
    Rscript "$OUTPUT_DIR/run_asr.R" $i $NTAXA $MAX_TIME "$SIM_DIR"
    
    # TreeTime
    cd "$SIM_DIR"
    mkdir -p treetime_results
    treetime_start=$(date +%s)
    treetime mugration --tree "tree_${i}.nwk" --states "metadata_${i}.csv" --attribute numeric_state --confidence --outdir "treetime_results" --verbose 0 > /dev/null 2>&1
    treetime_success=$?
    treetime_time=$(($(date +%s) - treetime_start))
    
    # PastML
    mkdir -p pastml_results
    pastml_start=$(date +%s)
    pastml --tree "tree_${i}.nwk" --data "metadata_${i}.txt" --columns numeric_state --work_dir "pastml_results" --data_sep "," > /dev/null 2>&1
    pastml_success=$?
    pastml_time=$(($(date +%s) - pastml_start))
    
    cd - > /dev/null
    
    # Save external method timing in simulation folder
    echo "simulation,method,time_seconds,success" > "$SIM_DIR/external_timing.csv"
    echo "$i,treetime,$treetime_time,$([[ $treetime_success -eq 0 ]] && echo true || echo false)" >> "$SIM_DIR/external_timing.csv"
    echo "$i,pastml,$pastml_time,$([[ $pastml_success -eq 0 ]] && echo true || echo false)" >> "$SIM_DIR/external_timing.csv"
done
