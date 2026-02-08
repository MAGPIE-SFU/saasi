#' Adapted from private `diversitree::make.tree.bisse` function with modifications
#' to incorporate sampling events and support three or more discrete states.
#'
#' @noRd
make.tree.bisse_modified <- function(pars, k, max.taxa = Inf, max.t = Inf, x0) {
  # Create matrix of possible state transitions (excluding self-transitions)
  to <- matrix(unlist(lapply(1 : k, function(i) (1 : k)[-i])), k, k - 1, TRUE)
  
  # Initialize tracking variables for lineage states
  extinct <- FALSE
  split   <- FALSE
  sam <- FALSE  # Track sampling events
  
  parent <- 0
  n.i <- rep(0, k)  # Number of lineages in each state
  r.i <- rowSums(pars)  # Total rate for each state
  len <- 0  # Edge lengths
  t <- 0  # Current time
  hist <- list()  # History of state changes
  
  # Initialize with single lineage at root state (always 1-based indexing)
  states <- x0
  n.taxa <- lineages <- n.i[x0] <- 1
  start <- 0
  
  while (n.taxa <= max.taxa && n.taxa > 0) {
    ## Calculate waiting time until next event
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    dt <- stats::rexp(1, r.tot)
    t <- t + dt
    
    # Check if simulation exceeds maximum time
    if (t > max.t) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }
    
    len[lineages] <- len[lineages] + dt
    
    # Determine event type and affected lineage
    # Event types: 1=speciation, 2=extinction, 3=sampling, >3=character change
    state <- sample(k, 1, FALSE, r.n / r.tot)
    
    ## Select random lineage in the chosen state
    j <- sample(n.i[state], 1)
    lineage <- lineages[states[lineages] == state][j]
    type <- sample(3 + (k - 1), 1, FALSE, pars[state, ])
    
    if (type == 1) {
      ## Speciation event
      if (n.taxa == max.taxa)
        break  # Reached maximum taxa limit
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      sam[new.i] <- FALSE
      states[new.i] <- state
      parent[new.i] <- lineage
      start[new.i] <- t
      len[new.i] <- 0
      
      n.i[state] <- n.i[state] + 1
      n.taxa <- n.taxa + 1
      
      lineages <- which(!split & !extinct & !sam)
    } else if (type == 2) {
      ## Extinction event: lineage goes extinct
      extinct[lineage] <- TRUE
      
      lineages <- which(!split & !extinct & !sam)
      
      n.i[state] <- n.i[state] - 1
      n.taxa <- n.taxa - 1
    } else if (type == 3) {
      ## Sampling event: lineage is sampled (becomes a tip in non-ultrametric tree)
      sam[lineage] <- TRUE
      n.taxa <- n.taxa - 1
      n.i[state] <- n.i[state] - 1
      lineages <- which(!split & !extinct & !sam)
    } else {
      ## Character state change: transition to a different state
      states[lineage] <- state.new <- to[state, type - 3]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1, -1)
      hist[[length(hist) + 1]] <- c(lineage, t, state, state.new)
    }
  }
  
  # Compile simulation results into structured data frame
  info <- data.frame(idx = seq_along(extinct), len = len, parent = parent,
                     start = start, state = states, extinct = extinct,
                     split = split, sam = sam)
  
  # Format state change history
  hist <- as.data.frame(do.call(rbind, hist))
  if (nrow(hist) == 0)
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0
  
  # Attach simulation metadata
  attr(info, "t") <- t
  attr(info, "hist") <- hist
  info
}

#' Convert simulation results to ape phylo format
#'
#' Adapted from private `diversitree::me.to.ape.bisse` function.
#' Converts internal simulation data structure to standard `ape` phylo object.
#'
#' Note: Some code style preferences are not applied to maintain consistency
#' with the original `diversitree` implementation.
#'
#' @noRd
me.to.ape.bisse <- function(x, root.state) {
  if (nrow(x) == 0)
    return(NULL)
  Nnode <- sum(!x$split) - 1
  n.tips <- sum(!x$split)
  
  # Create index mapping for nodes and tips
  x$idx2 <- NA
  x$idx2[!x$split] <- 1:n.tips
  x$idx2[x$split] <- order(x$idx[x$split]) + n.tips + 1
  
  # Map parent indices to phylo format
  i <- match(x$parent, x$idx)
  x$parent2 <- x$idx2[i]
  x$parent2[is.na(x$parent2)] <- n.tips + 1
  
  # Generate labels for tips (extinct vs sampled)
  tip.label <- ifelse(subset(x, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)
  
  x$name <- NA
  x$name[!x$split] <- tip.label
  x$name2 <- c(tip.label, node.label)[x$idx2]
  
  # Extract tip and node states from simulation
  tip.state <- x$state[match(1:n.tips, x$idx2)]
  names(tip.state) <- tip.label
  
  node.state <- x$state[match(1:Nnode + n.tips, x$idx2)]
  names(node.state) <- node.label
  node.state["nd1"] <- root.state
  
  # Process state change history
  hist <- attr(x, "hist")
  if (!is.null(hist)) {
    hist$idx2 <- x$idx2[match(hist$idx, x$idx)]
    hist$name2 <- x$name2[match(hist$idx, x$idx)]
    if (nrow(hist) > 0)
      hist <- hist[order(hist$idx2), ]
  }
  
  # Construct phylo object with all simulation metadata
  phy <- ape::reorder.phylo(structure(
    list(
      edge = cbind(x$parent2, x$idx2),
      Nnode = Nnode,
      tip.label = tip.label,
      tip.state = tip.state,
      node.label = node.label,
      node.state = node.state,
      edge.length = x$len,
      orig = x,
      hist = hist
    ),
    class = "phylo"
  ))
  phy$edge.state <- x$state[match(phy$edge[, 2], x$idx2)]
  phy
}

#' Simulate a birth/death/sampling tree with post-processing
#'
#' This tree can be passed to [saasi::saasi], and will include speciation,
#' extinction, sampling, and mutation events. The tree is post-processed to
#' remove tips at the present and ensure a minimum number of tips.
#'
#' @param x0 Natural number used as root state in the returned tree. Must be a
#' state declared in `params_df`.
#' @param max_taxa Maximum number of nodes allowed in the initial simulated tree.
#' @param max_t Maximum depth allowed in the returned tree.
#' @param include_extinct Boolean declaring whether extinct taxa are included in
#' the returned tree.
#' @param min_tip Minimum number of tips required in the post-processed tree.
#' If the processed tree has fewer tips, the simulation is repeated until this
#' condition is satisfied. Default is 1.
#' @return A `phylo` phylogenetic tree (`ape` format) with post-processing applied.
#' The tree includes `tip.state` attribute with character state labels.
#' @inheritParams saasi
#' @export
sim_bds_tree <- function(params_df, q_matrix, x0, max_taxa = 100, max_t = 100,
                         include_extinct = FALSE, min_tip = 1) {
  # Format parameters matrix for birth-death-sampling simulation
  k <- nrow(params_df)
  pars <- cbind(
    data.matrix(params_df[3:5]),
    matrix(t(q_matrix)[col(q_matrix) != row(q_matrix)], k, byrow = TRUE)
  )
  
  # Repeat simulation until we get a valid tree with minimum number of tips
  phy <- NULL
  n_tips <- 0
  
  while (is.null(phy) || n_tips < min_tip) {
    # Simulate initial tree
    info <- make.tree.bisse_modified(pars, k, max_taxa, max_t, x0)
    phy <- me.to.ape.bisse(info[-1, ], info$state[1])
    
    # Prune extinct lineages if requested
    if (!is.null(phy) && !include_extinct) {
      phy <- diversitree::prune(phy)
    }
    
    # Post-processing: extract history and identify tips to drop
    if (!is.null(phy)) {

      h <- diversitree::history.from.sim.discrete(phy, 1:k)
      

      original_node_labels <- phy$node.label
      original_node_state <- h$node.state  # All original node states
      
      node_depths <- ape::node.depth.edgelength(phy)
      tmrca <- max(node_depths)
      
      # Identify tips at the root (within 0.01 time units of TMRCA)
      tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
      
      # Drop tips if any exist
      if (length(tips_to_drop) > 0) {
        phy <- ape::drop.tip(phy, tips_to_drop)
        
        # tip.state: keep only tips that weren't dropped
        phy$tip.state <- phy$tip.state[setdiff(names(phy$tip.state), tips_to_drop)]
        
        phy$node.state <- original_node_state[match(phy$node.label, original_node_labels)]
        
      }

      phy$tip.state <- as.character(phy$tip.state)
      
      n_tips <- length(phy$tip.label)
      
      if (n_tips < min_tip) {
        phy <- NULL
      }
    }
  }
  
  return(phy)
}