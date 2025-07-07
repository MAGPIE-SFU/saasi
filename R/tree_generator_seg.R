# ============================================================================
# FIXED FUNCTIONS WITH ORIGINAL NAMES
# ============================================================================

#' Get psi values at specific simulation time using forward time (FIXED)
#'
#' @param sim_time Current simulation time (starts at 0, increases)
#' @param psi_df Data frame with time-dependent psi values in SIMULATION time
#' @param max_t Maximum simulation time (not used anymore, kept for compatibility)
#' @return Vector of psi values for each state at the given simulation time
get_psi_at_sim_time <- function(sim_time, psi_df, max_t) {
  # FIXED: NO time conversion - use simulation time directly
  
  time_col <- psi_df[, 1]
  psi_cols <- psi_df[, -1, drop = FALSE]
  
  # Sort by time if not already sorted
  if (is.unsorted(time_col)) {
    sort_order <- order(time_col)
    time_col <- time_col[sort_order]
    psi_cols <- psi_cols[sort_order, , drop = FALSE]
  }
  
  # Find which interval the sim_time falls into
  # Time intervals: [0, time1), [time1, time2), ..., [timeN, ∞)
  
  if (sim_time < time_col[1]) {
    # Before first time point - use first interval values
    result <- as.numeric(psi_cols[1, ])
  } else if (sim_time >= time_col[length(time_col)]) {
    # After last time point - use last interval values
    result <- as.numeric(psi_cols[nrow(psi_cols), ])
  } else {
    # Find the interval: sim_time is in [time_col[i], time_col[i+1])
    interval_idx <- max(which(time_col <= sim_time))
    result <- as.numeric(psi_cols[interval_idx, ])
  }
  
  return(result)
}

#' FIXED: Tree simulation with time-dependent psi (original function name)
#'
#' @param pars Parameter matrix (lambda, mu, psi, q_ij...)
#' @param k Number of states
#' @param max.taxa Maximum number of taxa
#' @param max.t Maximum simulation time
#' @param x0 Root state
#' @param psi_df Optional data frame with time-dependent psi values
make.tree.bisse_modified_time_psi <- function(pars, k, max.taxa = Inf, max.t = Inf, x0, psi_df = NULL) {
  # Matrix representation of state changes to consider
  to <- matrix(unlist(lapply(1 : k, function(i) (1 : k)[-i])), k, k - 1, TRUE)
  
  extinct <- FALSE
  split   <- FALSE
  sam <- FALSE
  
  parent <- 0
  n.i <- rep(0, k)
  
  # Initial rate calculation
  if (is.null(psi_df)) {
    r.i <- rowSums(pars)
  } else {
    # FIXED: Use simulation time directly (no conversion)
    initial_psi <- get_psi_at_sim_time(0, psi_df, max.t)
    pars_current <- pars
    pars_current[, 3] <- initial_psi  # Column 3 is psi
    r.i <- rowSums(pars_current)
  }
  
  len <- 0
  t <- 0  # Simulation time starts at 0 and increases
  hist <- list()
  
  states <- x0
  n.taxa <- lineages <- n.i[x0] <- 1
  start <- 0
  
  while (n.taxa <= max.taxa && n.taxa > 0) {
    # FIXED: Update psi values based on current simulation time (forward)
    if (!is.null(psi_df)) {
      current_psi <- get_psi_at_sim_time(t, psi_df, max.t)
      print(t)
      print(current_psi)
      pars_current <- pars
      pars_current[, 3] <- current_psi  # Update psi column
      r.i <- rowSums(pars_current)
      
      # Debug output
      if (t %% 10 < 0.1 || t == 0) {  # Print at start and every ~10 time units
        cat("Sim time:", round(t, 2), "Psi:", round(current_psi, 3), "\n")
      }
    } else {
      pars_current <- pars
    }
    
    ## When does an event happen?
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    
    if (r.tot <= 0) {
      break
    }
    
    dt <- rexp(1, r.tot)
    t <- t + dt
    
    if (t > max.t) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }
    
    len[lineages] <- len[lineages] + dt
    
    # Rest of the simulation logic unchanged
    state <- sample(k, 1, FALSE, r.n / r.tot)
    j <- sample(n.i[state], 1)
    lineage <- lineages[states[lineages] == state][j]
    type <- sample(3 + (k - 1), 1, FALSE, pars_current[state, ])
    
    if (type == 1) {
      ## Speciating:
      if (n.taxa == max.taxa)
        break
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
      ## Extinct
      extinct[lineage] <- TRUE
      lineages <- which(!split & !extinct & !sam)
      n.i[state] <- n.i[state] - 1
      n.taxa <- n.taxa - 1
    } else if (type == 3) {
      ## Sampling
      sam[lineage] <- TRUE
      n.taxa <- n.taxa - 1
      n.i[state] <- n.i[state] - 1
      lineages <- which(!split & !extinct & !sam)
    } else {
      ## Character switch:
      states[lineage] <- state.new <- to[state, type - 3]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1, -1)
      hist[[length(hist) + 1]] <- c(lineage, t, state, state.new)
    }
  }
  
  info <- data.frame(idx = seq_along(extinct), len = len, parent = parent,
                     start = start, state = states, extinct = extinct,
                     split = split, sam = sam)
  
  hist <- as.data.frame(do.call(rbind, hist))
  if (nrow(hist) == 0)
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0
  
  attr(info, "t") <- t
  attr(info, "hist") <- hist
  
  if (!is.null(psi_df)) {
    attr(info, "psi_df") <- psi_df
  }
  
  info
}


#' Copy/pasted from private `diversitree::me.to.ape.bisse` function.
#'
#' I did introduce minor changes to the original code to fix linting issues.
#' However, I did not change variable names to snake_case.
#'
#' @noRd
me.to.ape.bisse <- function(x, root.state) {
  if (nrow(x) == 0)
    return(NULL)
  Nnode <- sum(!x$split) - 1
  n.tips <- sum(!x$split)
  
  x$idx2 <- NA
  x$idx2[!x$split] <- 1:n.tips
  x$idx2[x$split] <- order(x$idx[x$split]) + n.tips + 1
  
  i <- match(x$parent, x$idx)
  x$parent2 <- x$idx2[i]
  x$parent2[is.na(x$parent2)] <- n.tips + 1
  
  tip.label <- ifelse(subset(x, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)
  
  x$name <- NA
  x$name[!x$split] <- tip.label
  ## More useful, but I don't want to clobber anything...
  x$name2 <- c(tip.label, node.label)[x$idx2]
  
  tip.state <- x$state[match(1:n.tips, x$idx2)]
  names(tip.state) <- tip.label
  
  node.state <- x$state[match(1:Nnode + n.tips, x$idx2)]
  names(node.state) <- node.label
  node.state["nd1"] <- root.state
  
  hist <- attr(x, "hist")
  if (!is.null(hist)) {
    hist$idx2 <- x$idx2[match(hist$idx, x$idx)]
    hist$name2 <- x$name2[match(hist$idx, x$idx)]
    if (nrow(hist) > 0)
      hist <- hist[order(hist$idx2), ]
  }
  
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


#' FIXED: Simulate tree with time-dependent psi (original function name)
#'
#' @param params_df Data frame with parameters
#' @param q_matrix Transition rate matrix  
#' @param x0 Root state
#' @param max_taxa Maximum number of taxa
#' @param max_t Maximum simulation time
#' @param include_extinct Include extinct lineages
#' @param psi_df Data frame with time-dependent psi in SIMULATION time
sim_bds_tree_time_psi <- function(params_df, q_matrix, x0, max_taxa = 100, max_t = 100,
                                  include_extinct = FALSE, psi_df = NULL) {
  k <- nrow(params_df)
  
  # Create parameter matrix
  if (!is.null(psi_df)) {
    # Use average psi values for initial matrix construction
    avg_psi <- colMeans(psi_df[, -1])
    temp_params_df <- params_df
    temp_params_df$psi <- avg_psi
    pars <- cbind(
      data.matrix(temp_params_df[3:5]),
      matrix(t(q_matrix)[col(q_matrix) != row(q_matrix)], k, byrow = TRUE)
    )
    cat("Using time-dependent psi (forward time)\n")
  } else {
    pars <- cbind(
      data.matrix(params_df[3:5]),
      matrix(t(q_matrix)[col(q_matrix) != row(q_matrix)], k, byrow = TRUE)
    )
    cat("Using constant psi\n")
  }
  
  phy <- NULL
  attempt <- 1
  max_attempts <- 10
  
  while (is.null(phy) && attempt <= max_attempts) {
    cat("Simulation attempt", attempt, "\n")
    
    tryCatch({
      info <- make.tree.bisse_modified_time_psi(pars, k, max_taxa, max_t, x0, psi_df)
      phy <- me.to.ape.bisse(info[-1, ], info$state[1])
      
      if (!is.null(phy)) {
        if (!is.null(psi_df)) {
          phy$psi_df <- psi_df
          phy$max_t <- max_t
        }
        
        if (!include_extinct) {
          phy <- diversitree::prune(phy)
        }
      }
    }, error = function(e) {
      cat("Simulation attempt", attempt, "failed:", e$message, "\n")
      phy <<- NULL
    })
    
    attempt <- attempt + 1
  }
  
  if (is.null(phy)) {
    stop("Failed to generate tree after ", max_attempts, " attempts")
  }
  
  return(phy)
}

#' FIXED: Create parameters for simulation (original function name)
#'
#' @param states Vector of state numbers
#' @param lambda Speciation rates for each state
#' @param mu Extinction rates for each state
#' @param psi_time_points Vector of SIMULATION time points where psi changes
#' @param psi_values Matrix where each row = psi values for all states at that time
#' @param prior Prior probabilities for each state
create_sim_params_time_psi <- function(states, lambda, mu, psi_time_points, psi_values, prior = NULL) {
  n_states <- length(states)
  
  if (is.null(prior)) {
    prior <- rep(1/n_states, n_states)
  }
  
  # FIXED: Create psi data frame with simulation time (no conversion)
  psi_df <- data.frame(
    time = psi_time_points,
    psi_values
  )
  names(psi_df) <- c("time", paste0("X", states))
  
  # Sort by time
  psi_df <- psi_df[order(psi_df$time), ]
  rownames(psi_df) <- NULL
  
  # Create params_df with average psi
  avg_psi <- colMeans(psi_values)
  
  params_df <- data.frame(
    state = states,
    prior = prior,
    lambda = lambda,
    mu = mu,
    psi = avg_psi
  )
  
  return(list(
    params_df = params_df,
    psi_df = psi_df
  ))
}

