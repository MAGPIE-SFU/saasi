# ============================================================================
# COMPLETE FUNCTIONS NEEDED FOR TIME-DEPENDENT TREE SIMULATION
# ============================================================================

#' Create time-dependent psi parameters using simple data frame approach
#'
#' @param time_points Vector of time points where psi changes
#' @param psi_matrix Matrix where each row corresponds to psi values at each time point
#'                   Columns represent states
#' @return Data frame ready to use with saasi
create_psi_stepwise <- function(time_points, psi_matrix) {
  # Validate inputs
  if (length(time_points) != nrow(psi_matrix)) {
    stop("Number of time points must match number of rows in psi_matrix")
  }
  
  if (!is.numeric(time_points)) {
    stop("time_points must be numeric")
  }
  
  if (!is.numeric(psi_matrix)) {
    stop("psi_matrix must be numeric")
  }
  
  # Create the data frame
  psi_df <- data.frame(
    time = time_points,
    psi_matrix
  )
  
  # Sort by time
  psi_df <- psi_df[order(psi_df$time), ]
  
  # Reset row names
  rownames(psi_df) <- NULL
  
  return(psi_df)
}

#' Get psi values at specific time using step function (constant within intervals)
#'
#' @param time_point Numeric value representing the time point
#' @param psi_df Data frame with time-dependent psi values. 
#'               First column: time, remaining columns: psi values for each state
#' @return Vector of psi values for each state at the given time point
get_psi_at_time_stepwise <- function(time_point, psi_df) {
  # Input validation
  if (!is.numeric(time_point)) {
    stop("time_point must be numeric, got: ", class(time_point))
  }
  
  if (!is.data.frame(psi_df)) {
    stop("psi_df must be a data frame, got: ", class(psi_df))
  }
  
  if (nrow(psi_df) == 0) {
    stop("psi_df is empty")
  }
  
  # Get time column (first column) and sort by time
  time_col <- psi_df[, 1]
  psi_cols <- psi_df[, -1, drop = FALSE]
  
  # Sort by time if not already sorted
  if (is.unsorted(time_col)) {
    sort_order <- order(time_col)
    time_col <- time_col[sort_order]
    psi_cols <- psi_cols[sort_order, , drop = FALSE]
  }
  
  # Find which interval the time_point falls into
  # Time intervals: [0, time1), [time1, time2), ..., [timeN, max_tree_age)
  
  if (time_point < time_col[1]) {
    # Before first time point - use first interval values
    result <- as.numeric(psi_cols[1, ])
  } else if (time_point >= time_col[length(time_col)]) {
    # After last time point - use last interval values
    result <- as.numeric(psi_cols[nrow(psi_cols), ])
  } else {
    # Find the interval: time_point is in [time_col[i], time_col[i+1])
    interval_idx <- max(which(time_col <= time_point))
    result <- as.numeric(psi_cols[interval_idx, ])
  }
  
  return(result)
}

#' Get psi values at specific simulation time using step function
#'
#' @param sim_time Current simulation time (forward time from 0)
#' @param psi_df Data frame with time-dependent psi values
#' @param max_t Maximum simulation time
#' @return Vector of psi values for each state at the given simulation time
get_psi_at_sim_time <- function(sim_time, psi_df, max_t) {
  # Convert simulation time to tree time (simulation runs forward, tree time runs backward)
  tree_time <- max_t - sim_time
  
  # Use the existing stepwise function
  return(get_psi_at_time_stepwise(tree_time, psi_df))
}

#' MODIFIED: Tree simulation with time-dependent psi
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
  
  # MODIFIED: Initial rate calculation
  if (is.null(psi_df)) {
    # Use constant psi from pars matrix
    r.i <- rowSums(pars)
  } else {
    # Use time-dependent psi - start with psi at time 0
    initial_psi <- get_psi_at_sim_time(0, psi_df, max.t)
    print(initial_psi)
    
    # Replace psi column in pars with current psi values
    pars_current <- pars
    pars_current[, 3] <- initial_psi  # Column 3 is psi
    r.i <- rowSums(pars_current)
  }
  
  len <- 0
  t <- 0
  hist <- list()
  
  # Always single lineage in our code
  states <- x0
  # We use 1-based states
  n.taxa <- lineages <- n.i[x0] <- 1
  start <- 0
  
  while (n.taxa <= max.taxa && n.taxa > 0) {
    # MODIFIED: Update psi values based on current simulation time
    if (!is.null(psi_df)) {
      current_psi <- get_psi_at_sim_time(t, psi_df, max.t)
      print(current_psi)
      pars_current <- pars
      pars_current[, 3] <- current_psi  # Update psi column
      r.i <- rowSums(pars_current)
    } else {
      pars_current <- pars
    }
    
    ## When does an event happen?
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    
    if (r.tot <= 0) {
      # No more events possible
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
    
    # type: 1: speciation, 2: extinction, 3: sampling, >3: char change
    state <- sample(k, 1, FALSE, r.n / r.tot)
    
    ## Pick a lineage for that state:
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
      ## Sampling; this is new
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
  
  # ADDED: Store the psi function used for simulation
  if (!is.null(psi_df)) {
    attr(info, "psi_df") <- psi_df
  }
  
  info
}

#' You'll need your existing me.to.ape.bisse function here
#' (Copy it from your original tree_generator.R file)

#' MODIFIED: Simulate a birth/death/sampling tree with time-dependent psi
#'
#' @param params_df Data frame with parameters. Can have time-dependent psi.
#' @param q_matrix Transition rate matrix
#' @param x0 Root state
#' @param max_taxa Maximum number of taxa
#' @param max_t Maximum simulation time
#' @param include_extinct Include extinct lineages
#' @param psi_df Optional: Data frame with time-dependent psi values.
#' @return A phylo object with time-dependent sampling history
sim_bds_tree_time_psi <- function(params_df, q_matrix, x0, max_taxa = 100, max_t = 100,
                                  include_extinct = FALSE, psi_df = NULL) {
  k <- nrow(params_df)
  
  # Extract psi values
  if (is.null(psi_df)) {
    # Check if params_df has time-dependent psi
    if (is.list(params_df$psi) && is.data.frame(params_df$psi[[1]])) {
      psi_df <- params_df$psi[[1]]
      cat("Using time-dependent psi from params_df\n")
    } else {
      cat("Using constant psi from params_df\n")
    }
  } else {
    cat("Using provided psi_df for time-dependent psi\n")
  }
  
  # Create parameter matrix - use average psi for initial setup
  if (!is.null(psi_df)) {
    # Use average psi values for matrix construction
    avg_psi <- colMeans(psi_df[, -1])  # Average across time, exclude time column
    temp_params_df <- params_df
    temp_params_df$psi <- avg_psi
    pars <- cbind(
      data.matrix(temp_params_df[3:5]),
      matrix(t(q_matrix)[col(q_matrix) != row(q_matrix)], k, byrow = TRUE)
    )
  } else {
    pars <- cbind(
      data.matrix(params_df[3:5]),
      matrix(t(q_matrix)[col(q_matrix) != row(q_matrix)], k, byrow = TRUE)
    )
  }
  
  phy <- NULL
  attempt <- 1
  max_attempts <- 50
  
  while (is.null(phy) && attempt <= max_attempts) {
    cat("Simulation attempt", attempt, "\n")
    
    tryCatch({
      info <- make.tree.bisse_modified_time_psi(pars, k, max_taxa, max_t, x0, psi_df)
      phy <- me.to.ape.bisse(info[-1, ], info$state[1])
      
      if (!is.null(phy)) {
        # Add psi information to the tree
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

#' Create parameters for simulation with time-dependent psi
#'
#' @param states Vector of state numbers
#' @param lambda Speciation rates for each state
#' @param mu Extinction rates for each state  
#' @param psi_time_points Vector of time points where psi changes
#' @param psi_values Matrix where each row = psi values for all states at that time
#' @param prior Prior probabilities for each state
#' @return List with params_df and psi_df ready for simulation
create_sim_params_time_psi <- function(states, lambda, mu, psi_time_points, psi_values, prior = NULL) {
  n_states <- length(states)
  
  if (is.null(prior)) {
    prior <- rep(1/n_states, n_states)
  }
  
  # Create psi data frame
  psi_df <- create_psi_stepwise(psi_time_points, psi_values)
  
  # Create params_df with average psi (for compatibility)
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
