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

#' @param child_likelihoods Child likelihood values
#' @param t0 Child age
#' @param tf Parent age 
#' @param params_df Parameters data frame
#' @param q_matrix Transition rate matrix
#' @return Vector of likelihoods at parent
#' @noRd
backwards_likelihoods_helper <- function(child_likelihoods,
                                         t0, 
                                         tf, 
                                         params_df, 
                                         q_matrix) {
  nstate <- nrow(params_df)
  
  # Extract psi data and determine if it's time-dependent
  psi_info <- extract_psi_data(params_df)
  
  if (!psi_info$is_time_dependent) {
    # Constant sampling rate - use original approach
    return(backwards_likelihoods_helper_constant(child_likelihoods, t0, tf, params_df, q_matrix, psi_info$psi_data))
  } else {
    # Time-dependent sampling rate - use segmented approach
    return(backwards_likelihoods_helper_segmented(child_likelihoods, t0, tf, params_df, q_matrix, psi_info$psi_data))
  }
}

#' Backwards helper for constant sampling rate  
#'
#' @noRd
backwards_likelihoods_helper_constant <- function(child_likelihoods,
                                                  t0, 
                                                  tf, 
                                                  params_df, 
                                                  q_matrix, 
                                                  psi_values) {
  # Number of states == n
  nstate <- nrow(params_df)
  
  func <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      # Use constant psi throughout
      current_psi <- psi_values
      
      # States 1...n (D equations)
      dD_equations_list <- lapply(seq_len(nstate), function(i) {
        # With i =/= j:
        # * psi[-i] is psi[j]
        # * q[i,][-i] is q[i, j]
        # * y[i + nstate] is Ei
        # * y[1:nstate][-i] is Dj
        # See derivation sxn in doi:10.1111/j.2041-210X.2012.00234.x for details
        return(
          -(lambda[i] + mu[i] + current_psi[i] + sum(q[i, ][-i])) * y[i]
          + 2 * lambda[i] * y[i + nstate] * y[i]
          + sum(q[i, ][-i] * y[1:nstate][-i])
        )
      })
      
      # States 1...n (E equations)
      dE_equations_list <- lapply(seq_len(nstate), function(i) {
        return(
          mu[i] - (lambda[i] + mu[i] + current_psi[i] + sum(q[i, ][-i]))
          * y[i + nstate] + lambda[i] * y[i + nstate]^2
          + sum(q[i][-i] * y[nstate + 1:nstate][-i])
        )
      })
      
      return(list(c(dD_equations_list, dE_equations_list)))
    })
  }
  
  # D1...Dn are child likelihoods, and E1...En are 1
  y <- c(child_likelihoods, rep(1, nstate))
  # Need to explicitly name index or events_df does not work
  names(y) <- seq_len(nstate * 2)
  
  times <- seq(0, tf , by = (tf) / 100)
  
  # Parameters
  parms <- list(
    lambda = params_df$lambda,
    mu = params_df$mu,
    current_psi=psi_values,
    q = q_matrix,
    nstate = nstate
  )
  
  # Events to set initial conditions
  events_df <- data.frame(
    var = seq_len(nstate),
    time = rep(t0),
    value = child_likelihoods,
    method = rep("replace", nstate)
  )
  suppressWarnings(
    sol <- deSolve::ode(y, times, func, parms, events = list(data = events_df))
  )
  
  ret <- utils::tail(sol, n = 1)[1 + 1:nstate]
  if (any(is.nan(ret))) {
    # if the value is nan, assign the likelihood to the previous likelihood,
    # because the value is too close so that the ode cannot tell the difference.
    return(child_likelihoods)
  }
  return(ret)
}

#' Backwards helper for time-dependent psi using chained constant segments
#'
#' @param child_likelihoods Initial child likelihood values
#' @param t0 Child age 
#' @param tf Parent age 
#' @param params_df Parameters data frame
#' @param q_matrix Transition rate matrix
#' @param psi_df Time-dependent psi data frame
#' @return Vector of likelihoods at parent
#' @noRd
backwards_likelihoods_helper_segmented <- function(child_likelihoods, t0, tf, params_df, q_matrix, psi_df) {
  nstate <- nrow(params_df)
  
  
  #TODO: sampling function, change the argument
  
  
  # Get intervals that overlap with this branch
  intervals <- get_branch_intervals(t0, tf, psi_df)

    if (nrow(intervals) == 0) {
    warning("No overlapping intervals found for branch, returning child likelihoods")
    return(as.numeric(child_likelihoods))
  }

  # Start with the initial child likelihoods
  current_likelihoods <- as.numeric(child_likelihoods)
  
  # Process each segment sequentially
  for (i in 1:nrow(intervals)) {
    segment_start <- intervals$segment_start[i]
    segment_end <- intervals$segment_end[i]
    segment_duration <- intervals$duration[i]
    psi_values <- intervals$psi_values[[i]]
    
    # Skip zero-duration segments
    if (segment_duration <= 1e-10) {
      next
    }
    
    # Create temporary params_df with constant psi for this segment
    temp_params_df <- params_df
    temp_params_df$psi <- psi_values
    
    # Call backwards_likelihoods_helper_constant for this segment
    # Use current_likelihoods as "child_likelihoods" for this segment
    # The segment goes from relative time 0 to segment_duration
    
    current_likelihoods <- backwards_likelihoods_helper_constant(
      child_likelihoods = current_likelihoods,
      t0 = 0,  # Start of this segment (relative time)
      tf = segment_duration,  # End of this segment (relative time)
      params_df = temp_params_df,
      q_matrix = q_matrix,
      psi_values = psi_values
    )
    
    # Validate result before continuing
    if (!is.numeric(current_likelihoods) || any(!is.finite(current_likelihoods))) {
      warning("Invalid result in segment ", i, ", returning previous result")
      break
    }
  }
  
  return(current_likelihoods)
}

#' Modified forwards likelihood helper that works for both cases
#'
#' @param parent_state_probabilities Parent state probabilities
#' @param t0 Parent age 
#' @param tf Child age 
#' @param params_df Parameters data frame  
#' @param q_matrix Transition rate matrix
#' @return Vector of likelihoods at child
#' @noRd
get_forwards_likelihoods <- function(parent_state_probabilities,
                                     t0, tf, params_df, q_matrix) {
  # Extract psi data and determine if it's time-dependent
  psi_info <- extract_psi_data(params_df)
  
  if (!psi_info$is_time_dependent) {
    # Constant sampling rate 
    return(get_forwards_likelihoods_constant(parent_state_probabilities, t0, tf, params_df, q_matrix, psi_info$psi_data))
  } else {
    # Time-dependent sampling rate
    return(get_forwards_likelihoods_segmented(parent_state_probabilities, t0, tf, params_df, q_matrix, psi_info$psi_data))
  }
}

#' Forwards helper for constant sampling rate
#'
#' @noRd
get_forwards_likelihoods_constant <- function(parent_state_probabilities, t0, tf, params_df, q_matrix, psi_values) {
  nstate <- nrow(params_df)
  
  func <- function(x, y, parms) {
    with(as.list(c(y, parms)), {
      current_psi <- psi_values
      
      # States 1...n (D equations - only transitions)
      dD_equations_list <- lapply(seq_len(nstate), function(i) {
        return(
          -(sum(q[i, ][-i]) * y[i])
          + sum(t(q)[i, ][-i] * y[1:nstate][-i])
        )
      })
      
      # States 1...n (E equations)
      dE_equations_list <- lapply(seq_len(nstate), function(i) {
        return(
          mu[i] - (lambda[i] + mu[i] + current_psi[i] + sum(q[i, ][-i]))
          * y[i + nstate] + lambda[i] * y[i + nstate]^2
          + sum(q[i][-i] * y[nstate + 1:nstate][-i])
        )
      })
      
      return(list(c(dD_equations_list, dE_equations_list)))
    })
  }
  
  # Initial conditions
  y <- c(parent_state_probabilities, rep(1, nstate))
  
  # Time setup
  times <- seq(0, t0, by = t0 / 100)
  # Parameters
  parms <- list(
    lambda = params_df$lambda,
    mu = params_df$mu,
    q = q_matrix,
    nstate = nstate
  )
  suppressWarnings(
    sol <- deSolve::ode(y, times, func, parms, method = "ode45", rtol = 1e-6)
  )
  
  # Get final values
  # Closest index to tf
  closest_index <- which.min(abs(sol[, 1] - tf))
  likelihoods <- unname(sol[closest_index, 1 + 1:nstate])
  
  return(likelihoods)
}

#' Forwards helper for time-dependent psi using chained constant segments
#'
#' @param parent_state_probabilities Initial parent state probabilities
#' @param t0 Parent age 
#' @param tf Child age 
#' @param params_df Parameters data frame
#' @param q_matrix Transition rate matrix
#' @param psi_df Time-dependent psi data frame
#' @return Vector of likelihoods at child
#' @noRd
get_forwards_likelihoods_segmented <- function(parent_state_probabilities, t0, tf, params_df, q_matrix, psi_df) {
  nstate <- nrow(params_df)
  
  # Get intervals (note: t0 > tf in forwards direction, so we swap them)
  intervals <- get_branch_intervals(tf, t0, psi_df)
  #print(intervals)
  if (nrow(intervals) == 0) {
    warning("No overlapping intervals found for forward branch")
    return(as.numeric(parent_state_probabilities))
  }
  
  # Start with the initial parent probabilities
  current_probabilities <- as.numeric(parent_state_probabilities)
  #print(current_probabilities)
  
  # Process segments in reverse order (from parent to child)
  # This is because intervals are ordered from tf to t0, but we want t0 to tf
  for (i in nrow(intervals):1) {
    segment_start <- intervals$segment_start[i]
    segment_end <- intervals$segment_end[i]
    segment_duration <- intervals$duration[i]
    psi_values <- intervals$psi_values[[i]]
    
    # Skip zero-duration segments
    if (segment_duration <= 1e-10) {
      next
    }
    
    # Create temporary params_df with constant psi for this segment
    temp_params_df <- params_df
    temp_params_df$psi <- psi_values
    
    # Call get_forwards_likelihoods_constant for this segment
    # Use current_probabilities as "parent_state_probabilities" for this segment
    current_probabilities <- get_forwards_likelihoods_constant(
      parent_state_probabilities = current_probabilities,
      t0 = segment_start,  # Start of this segment (relative time)
      tf = segment_end,  # End of this segment (relative time)
      params_df = temp_params_df,
      q_matrix = q_matrix,
      psi_values = psi_values
    )
    
    # Validate result before continuing
    if (!is.numeric(current_probabilities) || any(!is.finite(current_probabilities))) {
      warning("Invalid result in forward segment ", i, ", returning previous result")
      break
    }
  }
  #print(current_probabilities)
  return(current_probabilities)
}


#' Helper function to detect and extract psi data
#'
#' @param params_df Parameters data frame
#' @return List with: is_time_dependent (logical), psi_data (data.frame or numeric)
#' @noRd
extract_psi_data <- function(params_df) {
  psi_col <- params_df$psi
  
  # Case 1: Simple numeric vector (constant psi)
  if (is.numeric(psi_col)) {
    return(list(
      is_time_dependent = FALSE,
      psi_data = psi_col
    ))
  }
  
  # Case 2: List column with data frames (time-dependent psi)
  if (is.list(psi_col) && length(psi_col) > 0 && is.data.frame(psi_col[[1]])) {
    return(list(
      is_time_dependent = TRUE,
      psi_data = psi_col[[1]]  # Extract the data frame
    ))
  }
  
  # Case 3: Direct data frame assignment (time-dependent psi)
  if (is.data.frame(psi_col)) {
    return(list(
      is_time_dependent = TRUE,
      psi_data = psi_col
    ))
  }
  
  # Fallback: treat as constant
  warning("Unrecognized psi format, treating as constant")
  return(list(
    is_time_dependent = FALSE,
    psi_data = as.numeric(psi_col)
  ))
}

#' Get psi values at specific time using step function (constant within intervals)
#'
#' @param time_point Numeric value representing the time point
#' @param psi_df Data frame with time-dependent psi values. 
#'               First column: time, remaining columns: psi values for each state
#' @return Vector of psi values for each state at the given time point
#' @noRd
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

#' Get time intervals that overlap with a branch
#'
#' @param t_start Start time of branch (absolute time)
#' @param t_end End time of branch (absolute time)  
#' @param psi_df Data frame with time intervals
#' @return Data frame with overlapping intervals and their durations
#' @noRd
get_branch_intervals <- function(t_start, t_end, psi_df) {
  # Ensure t_start < t_end for clarity
  if (t_start > t_end) {
    temp <- t_start
    t_start <- t_end
    t_end <- temp
  }
  
  time_col <- psi_df[, 1]
  
  # Sort time points
  if (is.unsorted(time_col)) {
    sort_order <- order(time_col)
    time_col <- time_col[sort_order]
    psi_df <- psi_df[sort_order, ]
  }
  
  # Create interval boundaries
  # Intervals: [0, time1), [time1, time2), ..., [timeN, Inf)
  interval_starts <- c(0, time_col)
  interval_ends <- c(time_col, Inf)
  
  # Find which intervals overlap with [t_start, t_end]
  overlapping_intervals <- data.frame(
    interval_idx = integer(0),
    segment_start = numeric(0),
    segment_end = numeric(0),
    duration = numeric(0),
    psi_values = I(list())
  )
  
  for (i in seq_along(interval_starts)) {
    interval_start <- interval_starts[i]
    interval_end <- interval_ends[i]
    
    # Check if branch overlaps with this interval
    if (t_start < interval_end && t_end > interval_start) {
      # Calculate overlapping segment
      segment_start <- max(t_start, interval_start)
      segment_end <- min(t_end, interval_end)
      duration <- segment_end - segment_start
      
      # Get psi values for this interval
      if (i <= nrow(psi_df)) {
        psi_values <- as.numeric(psi_df[i, -1])
      } else {
        # Last interval - use last row of psi values
        psi_values <- as.numeric(psi_df[nrow(psi_df), -1])
      }
      
      # Add to result
      overlapping_intervals <- rbind(overlapping_intervals, data.frame(
        interval_idx = i,
        segment_start = segment_start,
        segment_end = segment_end,
        duration = duration,
        psi_values = I(list(psi_values))
      ))
    }
  }
  
  return(overlapping_intervals)
}