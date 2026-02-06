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
        # See derivation sxn in doi:10.1111/j.2041-210X.2012.00234.x for details
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
  
  # Suppress warnings about t0 not in times
  suppressWarnings(
    sol <- deSolve::ode(y, times, func, parms, events = list(data = events_df),
                        atol = 1e-15, rtol = 1e-12)
  )
  
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
    sol <- deSolve::ode(y, times, func, parms, method = "ode45", rtol = 1e20)
  )
  
  # Closest index to tf
  closest_index <- which.min(abs(sol[, 1] - tf))
  likelihoods <- unname(sol[closest_index, 1 + 1:nstate])
  
  return(likelihoods)
}