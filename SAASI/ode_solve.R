library("deSolve")

options(warn=-1)

get_backwards_likelihoods <- function(left_likelihoods, right_likelihoods,
                                      left_t0, right_t0, tf,
                                      params_df, q_matrix) {

  left_sol <- get_backwards_likelihoods_helper(left_likelihoods,
                                               left_t0, tf,
                                               params_df, q_matrix)

  right_sol <- get_backwards_likelihoods_helper(right_likelihoods,
                                                right_t0, tf,
                                                params_df, q_matrix)

  likelihoods <- 2* params_df$lambda * left_sol * right_sol
  return(likelihoods)
}

get_backwards_likelihoods_helper <- function(child_likelihoods, 
                                             t0, tf,
                                             params_df, q_matrix) {
  # Number of states == n
  nstate <- nrow(params_df)
  
  func <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      # States 1...n
      dD_equations_list <- lapply(seq_len(nstate), function(i) ({
        # With i =/= j:
        # * Ψ[-i] is Ψ[j]
        # * q[i,][-i] is q[i, j]
        # * y[i + nstate] is Ei
        # * y[1:nstate][-i] is Dj
        # See derivation sxn in doi:10.1111/j.2041-210X.2012.00234.x for details
        return(
          -(λ[i] + μ[i] + Ψ[i] + sum(q[i,][-i])) * y[i]
          + 2 * λ[i] * y[i + nstate] * y[i]
          + sum(q[i,][-i] * y[1:nstate][-i])
        )
      }))
      
      # States 1...n
      dE_equations_list <- lapply(seq_len(nstate), function(i) ({
        # With i =/= j:
        # * Ψ[-i] is Ψ[j]
        # * q[i,][-i] is q[i, j]
        # * y[i + nstate] is Ei
        # * y[nstate + 1:nstate][-i] is Ej
        # See derivation sxn in doi:10.1111/j.2041-210X.2012.00234.x for details
        return(
          μ[i] - (λ[i] + μ[i] + Ψ[i]+ sum(q[i,][-i])) * y[i + nstate]
          + λ[i] * y[i + nstate]^2
          + sum(q[i][-i] * y[nstate + 1:nstate][-i])
        )
      }))
      
      return(list(c(dD_equations_list, dE_equations_list)))
    })
  }
  
  
  # D1...Dn are NA, and E1...En are 1
  y <- c(child_likelihoods, rep(1, nstate))
  
  # Need to explicitly name index or events_df does not work
  names(y) <- seq_len(nstate * 2)
  
  times <- seq(0, tf, by = tf / 100)
  parms <- list(λ = params_df$lambda,
                μ = params_df$mu,
                Ψ = params_df$psi,
                q = q_matrix,
                nstate = nstate)
  #print("y")
  #print(y)
  # Force D1...Dn at t0 to be same as children
  events_df <- data.frame(var = seq_len(nstate),
                          time = rep(t0),
                          value = child_likelihoods,
                          method = rep("replace", nstate))
  #browser()
  # Suppress warnings about t0 not in times
  # suppressWarnings(
  #   sol <- ode(y, times, func, parms, events = list(data = events_df),method="ode45",atol = 0, rtol = 1e-14)
  # )
  #print(y)
  #print("events_df")
  #print(events_df)
  suppressWarnings(
    sol <- ode(y, times, func, parms, events = list(data = events_df),atol=1e-2,rtol=1e-2)
  )
  if(is.nan((tail(sol, n = 1)[1 + 1:nstate][1]))){
    print("the value is nan")
    # if the value is nan, assign the likelihood to the previous likelihood,
    # because the value is too close so that the ode cannot tell the difference.
    #return(child_likelihoods)
    #sol <- ode(y, times, func, parms,method="ode45", atol= 1e-6,rtol = 1e-6)
  }

  if(is.nan((tail(sol, n = 1)[1 + 1:nstate][1]))){
    return(child_likelihoods)
  }
  else{return(tail(sol, n = 1)[1 + 1:nstate])}
}

get_forwards_likelihoods <- function(parent_state_probabilities, t0, tf,
                                     params_df, q_matrix) {
  # Number of states == n
  nstate <- nrow(params_df)
  
  func <- function(x, y, parms) {
    with(as.list(c(y, parms)), {
      # States 1...n
      dD_equations_list <- lapply(seq_len(nstate), function(i) ({
        
        return(
          -(sum(q[i,][-i]))* y[i]
          + sum(t(q)[i,][-i] * y[1:nstate][-i])
        )
      }))
      
      
      
      # States 1...n
      dE_equations_list <- lapply(seq_len(nstate), function(i) ({
        return(
          μ[i] - (λ[i] + μ[i] +Ψ[i]+ sum( q[i,][-i])) * y[i + nstate]
          + λ[i] * y[i + nstate]^2
          + sum(q[i][-i] * y[nstate + 1:nstate][-i])
        )
      }))
      
      return(list(c(dD_equations_list, dE_equations_list)))
    })
  }
  yini <- c(parent_state_probabilities, rep(1, nstate))
  
  # Increment time in the positive direction because otherwise the ode solver
  # can run into errors with negative numbers being smaller than machine min.
  times <- seq(0, t0, by = t0 / 100)
  parms <- list(λ = params_df$lambda,
                μ = params_df$mu,
                Ψ = params_df$psi,
                q = q_matrix,
                nstate = nstate)
  
  # Suppress warnings about initial conditions guessed as 0
  #print("doing forward")
  suppressWarnings(
    sol <- ode(yini, times, func, parms,method="ode45", rtol = 1e-20)
  )
  
  # Closest index to tf
  closest_index <- which.min(abs(sol[, 1] - tf))
  #print(closest_index)
  likelihoods <- unname(sol[closest_index, 1 + 1:nstate])
  #print(sum(likelihoods))
  #return(likelihoods / sum(likelihoods))
  return(likelihoods)
}

