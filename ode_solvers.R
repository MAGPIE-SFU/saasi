library("deSolve")
library("bvpSolve", warn.conflicts = FALSE)

get_backwards_likelihoods <- function(left_likelihoods, right_likelihoods,
                                      left_t0, right_t0, tf,
                                      params_df, q_matrix) {
  left_sol <- get_backwards_likelihoods_helper(left_likelihoods,
                                               left_t0, tf,
                                               params_df, q_matrix)
  right_sol <- get_backwards_likelihoods_helper(right_likelihoods,
                                                right_t0, tf,
                                                params_df, q_matrix)
  return(params_df$lambda * left_sol * right_sol)
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
          -(λ[i] + μ[i] + sum(Ψ[-i] + q[i,][-i])) * y[i]
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
          μ[i] - (λ[i] + μ[i] + sum(Ψ[-i] + q[i,][-i])) * y[i + nstate]
          + λ[i] * y[i + nstate]^2
          + sum(q[i][-i] * y[nstate + 1:nstate][-i])
        )
      }))

      return(list(c(dD_equations_list, dE_equations_list)))
    })
  }
  

  # D1...Dn are NA, and E1...En are 1
  y <- c(rep(NA, nstate), rep(1, nstate))
  # Need to explicitly name index or events_df does not work
  names(y) <- seq_len(nstate * 2)

  times <- seq(0, tf, by = tf / 100)
  parms <- list(λ = params_df$lambda,
                μ = params_df$mu,
                Ψ = params_df$psi,
                q = q_matrix,
                nstate = nstate)

  # Force D1...Dn at t0 to be same as children
  events_df <- data.frame(var = seq_len(nstate),
                          time = rep(t0),
                          value = child_likelihoods,
                          method = rep("replace", nstate))
  # browser()

  # Suppress warnings about t0 not in times
  suppressWarnings(
    sol <- ode(y, times, func, parms, events = list(data = events_df))
  )

  return(tail(sol, n = 1)[1 + 1:nstate])
}

get_state_probabilities <- function(parent_state_probabilities, t0, tf,
                                    params_df, q_matrix) {
  # Number of states == n
  nstate <- nrow(params_df)

  func <- function(x, y, parms) {
    with(as.list(c(y, parms)), {
      # States 1...n
      dD_equations_list <- lapply(seq_len(nstate), function(i) ({
        # With i =/= j:
        # * Ψ[-i] is Ψ[j]
        # * q[i,][-i] is q[i, j]
        # * y[i + nstate] is Ei
        # * y[1:nstate][-i] is Dj
        # See derivation sxn in doi:10.1111/j.2041-210X.2012.00234.x for details
        return(-(
          -(λ[i] + μ[i] + sum(Ψ[-i] + q[i,][-i])) * y[i]
          + 2 * λ[i] * y[i + nstate] * y[i]
          + sum(q[i,][-i] * y[1:nstate][-i])
        ))
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
          μ[i] - (λ[i] + μ[i] + sum(Ψ[-i] + q[i,][-i])) * y[i + nstate]
          + λ[i] * y[i + nstate]^2
          + sum(q[i][-i] * y[nstate + 1:nstate][-i])
        )
      }))

      return(list(c(dD_equations_list, dE_equations_list)))
    })
  }

  # D1...Dn are parent state probabilities, and E1...En are NA
  yini <- c(parent_state_probabilities, rep(NA, nstate))
  # D1...Dn are NA, and E1...En are 0
  yend <- c(rep(NA, nstate), rep(0, nstate))

  # Increment time in the positive direction because otherwise the ode solver
  # can run into errors with negative numbers being smaller than machine min.
  x <- seq(0, t0, by = t0 / 100)
  parms <- list(λ = params_df$lambda,
                μ = params_df$mu,
                Ψ = params_df$psi,
                q = q_matrix,
                nstate = nstate)

  # Suppress warnings about initial conditions guessed as 0
  suppressWarnings(
    # Run opposite directions because of positively increasing x. Should not
    # affect result.
    sol <- bvpshoot(yend, x, func, yini, parms)
  )

  # Closest index to tf
  closest_index <- which.min(abs(sol[, "x"] - tf))

  return(unname(sol[closest_index, 1 + 1:nstate]))
}

get_forwards_sol <- function(ancestral_state_1,
                             ancestral_state_2,
                             t0, tf,
                             λ, μ, Ψ, q, id) {
  func <- function(x, y, parms) {
    dd1 <- -(-(λ + μ + Ψ[2] + q[[1, 2]]) * y["d1"]
             + 2 * λ * y["e1"] * y["d1"] + q[[1, 2]] * y["d2"])
    dd2 <- -(-(λ + μ + Ψ[1] + q[[2, 1]]) * y["d2"]
             + 2 * λ * y["e2"] * y["d2"] + q[[2, 1]] * y["d1"])

    de1 <- (μ - (λ + μ + Ψ[2] + q[[1, 2]]) * y["e1"]
            + λ * y["e1"]^2 + q[[1, 2]] * y["e2"])
    de2 <- (μ - (λ + μ + Ψ[1] + q[[2, 1]]) * y["e2"]
            + λ * y["e2"]^2 + q[[2, 1]] * y["e1"])

    return(list(c(dd1, dd2, de1, de2)))
  }

  yini <- c(d1 = ancestral_state_1,
            d2 = ancestral_state_2,
            e1 = NA,
            e2 = NA)
  yend <- c(d1 = NA,
            d2 = NA,
            e1 = 0,
            e2 = 0)
  names(yini) <- c("d1", "d2", "e1", "e2")
  names(yend) <- c("d1", "d2", "e1", "e2")
  # Increment time in the positive direction because otherwise the ode solver
  # can run into errors with negative numbers being smaller than machine min.
  x <- seq(0, t0, by = t0 / 100)
  parms <- c(λ = λ, μ = μ, Ψ = Ψ, q = q)

  # Suppress warnings about initial conditions guessed as 0
  # TODO correct guess? Yexuan's program in Mathematica seems to guess 0 as well
  suppressWarnings(
    # Run opposite directions because of positively increasing x. Should not
    # affect result.
    sol <- bvpshoot(yend, x, func, yini, parms)
  )

  # TODO Is this right? Not technically at tf
  # Closest index to tf
  closest_index <- which.min(abs(sol[, "x"] - tf))
  tf_d1 <- sol[[closest_index, "d1"]]
  tf_d2 <- sol[[closest_index, "d2"]]

  return(c(forwards_state_1 = tf_d1, forwards_state_2 = tf_d2))
}
