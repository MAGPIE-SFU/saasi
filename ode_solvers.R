library("deSolve")
library("bvpSolve", warn.conflicts = FALSE)

get_backwards_sol <- function(backwards_state_1,
                              backwards_state_2,
                              t0, tf,
                              λ, μ, Ψ, q, id) {
  func <- function(t, y, parms) {
    dd1 <- (-(λ + μ + Ψ[2] + q[[1, 2]]) * y["d1"]
            + 2 * λ * y["e1"] * y["d1"] + q[[1, 2]] * y["d2"])
    dd2 <- (-(λ + μ + Ψ[1] + q[[2, 1]]) * y["d2"]
            + 2 * λ * y["e2"] * y["d2"] + q[[2, 1]] * y["d1"])

    de1 <- (μ - (λ + μ + Ψ[2] + q[[1, 2]]) * y["e1"]
            + λ * y["e1"]^2 + q[[1, 2]] * y["e2"])
    de2 <- (μ - (λ + μ + Ψ[1] + q[[2, 1]]) * y["e2"]
            + λ * y["e2"]^2 + q[[2, 1]] * y["e1"])

    return(list(c(dd1, dd2, de1, de2)))
  }

  y <- c(d1 = NA,
         d2 = NA,
         e1 = 1,
         e2 = 1)
  names(y) <- c("d1", "d2", "e1", "e2")
  times <- seq(0, tf, by = tf / 100)
  parms <- c(λ = λ, μ = μ, Ψ = Ψ, q = q)

  events_df <- data.frame(var = c("d1", "d2"),
                          time = c(t0, t0),
                          value = c(backwards_state_1, backwards_state_2),
                          method = c("replace", "replace"))

  # Suppress warnings about t0 not in times
  suppressWarnings(
    sol <- ode(y, times, func, parms, events = list(data = events_df))
  )

  # TODO Is this right? Not technically at tf
  tf_d1 <- tail(sol[, "d1"], n = 1)
  tf_d2 <- tail(sol[, "d2"], n = 1)

  return(c(backwards_state_1 = tf_d1, backwards_state_2 = tf_d2))
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
