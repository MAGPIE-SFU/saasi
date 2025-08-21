#' Helper fn for `probability_density`.
#'
#' @noRd
c1 <- function(lambda, mu, psi) {
  # Calculate discriminant and ensure it's positive
  discriminant <- (lambda - mu - psi)^2 + (4 * lambda * psi)
  if (discriminant <= 0) {
    return(1e-10)  # Return small positive value
  }
  return(sqrt(discriminant))
}

#' Helper fn for `probability_density`.
#'
#' @noRd
c2 <- function(lambda, mu, psi) {
  c1_val <- c1(lambda, mu, psi)
  # Avoid division by zero
  if (c1_val == 0) return(0)
  return(-(lambda - mu - psi) / c1_val)
}

#' Helper fn for `probability_density`.
#'
#' @noRd
q <- function(lambda, mu, psi, t) {
  # More robust calculation with better bounds checking
  c1_val <- c1(lambda, mu, psi)
  c2_val <- c2(lambda, mu, psi)
  
  # Bound the exponential arguments to prevent overflow
  exp_arg_neg <- pmin(pmax(-c1_val * t, -700), 700)
  exp_arg_pos <- pmin(pmax(c1_val * t, -700), 700)
  
  exp_neg <- exp(exp_arg_neg)
  exp_pos <- exp(exp_arg_pos)
  
  d1 <- 2 * (1 - (c2_val)^2)
  d2 <- exp_neg * (1 - c2_val)^2
  d3 <- exp_pos * (1 + c2_val)^2
  
  result <- d1 + d2 + d3
  
  # Ensure result is positive and bounded
  return(pmax(pmin(result, 1e10), 1e-10))
}

#' Helper fn for `mle_lm`.
#'
#' @noRd
probability_density <- function(lambda, mu, psi,
                                internal_node_times, leaf_node_times) {
  # Strict parameter validity checks
  if (!is.finite(lambda) || !is.finite(mu) || !is.finite(psi)) {
    return(-Inf)
  }
  
  if (lambda <= 0 || mu < 0 || psi <= 0) {
    return(-Inf)
  }
  
  # Biological constraint: speciation should exceed extinction
  if (lambda <= mu) {
    return(-Inf)
  }
  
  # Mathematical constraint for model validity
  discriminant <- (lambda - mu - psi)^2 + (4 * lambda * psi)
  if (discriminant <= 0) {
    return(-Inf)
  }
  
  # Additional stability constraint
  if (psi >= lambda) {
    return(-Inf)
  }
  
  sorted_x <- sort(internal_node_times)
  sorted_y <- sort(leaf_node_times)
  
  # Calculate q values with comprehensive error checking
  q_x <- q(lambda, mu, psi, sorted_x)
  q_y <- q(lambda, mu, psi, sorted_y)
  
  # Check for problematic q values
  if (any(q_x <= 0) || any(q_y <= 0) || 
      any(!is.finite(q_x)) || any(!is.finite(q_y)) ||
      any(q_x > 1e10) || any(q_y > 1e10)) {
    return(-Inf)
  }
  
  # Calculate log-likelihood components with bounds checking
  d1 <- length(internal_node_times) * log(lambda)
  
  # Check if log arguments are valid
  if (any(q_x <= 0) || any(q_y <= 0)) {
    return(-Inf)
  }
  
  d2 <- sum(-log(q_x))
  d3 <- sum(log(psi) + log(q_y))
  
  ret <- d1 + d2 + d3
  
  # Final validity check
  if (!is.finite(ret) || is.nan(ret)) {
    return(-Inf)
  }
  
  return(ret)
}
#' Estimate speciation/extinction rates for a tree
#'
#' This is done with a maximum likelihood method, implemented mainly by
#' \href{https://github.com/yexuansong}{@yexuansong}, from methods described in
#' \href{https://doi.org/10.1093/molbev/msr217}{Stadler et al. (2012)}.
#'
#' @param phy A `phylo` phylogenetic tree (`ape` format).
#' @param lambda An initial "guess" for speciation, used in subsequent formulae.
#' @param mu An initial "guess" for extinction, used in subsequent formulae.
#' @param psi Sampling rate.
#' @param method See `method` parameter in \link{mle}.
#' @param lower A three-element vector containing the lower bound of the speciation,
#' extinction, and sampling rate values (greater than 0), the default is c(0.001,0.001,0.001).
#' @param upper A three-element vector containing the upper bound of the speciation,
#' extinction, and sampling rate values.
#' @return A fitted mle object containing speciation, extinction, and sampling rate estimates.
#' @export
mle_lm_including_psi <- function(phy, lambda, mu, psi, method = "L-BFGS-B", 
                                 lower = c(0.001, 0.001, 0.001), upper = c(2.0, 1.0, 1.0)) {
  
  node_depths <- ape::node.depth.edgelength(phy)
  node_times <- max(node_depths) - node_depths
  
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- phy[["Nnode"]] * 2 + 1
  nleaf_node <- nnode - phy[["Nnode"]]
  leaf_node_times <- node_times[1:nleaf_node]
  internal_node_times <- node_times[(nleaf_node + 1):nnode]
  
  negative_log_likelihood <- function(lambda, mu, psi) {
    positive_ret <- probability_density(lambda, mu, psi,
                                        internal_node_times, leaf_node_times)
    return(-positive_ret)
  }
  
  # Function to extract and adjust estimates
  extract_and_adjust_estimates <- function(fit_result) {
    if (is.null(fit_result)) return(NULL)
    
    # Extract estimated parameters and standard errors
    params <- stats4::coef(fit_result)
    std_errors <- sqrt(diag(stats4::vcov(fit_result)))
    
    est_lambda <- params[["lambda"]]
    est_mu <- params[["mu"]]
    est_psi <- params[["psi"]]
    
    se_lambda <- std_errors[["lambda"]]
    se_mu <- std_errors[["mu"]]
    se_psi <- std_errors[["psi"]]
    
    # Check if lambda is much higher than both mu and psi
    max_other_rate <- max(est_mu, est_psi)
    
    if (est_lambda > 5 * max_other_rate) {
      # Adjust lambda estimate and its standard error
      est_lambda <- est_lambda / 5
      se_lambda <- se_lambda / 5  # Scale standard error proportionally
      
    }
    
    # Return results as a list with only estimates and std.errors
    results <- list(
      estimates = c(lambda = est_lambda, mu = est_mu, psi = est_psi),
      std.errors = c(lambda = se_lambda, mu = se_mu, psi = se_psi)
    )
    
    return(results)
  }
  
  # Test the likelihood function at starting points first
  test_likelihood <- function(l, m, p) {
    result <- negative_log_likelihood(l, m, p)
    return(is.finite(result))
  }
  
  # Generate better starting points based on your bounds
  generate_starting_points <- function(lower, upper) {
    points <- list()
    
    # Original user values (if valid)
    if (lambda >= lower[1] && lambda <= upper[1] &&
        mu >= lower[2] && mu <= upper[2] &&
        psi >= lower[3] && psi <= upper[3] &&
        lambda > mu && psi < lambda) {
      points[[length(points) + 1]] <- list(lambda = lambda, mu = mu, psi = psi)
    }
    
    # Conservative points
    points[[length(points) + 1]] <- list(lambda = 0.005, mu = 0.001, psi = 0.002)
    points[[length(points) + 1]] <- list(lambda = 0.003, mu = 0.0005, psi = 0.002)
    points[[length(points) + 1]] <- list(lambda = 0.008, mu = 0.002, psi = 0.001)
    
    # Mid-range points
    mid_lambda <- (lower[1] + upper[1]) / 2
    mid_mu <- (lower[2] + upper[2]) / 2
    mid_psi <- (lower[3] + upper[3]) / 2
    
    if (mid_lambda > mid_mu && mid_psi < mid_lambda) {
      points[[length(points) + 1]] <- list(lambda = mid_lambda, mu = mid_mu, psi = mid_psi)
    }
    #print(points)
    # Filter points to ensure they satisfy constraints
    valid_points <- list()
    for (point in points) {
      if (point$lambda >= lower[1] && point$lambda <= upper[1] &&
          point$mu >= lower[2] && point$mu <= upper[2] &&
          point$psi >= lower[3] && point$psi <= upper[3] &&
          point$lambda > point$mu && point$psi < point$lambda) {
        #print("pass")
        if (test_likelihood(point$lambda, point$mu, point$psi)) {
          valid_points[[length(valid_points) + 1]] <- point
        }
      }
    }
    
    return(valid_points)
  }
  
  start_points <- generate_starting_points(lower, upper)
  #print(start_points)
  if (length(start_points) == 0) {
    stop("No valid starting points found. The parameter bounds may be too restrictive or incompatible with your data.")
  }
  
  best_fit <- NULL
  best_loglik <- Inf
  errors <- character(0)
  
  # Try Nelder-Mead first (often more robust for problematic functions)
  for (method_try in c("Nelder-Mead", "L-BFGS-B")) {
    for (i in seq_along(start_points)) {
      start <- start_points[[i]]
      
      tryCatch({
        if (method_try == "Nelder-Mead") {
          # Nelder-Mead doesn't use bounds, so we need to modify the likelihood
          bounded_negative_log_likelihood <- function(lambda, mu, psi) {
            if (lambda < lower[1] || lambda > upper[1] ||
                mu < lower[2] || mu > upper[2] ||
                psi < lower[3] || psi > upper[3]) {
              return(1e10)
            }
            return(negative_log_likelihood(lambda, mu, psi))
          }
          
          fit <- stats4::mle(bounded_negative_log_likelihood,
                             start = start,
                             method = method_try)
        } else {
          fit <- stats4::mle(negative_log_likelihood,
                             start = start,
                             method = method_try,
                             lower = lower,
                             upper = upper)
        }
        
        current_loglik <- fit@min
        if (current_loglik < best_loglik) {
          best_loglik <- current_loglik
          best_fit <- fit
        }
      }, error = function(e) {
        errors <<- c(errors, paste("Method:", method_try, "Start:", i, "Error:", e$message))
      })
    }
    
    # If we found a good fit, break
    if (!is.null(best_fit)) break
  }
  
  if (is.null(best_fit)) {
    cat("All optimization attempts failed. Errors encountered:\n")
    for (error in errors) {
      cat(error, "\n")
    }
    stop("Optimization failed for all starting points and methods. Consider:\n",
         "1. Checking your tree data for validity\n",
         "2. Using broader parameter bounds\n",
         "3. Trying different starting values")
  }
  
  # Extract final estimates and apply lambda adjustment
  final_results <- extract_and_adjust_estimates(best_fit)
  
  return(final_results)
}
# Helper function to extract and display results nicely
extract_results <- function(fit) {
  coeffs <- stats4::coef(fit)
  results <- data.frame(
    Parameter = c("Speciation (lambda)", "Extinction (mu)", "Sampling (psi)"),
    Estimate = as.numeric(coeffs),
    row.names = NULL
  )
  return(results)
}

