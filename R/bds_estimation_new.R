c1 <- function(lambda, mu, psi){
  abs(sqrt((lambda - mu - psi)^2 + 4 * lambda * psi))
}

c2 <- function(lambda, mu, psi){
  (-lambda + mu + psi) / c1(lambda, mu, psi)
}

logsumexp <- function(a, b, c){
  m <- pmax(a, pmax(b, c))
  m + log(exp(a - m) + exp(b - m) + exp(c - m))
}

log_q <- function(t, lambda, mu, psi){
  c1v <- c1(lambda, mu, psi)
  c2v <- c2(lambda, mu, psi)
  
  term <- 2 * (1 - c2v^2)
  if(term <= 0){
    return(rep(-Inf, length(t)))
  } 
  A <- log(term)
  B <- (-c1v * t) + 2 * log(abs(1 - c2v))
  C <- ( c1v * t) + 2 * log(abs(1 + c2v))
  logsumexp(A, B, C)
}

log_ge <- function(phy, lambda, mu, psi){
  node_depths <- ape::node.depth.edgelength(phy)
  node_times <- max(node_depths) - node_depths
  nnode <- phy[["Nnode"]] * 2 + 1
  nleaf_node <- nnode - phy[["Nnode"]]
  leaf_node_times <- node_times[1:nleaf_node]
  internal_node_times <- node_times[(nleaf_node + 1):nnode]
  log_lik <- phy[["Nnode"]] * log(lambda) +
    (phy[["Nnode"]] + 1) * log(psi) +
    sum(log_q(leaf_node_times, lambda, mu, psi)) -
    sum(log_q(internal_node_times, lambda, mu, psi))
  return(log_lik)
}

estimate_a <- function(phy, trim = c(0.10, 0.50)){
  
  if(length(trim) != 2 || !is.numeric(trim) || any(trim < 0) || any(trim > 1) || trim[1] >= trim[2]){
    stop("'trim' must be a numeric vector of length 2, and 0 <= trim[1] < trim[2] <= 1.")
  }
  
  phy <- ape::reorder.phylo(phy, "cladewise")
  h <- ape::node.depth.edgelength(phy)
  ntip <- ape::Ntip(phy)
  internal <- (ntip + 1):(ntip + phy$Nnode)
  bt <- sort(h[internal])
  bt <- bt[bt > 0]
  L <- 1 + seq_along(bt)
  
  tmin <- quantile(bt, trim[1])
  tmax <- quantile(bt, trim[2])
  keep <- bt >= tmin & bt <= tmax
  fit <- lm(log(L[keep]) ~ bt[keep])
  unname(coef(fit)[2])
}

generate_starting_points <- function(
    n_starts,
    mu,
    psi_max = 7,
    psi_min = 0,
    r0_min = 1,
    r0_max = 5,
    removal_rate_min = 0,
    removal_rate_max = Inf
){
  starts <- list()
  set.seed(123)
  
  # Identify lambda_min
  lambda_min_from_r0 <- r0_min * removal_rate_min
  lambda_min_from_growth <- removal_rate_min + 0.001
  lambda_min <- max(lambda_min_from_r0, lambda_min_from_growth)
  
  # Identify lambda_max
  if(is.finite(removal_rate_max)){
    lambda_max <- r0_max * removal_rate_max
  }else{
    lambda_max <- r0_max * 100
  }
  
  if(lambda_min >= lambda_max){
    warning("Cannot generate valid lambda values. Lambda bounds depend on removal_rate_min and removal_rate_max.")
    return(list())
  }
  
  max_attempts <- n_starts * 100
  attempt <- 0
  while(length(starts) < n_starts && attempt < max_attempts){
    attempt <- attempt + 1
    lambda_ini <- runif(1, lambda_min, lambda_max)
    psi_ini <- runif(1, psi_min, psi_max)
    
    removal_rate <- mu + psi_ini
    r0_ini <- lambda_ini / removal_rate
    
    valid <- TRUE
    # Removal rate bounds
    if(removal_rate < removal_rate_min || removal_rate > removal_rate_max){
      valid <- FALSE
    }
    # Lambda > mu + psi (R0 > 1)
    if(lambda_ini <= removal_rate){
      valid <- FALSE
    }
    # R0 bounds
    if(r0_ini < r0_min || r0_ini > r0_max){
      valid <- FALSE
    }
    # psi bounds
    if(psi_ini < psi_min || psi_ini > psi_max){
      valid <- FALSE
    }
    if(valid){
      starts[[length(starts) + 1]] <- c(lambda = lambda_ini, psi = psi_ini)
    }
  }
  
  if(length(starts) < n_starts){
    warning("Could not generate all requested starting points. Consider relaxing constraints.")
  }
  return(starts)
}

# Estimating lambda and psi given fixed mu
fit_bd_fixed_mu <- function(
    phy,
    mu,
    n_starts = 10,
    psi_max = 7,
    r0_min = 1,
    r0_max = 5,
    infectious_period_min = NULL,
    infectious_period_max = NULL
){
  if(!is.null(infectious_period_max)){
    removal_rate_min <- 1/infectious_period_max
  }else{
    removal_rate_min <- 0
  }
  if(!is.null(infectious_period_min)){
    removal_rate_max <- 1/infectious_period_min
  }else{
    removal_rate_max <- Inf
  }
  
  psi_min <- max(0, removal_rate_min - mu)
  if(is.finite(removal_rate_max)){
    psi_max_from_period <- removal_rate_max - mu
  }else{
    psi_max_from_period <- Inf
  }
  psi_max_eff <- min(psi_max, psi_max_from_period)
  
  if(psi_max_eff < psi_min){
    stop(sprintf(
      "Infeasible constraints: psi must be in [%.4f, %.4f] which is empty. mu=%.4f is incompatible with infectious_period_max=%.4f. Either decrease mu or increase infectious_period_max.", 
      psi_min, psi_max_eff, mu, infectious_period_max
    ))
  }
  
  lambda_min_from_r0 <- r0_min * removal_rate_min
  lambda_min_from_growth <- removal_rate_min + 1e-6
  lambda_min <- max(lambda_min_from_r0, lambda_min_from_growth)
  
  if(is.finite(removal_rate_max)){
    lambda_max <- r0_max * removal_rate_max
  }else{
    lambda_max <- r0_max * 1000
  }
  
  if(lambda_min >= lambda_max){
    stop(sprintf(
      "Infeasible constraints: lambda must be in [%.4f, %.4f] which is empty. Check that r0_min (%.2f), r0_max (%.2f), and infectious_period constraints are compatible.",
      lambda_min, lambda_max, r0_min, r0_max
    ))
  }
  
  lower <- c(lambda_min, psi_min)
  upper <- c(lambda_max, psi_max_eff)
  
  starts <- generate_starting_points(
    n_starts = n_starts,
    mu = mu,
    psi_max = psi_max_eff,
    psi_min = psi_min,
    r0_min = r0_min,
    r0_max = r0_max,
    removal_rate_min = removal_rate_min,
    removal_rate_max = removal_rate_max
  )
  
  if(length(starts) == 0){
    stop("Could not generate any valid starting points. Try relaxing the constraints.")
  }
  
  # Likelihood function with quadratic penalty
  loglik <- function(par){
    lambda <- par[1]
    psi    <- par[2]
    
    if(lambda <= 0 || psi < 0) return(-1e10)
    
    r0 <- lambda / (mu + psi)
    
    # Compute base likelihood
    base_lik <- tryCatch({
      log_ge(phy, lambda, mu, psi)
    }, error = function(e){
      return(-1e10)
    })
    
    if(!is.finite(base_lik)){
      return(base_lik)
    }
    
    # Apply penalty for R0 constraint violations
    penalty <- 0
    if(r0 < r0_min){
      violation <- r0_min - r0
      penalty <- 1000 * violation^2
    }
    if(r0 > r0_max){
      violation <- r0 - r0_max
      penalty <- 1000 * violation^2
    }
    
    return(base_lik - penalty)
  }
  
  results <- list()
  log_liks <- numeric(length(starts))
  
  for(i in seq_along(starts)){
    start_i <- starts[[i]]
    fit <- tryCatch({
      optim(
        par     = start_i,
        fn      = loglik,
        method  = "L-BFGS-B",
        lower   = lower,
        upper   = upper,
        control = list(fnscale = -1, maxit = 2000)
      )
    }, error = function(e){
      list(value = -1e10, convergence = 999, par = start_i, message = as.character(e))
    })
    
    results[[i]] <- fit
    log_liks[i] <- fit$value
  }
  
  # Find best result
  best_idx <- which.max(log_liks)
  best_fit <- results[[best_idx]]
  
  # Check if we got reasonable results
  if(!is.finite(best_fit$value) || best_fit$value < -1e9){
    stop("All optimization attempts failed. Try relaxing constraints.")
  }
  
  list(
    lambda             = best_fit$par[1],
    psi                = best_fit$par[2],
    mu                 = mu,
    r0                 = best_fit$par[1] / (mu + best_fit$par[2]),
    infectious_period  = 1 / (mu + best_fit$par[2]),
    logLik             = best_fit$value,
    convergence        = best_fit$convergence,
    constraints        = list(
      psi_max = psi_max_eff,
      psi_min = psi_min,
      r0_min = r0_min,
      r0_max = r0_max,
      infectious_period_min = infectious_period_min,
      infectious_period_max = infectious_period_max,
      lambda_min = lambda_min,
      lambda_max = lambda_max
    )
  )
}

# Main estimation function
estimate_bds_parameters <- function(
    phy,
    mu,
    trim = c(0.10, 0.50),
    n_starts = 10,
    force_two_step = FALSE,
    psi_max = 7,
    r0_min = 1,
    r0_max = 5,
    infectious_period_min = NULL,
    infectious_period_max = NULL
){
  # Input validation
  if(missing(mu) || is.null(mu)){
    stop("Extinction rate 'mu' is required.")
  }
  
  if(!inherits(phy, "phylo")){
    stop("Input 'phy' must be a phylo object.")
  }
  
  # Early constraint compatibility validation
  if(!is.null(infectious_period_min) && !is.null(infectious_period_max)){
    removal_rate_min <- 1/infectious_period_max
    removal_rate_max <- 1/infectious_period_min
    
    # Check mu compatibility
    if(mu >= removal_rate_max){
      stop(sprintf(
        "mu (%.4f) >= 1/infectious_period_min (%.4f). Impossible to satisfy infectious period constraints. Either decrease mu or increase infectious_period_min.",
        mu, removal_rate_max
      ))
    }
    
    # Check lambda compatibility
    min_required_lambda <- r0_min * removal_rate_min
    max_allowed_lambda <- r0_max * removal_rate_max
    
    if(min_required_lambda > max_allowed_lambda){
      stop(sprintf(
        "Conflicting R0 constraints: r0_min requires lambda >= %.4f, but r0_max allows lambda <= %.4f. Either decrease r0_min or increase r0_max.",
        min_required_lambda, max_allowed_lambda
      ))
    }
    
    # Check psi range
    psi_min <- max(0, removal_rate_min - mu)
    psi_max_eff <- min(psi_max, removal_rate_max - mu)
    
    if(psi_min > psi_max_eff){
      stop(sprintf(
        "Infeasible psi range: [%.4f, %.4f]. mu=%.4f is incompatible with infectious period constraints. Either decrease mu or adjust infectious_period bounds.",
        psi_min, psi_max_eff, mu
      ))
    }
  }
  
  # Initialize variables
  two_step_successful <- FALSE
  estimated_psi <- NA
  estimated_lambda <- NA
  a <- NA
  b <- NA
  fixed_mu_fit <- NULL
  
  tryCatch({
    # Step 1: Fit birth and sampling rates given fixed mu (MLE)
    fixed_mu_fit <- fit_bd_fixed_mu(
      phy,
      mu = mu,
      n_starts = n_starts,
      psi_max = psi_max,
      r0_min = r0_min,
      r0_max = r0_max,
      infectious_period_min = infectious_period_min,
      infectious_period_max = infectious_period_max
    )
    
    # Step 2: Estimate net rate from lineage-through-time data
    a <- estimate_a(phy, trim = trim)
    
    # Step 3: Calculate ratio b = lambda/(mu+psi)
    b <- fixed_mu_fit$lambda / (fixed_mu_fit$psi + fixed_mu_fit$mu)
    
    # Step 4: Check numerical stability
    if(is.na(a)){
      warning("Failed to estimate net rate 'a' from lineage-through-time data. Using MLE estimates only.")
    }else if(abs(b - 1) < 0.02){
      warning(sprintf(
        "Ratio b = %.4f is too close to 1 (difference: %.4f). Two-step method numerically unstable. Using MLE estimates only.",
        b, abs(b - 1)
      ))
    }else{
      # Step 5: Solve for lambda and psi using two-step method
      estimated_psi <- (a + mu - b * mu) / (b - 1)
      estimated_lambda <- a + mu + estimated_psi
      
      # Step 6: Validate two-step estimates
      constraints_satisfied <- TRUE
      constraint_violations <- c()
      
      # Check positivity
      if(estimated_psi < 0){
        constraints_satisfied <- FALSE
        constraint_violations <- c(constraint_violations, sprintf("psi < 0 (%.4f)", estimated_psi))
      }
      
      if(estimated_lambda < 0){
        constraints_satisfied <- FALSE
        constraint_violations <- c(constraint_violations, sprintf("lambda < 0 (%.4f)", estimated_lambda))
      }
      
      # Check R0 bounds
      if(estimated_psi >= 0 && estimated_lambda >= 0){
        removal_rate <- mu + estimated_psi
        estimated_r0 <- estimated_lambda / removal_rate
        
        # Check growth (R0 > 1)
        if(estimated_lambda <= removal_rate){
          constraints_satisfied <- FALSE
          constraint_violations <- c(constraint_violations, sprintf("R0 <= 1 (%.4f)", estimated_r0))
        }
        
        # Check R0 bounds
        if(estimated_r0 < r0_min){
          constraints_satisfied <- FALSE
          constraint_violations <- c(constraint_violations, sprintf("R0 < r0_min (%.4f < %.4f)", estimated_r0, r0_min))
        }
        
        if(estimated_r0 > r0_max){
          constraints_satisfied <- FALSE
          constraint_violations <- c(constraint_violations, sprintf("R0 > r0_max (%.4f > %.4f)", estimated_r0, r0_max))
        }
        
        # Check psi bounds
        if(estimated_psi > psi_max){
          constraints_satisfied <- FALSE
          constraint_violations <- c(constraint_violations, sprintf("psi > psi_max (%.4f > %.4f)", estimated_psi, psi_max))
        }
        
        # Check infectious period bounds
        inf_period <- 1 / removal_rate
        if(!is.null(infectious_period_min) && inf_period < infectious_period_min){
          constraints_satisfied <- FALSE
          constraint_violations <- c(constraint_violations, sprintf("inf.period < min (%.4f < %.4f)", inf_period, infectious_period_min))
        }
        
        if(!is.null(infectious_period_max) && inf_period > infectious_period_max){
          constraints_satisfied <- FALSE
          constraint_violations <- c(constraint_violations, sprintf("inf.period > max (%.4f > %.4f)", inf_period, infectious_period_max))
        }
      }
      if(constraints_satisfied){
        two_step_successful <- TRUE
      }
    }
    
  }, error = function(e){
    warning(sprintf("Estimation failed: %s", e$message))
    stop("Cannot proceed. Try relaxing constraints or checking tree structure.")
  })
  
  # Determine final estimates
  if(two_step_successful && !force_two_step){
    final_lambda <- estimated_lambda
    final_psi <- estimated_psi
    estimation_method <- "two_step_constrained"
  }else if(force_two_step && !is.na(estimated_lambda) && !is.na(estimated_psi)){
    # User forced two-step even if constraints violated
    final_lambda <- estimated_lambda
    final_psi <- estimated_psi
    estimation_method <- "two_step_forced"
    warning("Using two-step estimates despite constraint violations because force_two_step=TRUE.")
  }else{
    final_lambda <- fixed_mu_fit$lambda
    final_psi <- fixed_mu_fit$psi
    estimation_method <- "mle_constrained"
  }
  
  final_mu <- mu
  final_inf_period <- 1 / (final_mu + final_psi)
  
  # Return results
  result <- list(
    lambda = final_lambda,
    psi = final_psi,
    mu = final_mu,
    r0 = final_lambda / (final_mu + final_psi),
    infectious_period = final_inf_period,
    method = estimation_method,
    a = a,
    b = b,
    n_tips = length(phy$tip.label),
    n_nodes = phy$Nnode,
    mle_fit = fixed_mu_fit,
    constraints = list(
      psi_max = psi_max,
      r0_min = r0_min,
      r0_max = r0_max,
      infectious_period_min = infectious_period_min,
      infectious_period_max = infectious_period_max,
      mu_fixed = TRUE
    )
  )
  return(result)
}

