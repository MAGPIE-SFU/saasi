#' Calculate c1 parameter for birth-death-sampling model
#'
#' @param lambda Speciation rate
#' @param mu Extinction rate
#' @param psi Sampling rate
#' @return Value of c1
#' @keywords internal
#' @noRd
c1 <- function(lambda, mu, psi){
  abs(sqrt((lambda - mu - psi)^2 + 4 * lambda * psi))
}

#' Calculate c2 parameter for birth-death-sampling model
#'
#' @param lambda Speciation rate
#' @param mu Extinction rate
#' @param psi Sampling rate
#' @return Value of c2
#' @keywords internal
#' @noRd
c2 <- function(lambda, mu, psi){
  (-lambda + mu + psi) / c1(lambda, mu, psi)
}

#' Numerically stable log-sum-exp function
#'
#' @param a First log value
#' @param b Second log value
#' @param c Third log value
#' @return log(exp(a) + exp(b) + exp(c))
#' @keywords internal
#' @noRd
logsumexp <- function(a, b, c){
  m <- pmax(a, pmax(b, c))
  m + log(exp(a - m) + exp(b - m) + exp(c - m))
}

#' Calculate log of q function for birth-death-sampling model
#'
#' @param t Time values
#' @param lambda Speciation rate
#' @param mu Extinction rate
#' @param psi Sampling rate
#' @return Log of q function values
#' @keywords internal
#' @noRd
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

#' Calculate log-likelihood under the general birth-death-sampling model
#'
#' @param phy A phylo object
#' @param lambda Speciation rate
#' @param mu Extinction rate
#' @param psi Sampling rate
#' @return Log-likelihood value
#' @keywords internal
#' @noRd
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


#' Estimate net diversification rate from lineage-through-time data
#'
#' Estimates the net diversification rate (lambda - mu - psi) by fitting a linear
#' model with the branch lengths, computed using a trimmed portion of the internal node times,
#' as the independent variable and the log-transformed lineage-through-time counts as the dependent variable.
#'
#' @param phy A phylo object with branch lengths
#' @param trim Numeric vector of length 2 specifying quantiles for trimming.
#'   Default is c(0.10, 0.50), using the middle 40% of node times.
#' @return Estimated net diversification rate (a = lambda - mu - psi)
#' @keywords internal
#' @noRd
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
  
  tmin <- stats::quantile(bt, trim[1])
  tmax <- stats::quantile(bt, trim[2])
  keep <- bt >= tmin & bt <= tmax
  fit <- stats::lm(log(L[keep]) ~ bt[keep])
  unname(stats::coef(fit)[2])
}

#' Generate starting points for optimization
#'
#' Creates valid starting points for lambda and psi that satisfy all constraints
#' including R0 bounds, removal rate bounds, and growth requirements.
#' 
#' Attempts to generate n_starts values for lambda and psi.
#' However if more than 100*n_starts attempts are made, then fewer than n_starts values are returned
#' and a warning message is displayed to the user.
#'
#' @param n_starts Number of starting points to generate
#' @param mu Extinction rate (fixed)
#' @param psi_max Maximum sampling rate per unit time
#' @param psi_min Minimum sampling rate
#' @param r0_min Minimum basic reproductive number (lambda/(mu+psi))
#' @param r0_max Maximum basic reproductive number
#' @param removal_rate_min Minimum removal rate (1/(mu+psi))
#' @param removal_rate_max Maximum removal rate (1/(mu+psi))
#' @return List of starting point vectors, each with a valid lambda and psi
#' @keywords internal
#' @noRd
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
    lambda_ini <- stats::runif(1, lambda_min, lambda_max)
    psi_ini <- stats::runif(1, psi_min, psi_max)
    
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
    warning("Fewer than n_starts starting points could be generated. Consider relaxing parameter constraints.")
  }
  return(starts)
}

#' Fit birth-death-sampling model with fixed extinction rate
#'
#' Estimates lambda (speication rate) and psi (sampling rate) using maximum likelihood
#' with a fixed mu (extinction rate), subject to some epidemiological constraints.
#'
#' @param phy A phylo object with branch lengths
#' @param mu Fixed extinction rate
#' @param n_starts Number of starting points
#' @param psi_max Maximum sampling rate 
#' @param r0_min Minimum R0 
#' @param r0_max Maximum R0
#' @param infectious_period_min Minimum infectious period (1/(mu+psi))
#' @param infectious_period_max Maximum infectious period (1/(mu+psi))
#' @return List containing estimated parameters and diagnostics
#' @keywords internal
#' @noRd
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
    stop("Failed to generate any valid starting points. Relax the parameter constraints.")
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
      stats::optim(
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

#' Estimate rate parameters of a birth-death-sampling process
#'
#' `estimate_bds_parameters()` estimates the speciation rate \eqn{\lambda} and sampling rate \eqn{\psi} of a birth-death-sampling process given a phylogenetic tree and death rate \eqn{\mu}.
#' The function uses user-specified epidemiological constraints on the sampling rate, basic reproduction number, and the duration of the infectious period to improve estimability.
#'
#' Estimation is carried out by first estimating the reproduction number \eqn{\frac{\lambda}{\mu+\psi}} and the net diversification rate \eqn{\lambda-\mu-\psi},
#' then solving for the unknown \eqn{\lambda} and \eqn{\psi}. The net diversification rate is estimated by
#' performing a lineages-through-time (LTT) regression and the reproduction number is estimated by first computing maximum likelihood 
#' estimates (MLEs) of \eqn{\lambda} and \eqn{\psi} by numerically optimizing a likelihood and then compute the reproduction number.
#' The stability of the MLEs is assessed by checking if the estimate of \eqn{R_0} is within 0.02 of 1.
#' 
#' If the estimates of \eqn{\lambda} and \eqn{\psi} do not satisfy all of the constraints specified by the `psi_max`, `r0_min`, `r0_max`, `infectious_period_min`,
#' and `infectious_period_max` parameters, then the MLEs of \eqn{\lambda} and \eqn{\psi} are returned, unless `force_two_step=TRUE`
#' in which case the original estimates of \eqn{\lambda} and \eqn{\psi} are returned regardless of 
#' whether or not the constraints are satisfied.
#' If LTT regression fails then the MLEs of \eqn{\lambda} and \eqn{\psi} will be returned. 
#' 
#' The MLEs are obtained by running numerical optimization (see \link[stats::optim()]{optim}) `n_starts` times from randomly chosen starting points
#' and selecting the best maximizer of all attempts. Due to the user-imposed constraints, not all randomly generated starting
#' points are valid, so `100*n_starts` attempts are made to generate valid starting points. This may result in fewer than `n_starts` optimizations
#' being performed.
#'
#' If the two-step method produces estimates that violate constraints or is
#' numerically unstable, the function falls back to the MLE estimates from step 1.
#'
#' @param phy An object of class `phylo` with branch lengths. Must be rooted and binary.
#' @param mu The extinction rate. This must be non-negative.
#' @param trim A numeric vector of length 2 specifying the quantiles for trimming
#'   node times in the lineages-through-time regression (see Details). Default value is `c(0.10, 0.50)`.
#' @param n_starts The number of starting points for multi-start optimization. Default value is `100`.
#' @param force_two_step Logical. If `TRUE`, uses two-step estimates even if they
#'   violate constraints (see Details). Default value is `FALSE`.
#' @param psi_max Maximum sampling rate per unit time. Must be strictly positive. Default value is `7`.
#' @param r0_min Minimum basic reproductive number (R0 = lambda/(mu+psi)). Must be strictly positive. Default 1.
#' @param r0_max Maximum basic reproductive number. Must be strictly positive. Default 5.
#' @param infectious_period_min Minimum infectious period (1/(mu+psi)). Must be strictly positive. Default NULL.
#' @param infectious_period_max Maximum infectious period (1/(mu+psi)). Must be strictly positive. Default NULL.
#'
#' @return An object of class `bds_estimate` with the following attributes:
#' \itemize{
#'   \item `lambda`: Estimated speciation rate
#'   \item `psi`: Estimated sampling rate
#'   \item `mu`: Fixed death rate (input)
#'   \item `r0`: Estimated basic reproductive number (lambda/(mu+psi))
#'   \item `infectious_period`: Estimated infectious period
#'   \item `method`: Estimation method used (`"two_step_constrained"`, `"mle_constrained"`, or `"two_step_forced"`)
#'   \item `a`: Net diversification rate (lambda - mu - psi) from LTT
#'   \item `b`: Ratio lambda/(mu + psi) from MLE
#'   \item `n_tips`: Number of tips in the tree
#'   \item `n_nodes`: Number of internal nodes
#'   \item `mle_fit`: Full MLE fit object from step 1
#'   \item `constraints`: List of constraint parameters used
#' }
#' @seealso [estimate_transition_rates()]
#' @examples
#' data(ebola_tree)
#' 
#' # Use the following information about the 2013 Ebola outbreak to obtain estimates
#' # - Average removal rate of 5
#' # - Infectious periods range from 20 to 40 days
#' # - An upper bound on the sampling rate of 15
#' # - Plausible values for the basic reproduction number are between 1.5 and 3
#' mu <- 5
#' BDS_fit <- estimate(ebola_tree, mu = 5,
#'                                 psi_max = 15,
#'                                 infectious_period_min = 20/365,
#'                                 infectious_period_max = 40/365,
#'                                 r0_min = 1.5,
#'                                 r0_max = 3)
#'
#' @export
estimate_bds_parameters <- function(
    phy,
    mu,
    trim = c(0.10, 0.50),
    n_starts = 100,
    force_two_step = FALSE,
    psi_max = 7,
    r0_min = 1,
    r0_max = 5,
    infectious_period_min = NULL,
    infectious_period_max = NULL
){
  # Input validation
  if(missing(mu) || is.null(mu)){
    stop("Extinction rate mu must have a numeric value.")
  }
  
  if(!inherits(phy, "phylo")){
    stop("Input phy must be an object with class phylo.")
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
      warning("Failed to estimate the net rate from lineage-through-time data. Using maximum likelihood estimates only.")
    }else if(abs(b - 1) < 0.02){
      warning(sprintf(
        "Ratio lambda/(mu + psi) = %.4f is too close to 1 (difference: %.4f). This suggests that the two-step method is numerically unstable, so maximum likelihood estimates will be used.",
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

