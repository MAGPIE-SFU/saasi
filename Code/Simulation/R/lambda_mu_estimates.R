# An implementation from Stadler et al 2012 & Stadler 2010
library(stats4)
# inputs: node time (yi) and tip time (xi), maybe psi
# outputs: ML estimates for speciation and extinction rate

# Here we calculate the probability density of the tree conditioned on 
# 1 extand sampled individuals.

# See eq.9 , Theorem 3.11 for detail.

# Note that in our model we do not consider sampled individual with 
# sampled descendants, therefore k = 0.

#c1
c1 <- function(lambda,mu,psi){
  return(abs(sqrt((lambda-mu-psi)^2+4*lambda*psi)))
}
#c2
c2<- function(lambda,mu,psi){
  return(-(lambda-mu-psi)/(c1(lambda,mu,psi)))
}

p0 <- function(lambda,mu,psi,t){
  d1 <- lambda + mu + psi
  d2 <- 2 * lambda
  d3 <- exp(-c1(lambda,mu,psi)*t)*(1-c2(lambda,mu,psi)) - (1+c2(lambda,mu,psi))
  d4 <- exp(-c1(lambda,mu,psi)*t)*(1-c2(lambda,mu,psi)) + (1+c2(lambda,mu,psi))
  d5 <- (c1(lambda,mu,psi)*d3)/d4
  return((d1+d5)/d2)
}

q <- function(lambda,mu,psi,t){
  d1 <- 2 * (1-(c2(lambda,mu,psi))^2)
  d2 <- exp(-c1(lambda,mu,psi)*t)*(1-c2(lambda,mu,psi))^2
  d3 <- exp(c1(lambda,mu,psi)*t)*(1+c2(lambda,mu,psi))^2
  return(d1+d2+d3)
}


# probability density 

probability_density <- function(lambda,mu,psi,m,X,Y){
  # Add checks for invalid inputs
  sortX <- sort(X)
  sortY <- sort(Y)
  d1 <- m*log(lambda)
  d2 <- sum(-log(q(lambda,mu,psi,sortX)))
  d3 <- sum(log(psi)+log(q(lambda,mu,psi,sortY)))
  result <- d1+d2+d3
  return(result)
}



##### 
extract_tree_times <- function(tree) {
  # Get total tree height (max path length from root to tip)
  tree_height <- max(node.depth.edgelength(tree))
  
  # Get all node depths (distance from root)
  all_nodes_depth <- node.depth.edgelength(tree)
  
  # Convert depths to times (time from present)
  all_times <- tree_height - all_nodes_depth
  
  # Separate internal nodes and tips
  n_tips <- length(tree$tip.label)
  
  # Extract tip times
  tip_times <- all_times[1:n_tips]
  names(tip_times) <- tree$tip.label
  
  # Extract internal node times
  internal_times <- all_times[(n_tips + 1):length(all_times)]
  names(internal_times) <- (n_tips + 1):length(all_times)
  
  # Create data frames
  internal_df <- data.frame(
    node = names(internal_times),
    time = as.numeric(internal_times),
    row.names = NULL
  )
  
  tip_df <- data.frame(
    tip = names(tip_times),
    time = as.numeric(tip_times),
    row.names = NULL
  )
  
  return(list(
    internal_nodes = internal_df,
    tips = tip_df
  ))
}


branching_time <- extract_tree_times(phy)
X <- as.vector(branching_time$internal_nodes)
X <- X$time
Y <- as.vector(branching_time$tips)
Y <- Y$time

negative_log_likelihood <- function(params, m, X, Y) {
  lambda <- params[1]
  mu <- params[2]
  psi <- params[3]
  
  # Ensure parameters are positive
  if (lambda <= 0 || mu <= 0 || psi <= 0) return(Inf)
  
  # Compute the probability density
  likelihood <- probability_density(lambda, mu, psi, m, X, Y)
  likelihood
  # Return negative log-likelihood
  return(-likelihood)
}


#Find MLE with multiple starting points


# Define the negative log likelihood function
nLL <- function(lambda, mu,psi) {
  # Return negative log likelihood
  # We use negative because mle() minimizes
  -(probability_density(lambda, mu,psi, m=phy$Nnode, X=X, Y=Y))
}

# Find MLE
fit <- mle(nLL, 
           start = list(lambda = 1, mu = 1,psi=0.2),
           method = "L-BFGS-B",
           lower = c(0.0001, 0.0001,0.0001),
           upper = c(5, 5,5))

# View results
coef(fit)  # Get parameter estimates
