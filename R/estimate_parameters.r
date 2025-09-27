estimate_parameters <- function(phylobject, par_df, n_state, n_node, post_order_edges, topology) {
    # starting values
    temp <- par_df
    lambda <- rep(0.01, n_state)
    mu <- rep(0.005, n_state)
    Q <- rep(0.01, n_state^2)
    initial_state <- c(lambda, mu, Q)

    # lower bounds
    lb <- rep(0, length(initial_state))

    compute_phylobject <- function(x) {
        temp$lambda <- x[1:n_state]
        temp$mu <- x[1:n_state+n_state]
        Q <- matrix(x[1:n_state^2+2*n_state], nrow=n_state)
        logliks <- get_backwards_likelihoods_list(phylobject, temp, Q, n_state, n_nodes, post_order_edges, topology)
        root_loglik <- logliks[[phylobject[["Nnode"]]+2]]
        return(sum(root_loglik*par_df$prior))
    }

    MLE <- optim(initial_state, compute_phylobject, method="L-BFGS-B", lower=lb, control=list(fnscale=-1))
    return(list(lambda=MLE$par[1:n_state], 
                mu=MLE$par[1:n_state+n_state], 
                Q=matrix(MLE$par[1:n_state^2+2*n_state], nrow=n_state)))
}
