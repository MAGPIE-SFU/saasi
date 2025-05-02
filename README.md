<!-- README.md is generated from README.Rmd. Please edit that file -->

# Sampling Aware Ancestral State Inference (saasi)

<!-- badges: start -->
<!-- badges: end -->

Saasi is an ancestral state reconstruction method that accounts for
variation in sampling rates among locations or traits.

## Installation

You can install the development version of saasi from
[GitHub](https://github.com/) with:

     # install.packages("remotes")
     remotes::install_github("MAGPIE-SFU/saasi")

## The saasi package

This is a demo showing how to use the saasi package.

To run saasi, we need (1) a phylogenetic tree (class `phylo`), (2)
speciation, extinction and sampling rates (class `data.frame`), and (3)
a transition rate matrix (class `matrix`). The output will be a data
frame that containing the probabilities of each state for each internal
node of the phylogenetic tree.

## A simulation

We will first simulate a tree with known rates, and known internal node
states, to illustrate saasi’s ancestral state inferences.

We simulate a birth-death-sampling tree, for which we need to specify
speciation, extinction, sampling rates and transition rates.

    pars <- data.frame(state=c(1,2),prior=c(0.5,0.5),lambda=c(3,3),mu=c(0.05,0.05),psi=c(.1,1))

    qij_matrix <- function(k) {
      mat <- matrix(0.15, nrow = k, ncol = k)  
      diag(mat) <- NA  
      return(mat)
    }
    q_matrix = qij_matrix(2)

    set.seed(1)

    phy <- sim_bds_tree(pars, q_matrix, x0=1, max_taxa = 300, max_t = 300,
                 include_extinct = FALSE)

    plot(phy)

<img src="man/figures/README-parameters-1.png" width="100%" />

## Modifying tree

You might notice that the tree includes all the tips at the present day.
This is because the simulation stopped at the maximum allowed time. In
pathogen phylogenetics and phylogeography applications, we typically do
not have heterochronous sequences (from the present). In this
simulation, we therefore drop the tips at the present day.


    k=2
    phy <- prune(phy)
    h <- history.from.sim.discrete(phy, 1:k)

    true_phy_info <- as_tibble(phy)
    phy_data <- c(factor(h$tip.state),factor(h$node.state))
    true_phy_info$State <- phy_data
    true_phy <- ggtree(phy)
    true_phy <- true_phy  %<+% true_phy_info + geom_point(aes(color=State),size=2) +
      ggtitle("True Phylogeny") +
      theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
    true_phy

<img src="man/figures/README-modified-1.png" width="100%" />

    node_depths <- node.depth.edgelength(phy)
    tmrca <- max(node_depths)
    tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
    new_phy <- drop.tip(phy, tips_to_drop)
    phy_our <- new_phy
    true_phy_info_new <- as_tibble(new_phy)
    phy_data <- c(factor(h$tip.state),factor(h$node.state))
    true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
    new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
    phy <- new_phy
    true_phy_new <- ggtree(phy)
    true_phy_new <- true_phy_new  %<+% true_phy_info_new + geom_point(aes(color=State),size=2) +
      ggtitle("True Phylogeny - without present day tips") +
      theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
    true_phy_new

<img src="man/figures/README-modified-2.png" width="100%" />

## Ancestral state inference

Now we can use the simulated tree to do ancestral state inference. The
function `ace` in the `ape` package does ancestral character (here,
state) estimation without considering sampling rates, and it is a
natural comparison for saasi since it is widely used in large-scale
phylogeographic reconstructions. Comparing `saasi`’s reconstructions to
`ace`’s illustrates the impact of adjusting for sampling differences.

    library(ape)
    ace_phy <- phy
    ace_phy$node.label <- NULL
    # Note: Do not have this problem if use earlier version `ape`
    # Error in names(obj$ace) <- phy$node.label : 
    # attempt to set an attribute on NULL

    ace_phy$tip.state <- ace_phy$tip.state[setdiff(names(ace_phy$tip.state), tips_to_drop)]
    asr_ace<-ace(ace_phy$tip.state, ace_phy,type = "discrete", model="ER")

    ace_node_lik <- as.data.frame(asr_ace$lik.anc)
    ace_node_lik$node <- 1:new_phy$Nnode + Ntip(new_phy)

    ace_pie <- nodepie(ace_node_lik,cols=1:k)

    p1 <- ggtree(ace_phy)
    p1 <- p1 %<+% true_phy_info_new + geom_tippoint(aes(color=State),size=2)+
      ggtitle("ace") +
      theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
    p1 <- inset(p1, ace_pie,width = 0.07,height = 0.07,hjust=0.005)
    p1

<img src="man/figures/README-ace-1.png" width="100%" />

We see that `ace` would infer that most of the internal nodes are in
State 2 instead of State 1.

Now let’s try `saasi`.

    result <- saasi(phy,pars,q_matrix)
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    node_result <-  result[-(1:22), ]
    node_result
    #>                  1            2
    #> nd1   9.823748e-01 1.762524e-02
    #> nd2   9.994883e-01 5.116896e-04
    #> nd12  9.987111e-01 1.288924e-03
    #> nd21  9.787996e-01 2.120035e-02
    #> nd5   9.898658e-01 1.013425e-02
    #> nd6   9.952304e-01 4.769618e-03
    #> nd10  9.327190e-01 6.728104e-02
    #> nd57  9.733160e-01 2.668404e-02
    #> nd59  9.999765e-01 2.352009e-05
    #> nd60  8.618374e-03 9.913816e-01
    #> nd74  3.123013e-06 9.999969e-01
    #> nd137 1.984286e-04 9.998016e-01
    #> nd3   8.695732e-01 1.304268e-01
    #> nd14  9.903516e-01 9.648396e-03
    #> nd31  8.159245e-01 1.840755e-01
    #> nd27  2.219324e-03 9.977807e-01
    #> nd38  4.649731e-05 9.999535e-01
    #> nd29  5.227208e-03 9.947728e-01
    #> nd81  2.417753e-03 9.975822e-01
    #> nd112 1.710792e-05 9.999829e-01
    #> nd193 4.968901e-05 9.999503e-01
    node_result$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
    our_pie <- nodepie(node_result,cols=1:k)

    p2 <- ggtree(ace_phy)
    p2 <- p2 %<+% true_phy_info_new + geom_tippoint(aes(color=State),size=2)+
      ggtitle("SAASI") +
      theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
    p2 <- inset(p2, our_pie,width = 0.07,height = 0.07,hjust=0.005)
    p2

<img src="man/figures/README-saasi-1.png" width="100%" />

Due to accounting for the sampling differences `saasi` infers most of
the internal nodes correctly.

## Parameter estimations - Speciation and Extinction

Suppose we do not know the speciation and extinction rates for each
state, but we have some knowledge about the sampling rates (e.g. per
year/month, this should algin with your tree time). We estimate the
speciation and extinction rates using the method described in .


    estimates <- mle_lm(phy,lambda = 2, mu = 0.1, psi = 1,lower = c(0.001,0.001), upper = c(5,5))

    estimates[1]
    #> [1] 2.086275
    estimates[2]
    #> [1] 0.001

This MLE approach should be robust to different initial guesses, for
example:


    estimates <- mle_lm(phy,lambda = 100, mu = 100, psi = 1,lower = c(0.001,0.001), upper = c(500,500))

    estimates[1]
    #> [1] 2.086275
    estimates[2]
    #> [1] 0.001

However, sometimes the algorithm might cause an error: Error: L-BFGS-B
needs finite values of \`fn’ of Complex Objective Function. This is due
to a large branch length (the value exceed the .Machine$double.xmax), we
need to set the upper bound smaller.

TODO: for future versions, will find a more consistent way of finding
the estimates.

## Parameter estimations - Transition

If the transition rates are also unknown, one easy way of estimating
transition rates is using `ace`:


    # a function that convert the ace estimates to a matrix (one of the inputs in saasi)
    replace_matrix_with_vector <- function(matrix, vector) {
      for (i in 1:nrow(matrix)) {
        for (j in 1:ncol(matrix)) {
          matrix[i,j] <- vector[matrix[i,j]]
        }
      }
      return(matrix)
    }

    # TODO: will add this to the utility function

    q_matrix_est <- replace_matrix_with_vector(asr_ace$index.matrix,asr_ace$rates)
    q_matrix_est
    #>           [,1]      [,2]
    #> [1,]        NA 0.4950541
    #> [2,] 0.4950541        NA

Now we rerun `saasi` with estimated parameters.

    pars_est <- data.frame(state=c(1,2),
                       prior=c(0.5,0.5),
                       lambda=rep(estimates[1],2),
                       mu=rep(estimates[2],2),
                       psi=c(.1,1))

    result <- saasi(phy,pars_est,q_matrix_est)
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1
    #> [1] 0.1
    #> [1] 1
    #> [1] 1

    # now draw the plot

    node_result <-  result[-(1:22), ]
    node_result
    #>                  1           2
    #> nd1   0.9160470232 0.083952977
    #> nd2   0.9742168326 0.025783167
    #> nd12  0.9481213625 0.051878637
    #> nd21  0.8426787951 0.157321205
    #> nd5   0.8937366628 0.106263337
    #> nd6   0.9080256571 0.091974343
    #> nd10  0.7577604048 0.242239595
    #> nd57  0.8924616573 0.107538343
    #> nd59  0.9989037047 0.001096295
    #> nd60  0.0579786364 0.942021364
    #> nd74  0.0001308276 0.999869172
    #> nd137 0.0049735937 0.995026406
    #> nd3   0.6781574304 0.321842570
    #> nd14  0.8850996242 0.114900376
    #> nd31  0.6209197684 0.379080232
    #> nd27  0.0509743572 0.949025643
    #> nd38  0.0012343549 0.998765645
    #> nd29  0.0723530922 0.927646908
    #> nd81  0.0343703316 0.965629668
    #> nd112 0.0006489199 0.999351080
    #> nd193 0.0014922413 0.998507759
    node_result$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
    our_pie <- nodepie(node_result,cols=1:k)

    p3 <- ggtree(ace_phy)
    p3 <- p3 %<+% true_phy_info_new + geom_tippoint(aes(color=State),size=2)+
      ggtitle("SAASI - using estimated parameters") +
      theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
    p3 <- inset(p3, our_pie,width = 0.07,height = 0.07,hjust=0.005)
    p3

<img src="man/figures/README-rerun saasi with estimated parameters-1.png" width="100%" />

The result is slightly different than using the true parameters.
