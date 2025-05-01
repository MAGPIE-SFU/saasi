<!-- README.md is generated from README.Rmd. Please edit that file -->

# saasi

<!-- badges: start -->
<!-- badges: end -->

Ancestral state reconstruction method that accounts for variation in
sampling rate among locations.

## Installation

You can install the development version of saasi from
[GitHub](https://github.com/) with:

     # install.packages("remotes")
     remotes::install_github("ivansg44/saasi")

## saasi package

This is a demo showing how to use the saasi package.

To run saasi, we need (1). a phylogenetic tree (class phylo) (2).
speciation, extinction and sampling rates (class data.frame) and (3).
transition rate matrix (class matrix). The output will be a data frame
that showing the probabilities of each node in state *i*.

## simulation

To simulate a birth-death-sampling tree, we need to specify speciation,
extinction, sampling rates and transition rates.

    pars <- data.frame(state=c(1,2),freq=c(0.5,0.5),lambda=c(3,3),mu=c(0.05,0.05),psi=c(.1,1))

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

You might notice that the tree also include all the tips at the present
day. This is because the simulation stop at the max time. To drop the
tips at the present day, we need to


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

Now we can use the simulated tree to do ancestral state inference.
First, we test what does the ace do

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

We see that ace would infer that most of the internal nodes are in State
2 instead of State 1.

Now let’s try saasi

    result <- saasi(phy,pars,q_matrix)
    node_result <-  result[-(1:22), ]
    node_result
    #>                  1            2
    #> nd1   9.823747e-01 1.762526e-02
    #> nd2   9.994883e-01 5.116901e-04
    #> nd12  9.987111e-01 1.288930e-03
    #> nd21  9.787996e-01 2.120044e-02
    #> nd5   9.898657e-01 1.013425e-02
    #> nd6   9.952304e-01 4.769623e-03
    #> nd10  9.327189e-01 6.728112e-02
    #> nd57  9.733160e-01 2.668405e-02
    #> nd59  9.999765e-01 2.352009e-05
    #> nd60  8.618376e-03 9.913816e-01
    #> nd74  3.123013e-06 9.999969e-01
    #> nd137 1.984286e-04 9.998016e-01
    #> nd3   8.695730e-01 1.304270e-01
    #> nd14  9.903516e-01 9.648423e-03
    #> nd31  8.159241e-01 1.840759e-01
    #> nd27  2.219318e-03 9.977807e-01
    #> nd38  4.649720e-05 9.999535e-01
    #> nd29  5.227194e-03 9.947728e-01
    #> nd81  2.417752e-03 9.975822e-01
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

<img src="man/figures/README-saasi-1.png" width="100%" /> Now saasi
infers most of the internal nodes correctly.
