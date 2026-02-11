# set seed
set.seed(123)

# specifying bds parameters
demo_pars <- data.frame(state=c("1","2","3"), 
                   rootprior=c(1/3,1/3,1/3), 
                   lambda=c(1.5,1.5,1.5), 
                   mu=c(0.3,0.3,0.3), 
                   psi=c(0.1,0.5,0.5))

# specifying Q matrix
demo_Q <- matrix(0.3, 3, 3)
diag(demo_Q) <- -0.6

# adding names
rownames(demo_Q) <- c("1","2","3")
colnames(demo_Q) <- c("1","2","3")

# generating trees
demo_tree_prepared <- sim_bds_tree(demo_pars, demo_Q, x0=1, max_taxa=80, 
                    max_t=10, include_extinct=FALSE, min_tip=10)

# save tree and data

demo_tree <- demo_tree_prepared
demo_tree$tip.state <- NULL
demo_metadata <- data.frame(
  node = demo_tree_prepared$tip.label,
  state = demo_tree_prepared$tip.state
)

save(demo_tree_prepared, file = "data/demo_tree_prepared.rda")
save(demo_Q, file = "data/demo_Q.rda")
save(demo_pars, file = "data/demo_pars.rda")
save(demo_tree, file = "data/demo_tree.rda")
save(demo_metadata, file = "data/demo_metadata.rda")

