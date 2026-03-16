# uncomment this to run the full analysis, otherwise use rds file
# 
# source("Code/Simulation/R/simulation.R")
# 
# library(microbenchmark)
# library(ggplot2)
# library(dplyr)
# 
# # Function to generate parameters based on k
# generate_pars <- function(k) {
#   c(
#     rep(c(3, 1.5), length.out = k),      # speciation rates
#     rep(c(0.1, 0.05), length.out = k),   # extinction rates
#     rep(0.1, k),                          # another set of rates
#     rep(0.3, k**2-k)                         # transition rates
#   )
# }
# 
# # Function to generate q_matrix
# qij_matrix <- function(k) {
#   mat <- matrix(0.3, nrow = k, ncol = k)  
#   diag(mat) <- NA  
#   return(mat)
# }
# 
# # Initialize results storage
# timing_results <- data.frame(
#   k = integer(),
#   median_time = numeric(),
#   mean_time = numeric(),
#   min_time = numeric(),
#   max_time = numeric()
# )
# 
# # Set seed for reproducibility
# set.seed(1)
# 
# # Loop through different values of k
# for (k in 2:30) {
#   cat("Running benchmark for k =", k, "\n")
#   
#   # Generate parameters
#   pars <- generate_pars(k)
#   q_matrix <- qij_matrix(k)
#   
#   # Simulate tree
#   phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 10000, 10000)
#   
#   # Ensure the tree has at least 1000 tips
#   while(length(phy$tip.label) < 1000) {
#     phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 10000, 10000)
#   }
#   
#   phy <- prune(phy)
#   h <- history.from.sim.discrete(phy, 1:k)
#   
#   # Drop tips near TMRCA
#   node_depths <- node.depth.edgelength(phy)
#   tmrca <- max(node_depths)
#   tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
#   
#   if(length(tips_to_drop) > 0) {
#     phy <- drop.tip(phy, tips_to_drop)
#     phy$tip.state <- phy$tip.state[setdiff(names(phy$tip.state), tips_to_drop)]
#   }
#   
#   cat("Tree has", length(phy$tip.label), "tips and", phy$Nnode, "nodes\n")
#   
#   # Benchmark
#   timing <- microbenchmark(
#     saasi(pars, phy, q_matrix),
#     times = 1  # Adjust as needed
#   )
#   
#   # Store results (convert to seconds)
#   timing_results <- rbind(timing_results, data.frame(
#     k = k,
#     median_time = median(timing$time) / 1e9,  # Convert nanoseconds to seconds
#     mean_time = mean(timing$time) / 1e9,
#     min_time = min(timing$time) / 1e9,
#     max_time = max(timing$time) / 1e9
#   ))
#   
#   cat("Median time:", median(timing$time) / 1e9, "seconds\n\n")
# }
# 
# # Display results
# print(timing_results)
# 
# # Create visualization
# p1 <- ggplot(timing_results, aes(x = k, y = median_time)) +
#   geom_point(color = "#FFB6C1", size = 3) +
#   labs(
#     title = "SAASI Running Time vs Number of States",
#     x = "Number of States (k)",
#     y = "Time (seconds)"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 12, family = "serif"),
#     plot.title = element_text(size = 14, face = "bold")
#   )
# 
# saveRDS(timing_results,"running_time_wstates.rds")
# 
# # Save results
# write.csv(timing_results, "saasi_timing_results.csv", row.names = FALSE)

# RDS
d <- readRDS("Code/Results/suppfig12.rds")

p <- ggplot(d, aes(x = k, y = median_time)) +
  geom_point(color = "#FFB6C1", size = 3) +
  labs(
    title = "SAASI Running Time vs Number of States",
    x = "Number of States (k)",
    y = "Time (seconds)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12, family = "serif"),
    plot.title = element_text(size = 14, face = "bold")
  )