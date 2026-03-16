replace_matrix_with_vector <- function(matrix, vector) {
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- vector[matrix[i,j]]
    }
  }
  return(matrix)
}


tree_simulation_w_error <- function(pars,q_matrix){
  # first generate the tree
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 1000, 1000)
  # Ensure the tree has at least 30 tips
  while(length(phy$tip.label) < 20) {
    phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 1000, 1000)
  }
  phy <- prune(phy)
  print('make the tree')
  h2 <- history.from.sim.discrete(phy, 1:2)
  
  node_depths <- node.depth.edgelength(phy)
  tmrca <- max(node_depths)
  tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.001]
  new_phy <- drop.tip(phy, tips_to_drop)
  
  plot(new_phy)
  phy_our <- new_phy
  
  true_phy_info_new <- as_tibble(new_phy)
  phy_data <- c(factor(h2$tip.state),factor(h2$node.state))
  true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
  new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
  
  phy <- new_phy
  
  
  tip_states <- h2$tip.state
  tip_states <- tip_states[!names(tip_states) %in% tips_to_drop]
  #####################
  # perform ace
  #####################
  phy <- as.phylo(phy)
  phy$node.label <- NULL
  # Extract tip states
  #tip_states <- phy$tip.state
  # Perform ASR
  asr_ace <- ace(tip_states, phy, type="discrete", model="ER")
  asr_rates <- asr_ace$rates
  
  qij_matrix <- replace_matrix_with_vector(asr_ace$index.matrix,asr_ace$rates)
  
  print('finish ace')

  #####################
  # perform saasi
  #####################
  asr_our <- saasi(pars,phy,q_matrix)
  #print(asr_our)
  print('finish our correct')
  qij_matrix <- replace_matrix_with_vector(asr_ace$index.matrix,asr_ace$rates)
  
  asr_our_est <- saasi(pars,phy,qij_matrix)
  #print(asr_our_est)
  
  qij_matrix2 <- replace_matrix_with_vector(asr_ace$index.matrix,asr_ace$rates)
  qij_matrix2[,1] = qij_matrix2[,1]*2
  qij_matrix2[1,] = qij_matrix2[1,]/2

  asr_our_est_2 <- saasi(pars,phy,qij_matrix2)
  
  
  qij_matrix3 <- replace_matrix_with_vector(asr_ace$index.matrix,asr_ace$rates)
  qij_matrix3[,1] = qij_matrix3[,1]*3
  qij_matrix3[1,] = qij_matrix3[1,]/3

  asr_our_est_3 <- saasi(pars,phy,qij_matrix3)
  
  qij_matrix4 <- replace_matrix_with_vector(asr_ace$index.matrix,asr_ace$rates)
  qij_matrix4[,1] = qij_matrix4[,1]*5
  qij_matrix4[1,] = qij_matrix4[1,]/5

  asr_our_est_4 <- saasi(pars,phy,qij_matrix4)
  
  qij_matrix5 <- replace_matrix_with_vector(asr_ace$index.matrix,asr_ace$rates)
  qij_matrix5[,1] = qij_matrix5[,1]*10
  qij_matrix5[1,] = qij_matrix5[1,]/10

  asr_our_est_5 <- saasi(pars,phy,qij_matrix5)

  print('finish our estimate')
  
  #####################
  # perform accuracy check
  #####################  
  
  acc_asr <- accuracy_helper(asr_ace$lik.anc,true_phy_info_new$State,tip_states)
  #print("1")
  #acc_sim <- accuracy_helper(head(simmap_lik,n=length(h2$node.state)- length(tips_to_drop)),true_phy_info_new$State,tip_states)
  #print("2")
  acc_our <- accuracy_helper(asr_our[,1:2],true_phy_info_new$State,tip_states)
  acc_our_est <- accuracy_helper(asr_our_est[,1:2],true_phy_info_new$State,tip_states)
  acc_our_est2 <- accuracy_helper(asr_our_est_2[,1:2],true_phy_info_new$State,tip_states)
  acc_our_est3 <- accuracy_helper(asr_our_est_3[,1:2],true_phy_info_new$State,tip_states)
  acc_our_est4 <- accuracy_helper(asr_our_est_4[,1:2],true_phy_info_new$State,tip_states)
  acc_our_est5 <- accuracy_helper(asr_our_est_5[,1:2],true_phy_info_new$State,tip_states)
  
  #return the accuracies
  print("?")
  return(list(
    asr=acc_asr,
    our=acc_our,
    our_est=acc_our_est,
    our_est2=acc_our_est2,
    our_est3=acc_our_est3,
    our_est4=acc_our_est4,
    our_est5=acc_our_est5
    
  ))
}

pars <- c(1,  1,      
          .045, .045,
          .01,  .5, 
          .05,
          .05)  
qij_matrix <- function(k,val) {
  mat <- matrix(val, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2,0.05)
# generate 10 trees
results1 <- list()
for (i in 1:100){
  print(i)
  set.seed(i)
  results1[[i]] <- tree_simulation_w_error(pars,q_matrix)
}


####################
# accuracy plots
####################

# ace accuracy
out_acc_ace <- c()
out_prob_ace <- c()
ace_rate <- c()
for(i in 1:length(results1)){
  out_acc_ace=c(out_acc_ace,results1[[i]]$asr$acc)
  out_prob_ace=c(out_prob_ace, results1[[i]]$asr$prob_acc)
  ace_rate=c(ace_rate,results1[[i]]$asr_rates)
}

# our accuracy
out_acc_our <- c()
out_prob_our <- c()
for(i in 1:length(results1)){
  out_acc_our=c(out_acc_our,results1[[i]]$our$acc)
  out_prob_our=c(out_prob_our, results1[[i]]$our$prob_acc)
}


# our accuracy est
out_acc_our_est <- c()
out_prob_our_est <- c()
for(i in 1:length(results1)){
  out_acc_our_est=c(out_acc_our_est,results1[[i]]$our_est$acc)
  out_prob_our_est=c(out_prob_our_est, results1[[i]]$our_est$prob_acc)
}

out_acc_our_est2 <- c()
out_prob_our_est2 <- c()
for(i in 1:length(results1)){
  out_acc_our_est2=c(out_acc_our_est2,results1[[i]]$our_est2$acc)
  out_prob_our_est2=c(out_prob_our_est2, results1[[i]]$our_est2$prob_acc)
}

out_acc_our_est3 <- c()
out_prob_our_est3 <- c()
for(i in 1:length(results1)){
  out_acc_our_est3=c(out_acc_our_est3,results1[[i]]$our_est3$acc)
  out_prob_our_est3=c(out_prob_our_est3, results1[[i]]$our_est3$prob_acc)
}

out_acc_our_est4 <- c()
out_prob_our_est4 <- c()
for(i in 1:length(results1)){
  out_acc_our_est4=c(out_acc_our_est4,results1[[i]]$our_est4$acc)
  out_prob_our_est4=c(out_prob_our_est4, results1[[i]]$our_est4$prob_acc)
}


out_acc_our_est5 <- c()
out_prob_our_est5 <- c()
for(i in 1:length(results1)){
  out_acc_our_est5=c(out_acc_our_est5,results1[[i]]$our_est5$acc)
  out_prob_our_est5=c(out_prob_our_est5, results1[[i]]$our_est5$prob_acc)
}

plot_data_acc <- data.frame(
  ace = out_acc_ace,
  SAASI = out_acc_our,
  SAASI_est = out_acc_our_est,
  SAASI_est2 = out_acc_our_est2,
  SAASI_est3 = out_acc_our_est3,
  SAASI_est4 = out_acc_our_est4,
  SAASI_est5 = out_acc_our_est5
)

acc_data <- pivot_longer(plot_data_acc,
                         cols = c(ace,SAASI,SAASI_est,SAASI_est2,SAASI_est3,SAASI_est4,SAASI_est5),
                         names_to = "Metric",
                         values_to = "Accuracy")
acc_data$Metric <- factor(acc_data$Metric, levels = c("ace","SAASI","SAASI_est","SAASI_est2","SAASI_est3","SAASI_est4","SAASI_est5"))

p <- ggplot(acc_data, aes(x = Metric, y = Accuracy)) +
  geom_violin(width = 1.1, fill = "lightblue", color = "black",alpha=0.3) +  # Remove fill, keep border
  #geom_jitter(height = 0, width = 0.1) +
  #geom_boxplot(width = 0.2, fill = "white", color = "black") +  # Keep boxplot visible
  labs(x = "", y = "Accuracy", fill = "Methods") +
  ggtitle("") +
  scale_x_discrete(labels = c(
    "ace" = "ace",
    "SAASI" = "SAASI - true",
    "SAASI_est" = "SAASI - est qij",
    "SAASI_est2" = "SAASI - adj by 2",
    "SAASI_est3" = "SAASI - adj by 3",
    "SAASI_est4" = "SAASI - adj by 5",
    "SAASI_est5" = "SAASI - adj by 10")) +
  theme_minimal() +
  ylim(0.0, 1.0) +
  theme(text = element_text(size = 15, family = "serif"),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

p

#ggsave(file="figs3_updated.pdf", width = 210, height = 297, units = "mm")


