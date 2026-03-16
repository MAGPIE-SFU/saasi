library(irr)

######################
# This function returns 4 accuracy measurements 

accuracy_helper <- function(asr,true_states,tip_states){
  
  reconstructed_states <- apply(asr, 1, which.max)
  comparison_data <- data.frame(
    true = true_states,
    reconstructed = c(tip_states, reconstructed_states)
  )
  comparison_data <- tail(comparison_data,length(tip_states)+1)
  
  # Overall accuracy
  accuracy <- mean(comparison_data$true == comparison_data$reconstructed)
  
  prob_accuracy <- mean(sapply(1:nrow(asr), function(i) {
    asr[i, true_states[length(tip_states) + i]]
  }))
  return(list(
    acc=accuracy,
    prob_acc=prob_accuracy
  ))
}