# uncomment this to run the full analysis, otherwise use rds file

# library(ape)
# library(ggplot2)
# library(dplyr)
# 
# source("Code/Simulation/R/lost_info.R")
# 
# # NOTE: THIS SHOULD BE YOUR FILE DIR ONCE YOU RUN THE SHELL SCRIPT. THIS 
# # CODE GENERATES SUPPLMENTARY FIGURES 7 & 8.
# # Define your base directories
# base_dirs <- c(
#   "final_10x_drop/"
# )
# scenario_names <- c(
#   "baseline"
# )
# all_results_withoutds <- read.csv("final_10x/accuracy_results.csv")
# 
# 
# 
# saasi_accuracy_downsample <- function(tree,true_tree, true_tree_before_downsample,
#                                       true_tree_before_downsample_csv,asr, true.states) {
#   diff <- get_diff_hist(tree, true_tree,true_tree_before_downsample
#                         ,true_tree_before_downsample_csv)
#   reconstructed_states <- apply(asr, 1, which.max)
#   comparison_data <- data.frame(
#     true = tail(true.states$State, length(tree$tip.label) - 1),
#     reconstructed = reconstructed_states
#   )
#   
#   accuracy <- mean(comparison_data$true == comparison_data$reconstructed)
#   
#   prob_accuracy <- mean(sapply(1:nrow(asr), function(i) {
#     asr[i, true.states$State[length(tree$tip.label) + i]]
#   }))
#   correct_count <- sum(comparison_data$true == comparison_data$reconstructed)
#   accuracy_downsample <- correct_count / (length(tree$tip.label) - 1 + diff)
#   
#   prob_sum <- sum(sapply(1:nrow(asr), function(i) {
#     asr[i, true.states$State[length(tree$tip.label) + i]]
#   }))
#   prob_accuracy_downsample <- prob_sum / (nrow(asr) + diff)
#   
#   return(list(
#     acc = accuracy,
#     prob_acc = prob_accuracy,
#     acc_down = accuracy_downsample,
#     prob_acc_down = prob_accuracy_downsample
#   ))
# }
# 
# # Initialize combined results dataframe
# all_results <- data.frame()
# 
# # Loop through each base directory
# for(base_idx in seq_along(base_dirs)) {
#   base_dir <- base_dirs[base_idx]
#   scenario <- scenario_names[base_idx]
#   cat("\n=== Processing scenario:", scenario, "===\n")
#   
#   # Find all simulation directories
#   sim_dirs <- list.dirs(base_dir, recursive = FALSE)
#   sim_dirs <- sim_dirs[grep("sim_\\d+", sim_dirs)]
#   sim_nums <- as.numeric(gsub(".*sim_(\\d+)", "\\1", basename(sim_dirs)))
#   
#   # Added "20_0" for custom correction
#   treetime_corrections <- c("", "20_0")
#   treetime_labels <- c("1.0","Custom")
#   
#   # Process each simulation
#   for(i in seq_along(sim_dirs)) {
#     sim_dir <- sim_dirs[i]
#     sim_num <- sim_nums[i]
#     
#     cat("Processing simulation", sim_num, "\n")
#     
#     ####
#     
#     tryCatch({
#       # Read true states and tree
#       tree <- read.tree(file.path(sim_dir, paste0("tree_", sim_num, ".nwk")))
#       true_tree <- read.csv(file.path(sim_dir, paste0("true_tree_", sim_num, ".csv")))
#       true_tree_before_downsample <-  read.tree(file.path(sim_dir, paste0("tree_before_ds_", sim_num, ".nwk")))
#       true_tree_before_downsample_csv <-  read.csv(file.path(sim_dir, paste0("true_tree_before_ds_", sim_num, ".csv")))
#       
#       # Calculate tree size
#       tree_size <- tree$Nnode * 2 + 1
#       
#       # ACE
#       tryCatch({
#         ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
#         ace_acc <- saasi_accuracy_downsample(tree,
#                                              true_tree,
#                                              true_tree_before_downsample,
#                                              true_tree_before_downsample_csv,
#                                              ace_result[2:3], true_tree)
#         all_results <- rbind(all_results, data.frame(
#           scenario = scenario,
#           simulation = sim_num, 
#           method = "ACE", 
#           accuracy = ace_acc$acc, 
#           prob_accuracy = ace_acc$prob_acc,
#           accuracy_down = ace_acc$acc_down,
#           prob_accuracy_down = ace_acc$prob_acc_down,
#           tree_size = tree_size
#         ))
#       }, error = function(e) cat("ACE failed for sim", sim_num, "\n"))
#       
#       # SIMMAP
#       tryCatch({
#         simmap_result <- read.csv(file.path(sim_dir, paste0("simmap_", sim_num, ".csv")))
#         simmap_subset <- head(simmap_result[2:3], length(true_tree$State) - length(tree$tip.label))
#         simmap_acc <- saasi_accuracy_downsample(tree,
#                                                 true_tree,
#                                                 true_tree_before_downsample,
#                                                 true_tree_before_downsample_csv,
#                                                 simmap_subset, true_tree)
#         all_results <- rbind(all_results, data.frame(
#           scenario = scenario,
#           simulation = sim_num, 
#           method = "SIMMAP", 
#           accuracy = simmap_acc$acc, 
#           prob_accuracy = simmap_acc$prob_acc,
#           accuracy_down = simmap_acc$acc_down,
#           prob_accuracy_down = simmap_acc$prob_acc_down,
#           tree_size = tree_size
#         ))
#       }, error = function(e) cat("SIMMAP failed for sim", sim_num, "\n"))
#       
#       # SAASI
#       tryCatch({
#         saasi_result <- read.csv(file.path(sim_dir, paste0("saasi_", sim_num, ".csv")))
#         saasi_acc <- saasi_accuracy_downsample(tree,
#                                                true_tree,
#                                                true_tree_before_downsample,
#                                                true_tree_before_downsample_csv,
#                                                saasi_result[2:3], true_tree)
#         all_results <- rbind(all_results, data.frame(
#           scenario = scenario,
#           simulation = sim_num, 
#           method = "SAASI", 
#           accuracy = saasi_acc$acc, 
#           prob_accuracy = saasi_acc$prob_acc,
#           accuracy_down = saasi_acc$acc_down,
#           prob_accuracy_down = saasi_acc$prob_acc_down,
#           tree_size = tree_size
#         ))
#       }, error = function(e) cat("SAASI failed for sim", sim_num, "\n"))
#       
#       # SAASI with pars est
#       tryCatch({
#         saasi_pars_result <- read.csv(file.path(sim_dir, paste0("saasi_pars_est_", sim_num, ".csv")))
#         saasi_pars_acc <- saasi_accuracy_downsample(tree,
#                                                     true_tree,
#                                                     true_tree_before_downsample,
#                                                     true_tree_before_downsample_csv,
#                                                     saasi_pars_result[2:3], true_tree)
#         all_results <- rbind(all_results, data.frame(
#           scenario = scenario,
#           simulation = sim_num, 
#           method = "SAASI_pars", 
#           accuracy = saasi_pars_acc$acc, 
#           prob_accuracy = saasi_pars_acc$prob_acc,
#           accuracy_down = saasi_pars_acc$acc_down,
#           prob_accuracy_down = saasi_pars_acc$prob_acc_down,
#           tree_size = tree_size
#         ))
#       }, error = function(e) cat("SAASI_pars failed for sim", sim_num, "\n"))
#       
#       # PastML
#       tryCatch({
#         pastml_file <- file.path(sim_dir, "pastml_results", 
#                                  "marginal_probabilities.character_numeric_state.model_F81.tab")
#         if(file.exists(pastml_file)) {
#           pastml_result <- read.table(pastml_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#           pastml_subset <- pastml_result[grep("^(root|nd\\d+)$", pastml_result$node), ]
#           ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
#           pastml_matched <- pastml_subset[match(ace_result$X, pastml_subset$node), ]
#           pastml_acc <- saasi_accuracy_downsample(tree,
#                                                   true_tree,
#                                                   true_tree_before_downsample,
#                                                   true_tree_before_downsample_csv,
#                                                   pastml_matched[2:3], true_tree)
#           all_results <- rbind(all_results, data.frame(
#             scenario = scenario,
#             simulation = sim_num, 
#             method = "PastML", 
#             accuracy = pastml_acc$acc, 
#             prob_accuracy = pastml_acc$prob_acc,
#             accuracy_down = pastml_acc$acc_down,
#             prob_accuracy_down = pastml_acc$prob_acc_down,
#             tree_size = tree_size
#           ))
#         }
#       }, error = function(e) cat("PastML failed for sim", sim_num, "\n"))
#       
#       for(j in seq_along(treetime_corrections)) {
#         correction_suffix <- treetime_corrections[j]
#         correction_label <- treetime_labels[j]
#         
#         tryCatch({
#           dir_name <- if(correction_suffix == "") {
#             "treetime_results"
#           } else {
#             paste0("treetime_results_", correction_suffix)
#           }
#           treetime_file <- file.path(sim_dir, dir_name, "confidence.csv")
#           if(file.exists(treetime_file)) {
#             treetime_result <- read.csv(treetime_file)
#             treetime_subset <- treetime_result[grep("^(nd\\d+)$", treetime_result$X.name), ]
#             ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
#             treetime_matched <- treetime_subset[match(ace_result$X, treetime_subset$X.name), ]
#             treetime_acc <- saasi_accuracy_downsample(tree,
#                                                       true_tree,
#                                                       true_tree_before_downsample,
#                                                       true_tree_before_downsample_csv,
#                                                       treetime_matched[2:3], true_tree)
#             all_results <- rbind(all_results, data.frame(
#               scenario = scenario,
#               simulation = sim_num, 
#               method = paste0("TreeTime_", correction_label), 
#               accuracy = treetime_acc$acc, 
#               prob_accuracy = treetime_acc$prob_acc,
#               accuracy_down = treetime_acc$acc_down,
#               prob_accuracy_down = treetime_acc$prob_acc_down,
#               tree_size = tree_size
#             ))
#           }
#         }, error = function(e) cat("TreeTime", correction_label, "failed for sim", sim_num, "\n"))
#       }
#       
#     }, error = function(e) {
#       cat("Failed to process simulation", sim_num, ":", e$message, "\n")
#     })
#   }
# }
# 
# # Optional: Filter by tree size if needed
# all_results <- all_results[all_results$tree_size >= 50, ]
# 
# # Save combined results
# # write.csv(all_results, "combined_accuracy_results.csv", row.names = FALSE)
# 
# #read csv file
# all_results_withoutds$scenario <- "withoutds"
# all_results_withoutds$accuracy_down <- NA
# all_results_withoutds$prob_accuracy_down <- NA
# 
# all_results_total <- rbind(all_results,all_results_withoutds)
# 
# p1 <- ggplot(all_results_total, aes(x = method)) +
#   geom_violin(data = subset(all_results_total, scenario != "baseline"),
#               aes(y = accuracy, fill = "No downsampling"), 
#               alpha = 0.5,
#               position = position_dodge(width = 0.8)) +
#   geom_violin(data = subset(all_results_total, method %in% c("ACE", "SIMMAP","TreeTime_1.0","PastML") & scenario != "withoutds"),
#               aes(y = accuracy_down, fill = "Downsampled with \nmissing transitions"), 
#               alpha = 0.5,
#               position = position_nudge(x = -0.25)) +
#   geom_violin(data = subset(all_results_total, method %in% c("ACE", "SIMMAP","TreeTime_1.0","PastML") & scenario != "withoutds"),
#               aes(y = accuracy, fill = "Downsampled only"), 
#               alpha = 0.5,
#               position = position_nudge(x = 0.25)) +
#   scale_x_discrete(limits = c("ACE", "PastML", "SIMMAP","TreeTime_1.0", "TreeTime_Custom",
#                               "SAASI", "SAASI_pars"),
#                    labels = c("ace", "PastML", "simmap","TreeTime",  "TreeTime\n(Correction)",
#                               "SAASI", "SAASI*")) +
#   scale_fill_manual(values = c("No downsampling" = "lightblue",
#                                "Downsampled with \nmissing transitions" = "lightgreen",
#                                "Downsampled only" = "lightcoral"),
#                     breaks = c("Downsampled only", 
#                                "No downsampling", 
#                                "Downsampled with \nmissing transitions"),
#                     name = "Scenario") +
#   labs(title = "a",
#        x = "", y = "") +
#   theme_minimal() +
#   theme(text = element_text(size = 15, family = "serif"),
#         plot.title = element_text(size = 15),
#         axis.text.x = element_blank(),  # Remove x-axis labels
#         axis.title.y = element_text(hjust = 0.5),
#         legend.title = element_text(size = 15, face = "bold")) +
#   ylim(0, 1.0001)
# 
# p2 <- ggplot(all_results_total, aes(x = method)) +
#   geom_violin(data = subset(all_results_total, scenario != "baseline"),
#               aes(y = prob_accuracy, fill = "No downsampling"), 
#               alpha = 0.5,
#               position = position_dodge(width = 0.8)) +
#   geom_violin(data = subset(all_results_total, method %in% c("ACE", "SIMMAP","TreeTime_1.0","PastML") & scenario != "withoutds"),
#               aes(y = prob_accuracy_down, fill = "Downsampled with \nmissing transitions"), 
#               alpha = 0.5,
#               position = position_nudge(x = -0.25)) +
#   geom_violin(data = subset(all_results_total, method %in% c("ACE", "SIMMAP","TreeTime_1.0","PastML") & scenario != "withoutds"),
#               aes(y = prob_accuracy, fill = "Downsampled only"), 
#               alpha = 0.5,
#               position = position_nudge(x = 0.25)) +
#   scale_x_discrete(limits = c("ACE", "PastML", "SIMMAP","TreeTime_1.0",  "TreeTime_Custom",
#                               "SAASI", "SAASI_pars"),
#                    labels = c("ace", "PastML", "simmap","TreeTime",  "TreeTime\n(Correction)",
#                               "SAASI", "SAASI*")) +
#   scale_fill_manual(values = c("No downsampling" = "lightblue",
#                                "Downsampled with \nmissing transitions" = "lightgreen",
#                                "Downsampled only" = "lightcoral"),
#                     breaks = c("Downsampled only", 
#                                "No downsampling", 
#                                "Downsampled with \nmissing transitions"),
#                     name = "Scenario") +
#   labs(title = "b",
#        x = "", y = "") +
#   theme_minimal() +
#   theme(text = element_text(size = 15, family = "serif"),
#         plot.title = element_text(size = 15),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title.y = element_text(hjust = 0.5),
#         axis.title.x = element_text(hjust = 0.5),
#         legend.title = element_text(size = 15, face = "bold")) +
#   ylim(0, 1.0001)
# 
# p3 <- ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = "bottom")
# 
# 
# # Summary statistics including tree size by scenario
# summary_stats <- all_results %>%
#   group_by(scenario, method) %>%
#   summarise(
#     n = n(),
#     mean_accuracy = mean(accuracy, na.rm = TRUE),
#     sd_accuracy = sd(accuracy, na.rm = TRUE),
#     mean_prob_accuracy = mean(prob_accuracy, na.rm = TRUE),
#     sd_prob_accuracy = sd(prob_accuracy, na.rm = TRUE),
#     mean_tree_size = mean(tree_size, na.rm = TRUE),
#     sd_tree_size = sd(tree_size, na.rm = TRUE),
#     min_tree_size = min(tree_size, na.rm = TRUE),
#     max_tree_size = max(tree_size, na.rm = TRUE),
#     .groups = 'drop'
#   )
# 
# print(summary_stats)
# #write.csv(summary_stats, "combined_accuracy_summary.csv", row.names = FALSE)
# 
# # Overall summary across all scenarios
# overall_summary <- all_results %>%
#   group_by(method) %>%
#   summarise(
#     n = n(),
#     mean_accuracy = mean(accuracy, na.rm = TRUE),
#     sd_accuracy = sd(accuracy, na.rm = TRUE),
#     mean_prob_accuracy = mean(prob_accuracy, na.rm = TRUE),
#     sd_prob_accuracy = sd(prob_accuracy, na.rm = TRUE),
#     mean_tree_size = mean(tree_size, na.rm = TRUE),
#     sd_tree_size = sd(tree_size, na.rm = TRUE),
#     min_tree_size = min(tree_size, na.rm = TRUE),
#     max_tree_size = max(tree_size, na.rm = TRUE),
#     .groups = 'drop'
#   )
# 
# 
# 
# 
# overall_summary <- all_results %>%
#   group_by(method) %>%
#   summarise(
#     n = n(),
#     mean_accuracy = mean(accuracy, na.rm = TRUE),
#     sd_accuracy = sd(accuracy, na.rm = TRUE),
#     mean_prob_accuracy = mean(prob_accuracy, na.rm = TRUE),
#     sd_prob_accuracy = sd(prob_accuracy, na.rm = TRUE),
#     mean_accuracy_down = mean(accuracy_down, na.rm = TRUE),
#     sd_accuracy_down = sd(accuracy_down, na.rm = TRUE),
#     mean_prob_accuracy_down = mean(prob_accuracy_down, na.rm = TRUE),
#     sd_prob_accuracy_down = sd(prob_accuracy_down, na.rm = TRUE),
#     mean_tree_size = mean(tree_size, na.rm = TRUE),
#     sd_tree_size = sd(tree_size, na.rm = TRUE),
#     min_tree_size = min(tree_size, na.rm = TRUE),
#     max_tree_size = max(tree_size, na.rm = TRUE),
#     .groups = 'drop'
#   )
# overall_summary[,c(1,3:10)]


# Save summary tables
# write.csv(overall_summary, "summary_statistics_accuracies_down_onefourth.csv", row.names = FALSE)
# saveRDS(all_results_total,"all_results_withdownsample_onefourth.rds")


# RDS
d <- readRDS("Code/Results/all_results_withdownsample_onefourth.rds")

p1 <- ggplot(d, aes(x = method)) +
  geom_violin(data = subset(d, scenario != "baseline"),
              aes(y = accuracy, fill = "No downsampling"),
              alpha = 0.5,
              position = position_dodge(width = 0.8)) +
  geom_violin(data = subset(d, method %in% c("ACE", "SIMMAP","TreeTime_1.0","PastML") & scenario != "withoutds"),
              aes(y = accuracy_down, fill = "Downsampled with \nmissing transitions"),
              alpha = 0.5,
              position = position_nudge(x = -0.25)) +
  geom_violin(data = subset(d, method %in% c("ACE", "SIMMAP","TreeTime_1.0","PastML") & scenario != "withoutds"),
              aes(y = accuracy, fill = "Downsampled only"),
              alpha = 0.5,
              position = position_nudge(x = 0.25)) +
  scale_x_discrete(limits = c("ACE", "PastML", "SIMMAP","TreeTime_1.0", "TreeTime_Custom",
                              "SAASI", "SAASI_pars"),
                   labels = c("ace", "PastML", "simmap","TreeTime",  "TreeTime\n(Correction)",
                              "SAASI", "SAASI*")) +
  scale_fill_manual(values = c("No downsampling" = "lightblue",
                               "Downsampled with \nmissing transitions" = "lightgreen",
                               "Downsampled only" = "lightcoral"),
                    breaks = c("Downsampled only",
                               "No downsampling",
                               "Downsampled with \nmissing transitions"),
                    name = "Scenario") +
  labs(title = "a",
       x = "", y = "Accuracy") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "serif"),
        plot.title = element_text(size = 15),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.title.y = element_text(hjust = 0.5),
        legend.title = element_text(size = 15, face = "bold")) +
  ylim(0, 1.0001)

p2 <- ggplot(d, aes(x = method)) +
  geom_violin(data = subset(d, scenario != "baseline"),
              aes(y = prob_accuracy, fill = "No downsampling"),
              alpha = 0.5,
              position = position_dodge(width = 0.8)) +
  geom_violin(data = subset(d, method %in% c("ACE", "SIMMAP","TreeTime_1.0","PastML") & scenario != "withoutds"),
              aes(y = prob_accuracy_down, fill = "Downsampled with \nmissing transitions"),
              alpha = 0.5,
              position = position_nudge(x = -0.25)) +
  geom_violin(data = subset(d, method %in% c("ACE", "SIMMAP","TreeTime_1.0","PastML") & scenario != "withoutds"),
              aes(y = prob_accuracy, fill = "Downsampled only"),
              alpha = 0.5,
              position = position_nudge(x = 0.25)) +
  scale_x_discrete(limits = c("ACE", "PastML", "SIMMAP","TreeTime_1.0",  "TreeTime_Custom",
                              "SAASI", "SAASI_pars"),
                   labels = c("ace", "PastML", "simmap","TreeTime",  "TreeTime\n(Correction)",
                              "SAASI", "SAASI*")) +
  scale_fill_manual(values = c("No downsampling" = "lightblue",
                               "Downsampled with \nmissing transitions" = "lightgreen",
                               "Downsampled only" = "lightcoral"),
                    breaks = c("Downsampled only",
                               "No downsampling",
                               "Downsampled with \nmissing transitions"),
                    name = "Scenario") +
  labs(title = "b",
       x = "", y = "Accuracy") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "serif"),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(hjust = 0.5),
        axis.title.x = element_text(hjust = 0.5),
        legend.title = element_text(size = 15, face = "bold")) +
  ylim(0, 1.0001)

p3 <- ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = "bottom")
