# uncomment this to run the full analysis, otherwise use rds file
# 
# library(ape)
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(ggpubr)
# library(rstatix)
# 
# source("Code/Simulation/R/lost_info.R")
# 
# # NOTE: THIS SHOULD BE YOUR FILE DIR ONCE YOU RUN THE SHELL SCRIPT. THIS 
# # CODE SUPPLMENTARY FIGURE 9 & 10.
# base_dir <- "final_equal_sampling/"  
# 
# # Accuracy function
# saasi_accuracy <- function(tree, asr, true.states) {
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
#   
#   return(list(
#     acc = accuracy,
#     prob_acc = prob_accuracy
#   ))
# }
# 
# # Find all simulation directories
# sim_dirs <- list.dirs(base_dir, recursive = FALSE)
# sim_dirs <- sim_dirs[grep("sim_\\d+", sim_dirs)]
# sim_nums <- as.numeric(gsub(".*sim_(\\d+)", "\\1", basename(sim_dirs)))
# 
# all_results <- data.frame()
# all_timing <- data.frame()
# 
# treetime_corrections <- c("1_0", "2_0", "5_0", "10_0")
# treetime_labels <- c("1.0", "2.0", "5.0", "10.0")
# 
# for(i in seq_along(sim_dirs)) {
#   sim_dir <- sim_dirs[i]
#   sim_num <- sim_nums[i]
#   
#   cat("Processing simulation", sim_num, "\n")
#   
#   tryCatch({
#     true_tree <- read.csv(file.path(sim_dir, paste0("true_tree_", sim_num, ".csv")))
#     tree <- read.tree(file.path(sim_dir, paste0("tree_", sim_num, ".nwk")))
#     
#     tree_size <- tree$Nnode * 2 + 1
#     
#     timing_file_1 <- file.path(sim_dir, paste0("timing_", sim_num, ".csv"))
#     if(file.exists(timing_file_1)) {
#       timing_data_1 <- read.csv(timing_file_1)
#       timing_data_1$simulation <- sim_num
#       timing_data_1$tree_size <- tree_size
#       all_timing <- rbind(all_timing, timing_data_1)
#     }
#     
#     timing_file_2 <- file.path(sim_dir, "external_timing.csv")
#     if(file.exists(timing_file_2)) {
#       timing_data_2 <- read.csv(timing_file_2)
#       timing_data_2$tree_size <- tree_size
#       all_timing <- rbind(all_timing, timing_data_2)
#     }
#     
#     # ACE
#     tryCatch({
#       ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
#       ace_acc <- saasi_accuracy(tree, ace_result[2:3], true_tree)
#       all_results <- rbind(all_results, data.frame(
#         simulation = sim_num, method = "ACE", 
#         accuracy = ace_acc$acc, prob_accuracy = ace_acc$prob_acc,
#         tree_size = tree_size
#       ))
#     }, error = function(e) cat("ACE failed for sim", sim_num, "\n"))
#     
#     # SIMMAP
#     tryCatch({
#       simmap_result <- read.csv(file.path(sim_dir, paste0("simmap_", sim_num, ".csv")))
#       simmap_subset <- head(simmap_result[2:3], length(true_tree$State) - length(tree$tip.label))
#       simmap_acc <- saasi_accuracy(tree, simmap_subset, true_tree)
#       all_results <- rbind(all_results, data.frame(
#         simulation = sim_num, method = "SIMMAP", 
#         accuracy = simmap_acc$acc, prob_accuracy = simmap_acc$prob_acc,
#         tree_size = tree_size
#       ))
#     }, error = function(e) cat("SIMMAP failed for sim", sim_num, "\n"))
#     
#     # SAASI
#     tryCatch({
#       saasi_result <- read.csv(file.path(sim_dir, paste0("saasi_", sim_num, ".csv")))
#       saasi_acc <- saasi_accuracy(tree, saasi_result[2:3], true_tree)
#       all_results <- rbind(all_results, data.frame(
#         simulation = sim_num, method = "SAASI", 
#         accuracy = saasi_acc$acc, prob_accuracy = saasi_acc$prob_acc,
#         tree_size = tree_size
#       ))
#     }, error = function(e) cat("SAASI failed for sim", sim_num, "\n"))
#     
#     # SAASI with pars est
#     tryCatch({
#       saasi_pars_result <- read.csv(file.path(sim_dir, paste0("saasi_pars_est_", sim_num, ".csv")))
#       saasi_pars_acc <- saasi_accuracy(tree, saasi_pars_result[2:3], true_tree)
#       all_results <- rbind(all_results, data.frame(
#         simulation = sim_num, method = "SAASI_pars", 
#         accuracy = saasi_pars_acc$acc, prob_accuracy = saasi_pars_acc$prob_acc,
#         tree_size = tree_size
#       ))
#     }, error = function(e) cat("SAASI_pars failed for sim", sim_num, "\n"))
#     
#     # PastML
#     tryCatch({
#       pastml_file <- file.path(sim_dir, "pastml_results", 
#                                "marginal_probabilities.character_numeric_state.model_F81.tab")
#       if(file.exists(pastml_file)) {
#         pastml_result <- read.table(pastml_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#         pastml_subset <- pastml_result[grep("^(root|nd\\d+)$", pastml_result$node), ]
#         ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
#         pastml_matched <- pastml_subset[match(ace_result$X, pastml_subset$node), ]
#         pastml_acc <- saasi_accuracy(tree, pastml_matched[2:3], true_tree)
#         all_results <- rbind(all_results, data.frame(
#           simulation = sim_num, method = "PastML", 
#           accuracy = pastml_acc$acc, prob_accuracy = pastml_acc$prob_acc,
#           tree_size = tree_size
#         ))
#       }
#     }, error = function(e) cat("PastML failed for sim", sim_num, "\n"))
#     
#     # TreeTime with different sampling bias corrections
#     for(j in seq_along(treetime_corrections)) {
#       correction_suffix <- treetime_corrections[j]
#       correction_label <- treetime_labels[j]
#       
#       tryCatch({
#         treetime_file <- file.path(sim_dir, "treetime_results", "confidence.csv")
#         if(file.exists(treetime_file)) {
#           treetime_result <- read.csv(treetime_file)
#           treetime_subset <- treetime_result[grep("^(nd\\d+)$", treetime_result$X.name), ]
#           ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
#           treetime_matched <- treetime_subset[match(ace_result$X, treetime_subset$X.name), ]
#           treetime_acc <- saasi_accuracy(tree, treetime_matched[2:3], true_tree)
#           all_results <- rbind(all_results, data.frame(
#             simulation = sim_num, 
#             method = paste0("TreeTime_", correction_label), 
#             accuracy = treetime_acc$acc, 
#             prob_accuracy = treetime_acc$prob_acc,
#             tree_size = tree_size
#           ))
#         }
#       }, error = function(e) cat("TreeTime", correction_label, "failed for sim", sim_num, "\n"))
#     }
#     
#   }, error = function(e) {
#     cat("Failed to process simulation", sim_num, ":", e$message, "\n")
#   })
# }
# 
# # Filter by tree size
# all_results <- all_results[all_results$tree_size >= 50, ]
# all_timing <- all_timing[all_timing$tree_size >= 50, ]
# all_timing <- all_timing[all_timing$time_seconds <= 50, ]
# 
# # Process timing data - standardize method names
# all_timing <- all_timing %>%
#   mutate(method = case_when(
#     grepl("treetime_1\\.0", method, ignore.case = TRUE) ~ "TreeTime_1.0",
#     grepl("treetime_2\\.0", method, ignore.case = TRUE) ~ "TreeTime_2.0",
#     grepl("treetime_5\\.0", method, ignore.case = TRUE) ~ "TreeTime_5.0",
#     grepl("treetime_10\\.0", method, ignore.case = TRUE) ~ "TreeTime_10.0",
#     grepl("^treetime$", method, ignore.case = TRUE) ~ "TreeTime_1.0",  
#     grepl("pastml", method, ignore.case = TRUE) ~ "PastML",
#     grepl("ace", method, ignore.case = TRUE) ~ "ACE",
#     grepl("simmap", method, ignore.case = TRUE) ~ "SIMMAP",
#     grepl("est", method, ignore.case = TRUE) ~ "SAASI_pars",
#     grepl("saasi", method, ignore.case = TRUE) ~ "SAASI",
#     TRUE ~ method
#   ))
# 
# all_results_filtered <- all_results %>% 
#   filter(!method %in% c("TreeTime_2.0", "TreeTime_5.0", "TreeTime_10.0"))
# 
# d1 <- all_results %>% filter(method %in% "ACE")
# d2 <- all_results %>% filter(method %in% "PastML")
# d3 <- all_results %>% filter(method %in% "SIMMAP")
# d4 <- all_results %>% filter(method %in% "TreeTime_1.0")
# d5 <- all_results %>% filter(method %in% "SAASI")
# d6 <- all_results %>% filter(method %in% "SAASI_pars")
# 
# ks.test(d5$accuracy,d1$accuracy)
# ks.test(d5$accuracy,d2$accuracy)
# ks.test(d5$accuracy,d3$accuracy)
# ks.test(d5$accuracy,d4$accuracy)
# ks.test(d5$prob_accuracy,d1$prob_accuracy)
# ks.test(d5$prob_accuracy,d2$prob_accuracy)
# ks.test(d5$prob_accuracy,d3$prob_accuracy)
# ks.test(d5$prob_accuracy,d4$prob_accuracy)
# 
# 
# ks.test(d6$accuracy,d1$accuracy)
# ks.test(d6$accuracy,d2$accuracy)
# ks.test(d6$accuracy,d3$accuracy)
# ks.test(d6$accuracy,d4$accuracy)
# ks.test(d6$prob_accuracy,d1$prob_accuracy)
# ks.test(d6$prob_accuracy,d2$prob_accuracy)
# ks.test(d6$prob_accuracy,d3$prob_accuracy)
# ks.test(d6$prob_accuracy,d4$prob_accuracy)
# 
# # Create comprehensive summary statistics
# create_summary_table <- function(data, value_col) {
#   methods <- c("ACE", "PastML", "SIMMAP", "TreeTime_1.0", "SAASI", "SAASI_pars")
#   
#   summary <- data %>%
#     filter(method %in% methods) %>%
#     group_by(method) %>%
#     summarise(
#       Median = round(median(!!sym(value_col)), 4),
#       Mean = round(mean(!!sym(value_col)), 4),
#       SD = round(sd(!!sym(value_col)), 4),
#       Min = round(min(!!sym(value_col)), 4),
#       Max = round(max(!!sym(value_col)), 4),
#       Q25 = round(quantile(!!sym(value_col), 0.25), 4),
#       Q75 = round(quantile(!!sym(value_col), 0.75), 4),
#       IQR = round(IQR(!!sym(value_col)), 4),
#       .groups = 'drop'
#     ) %>%
#     arrange(desc(Median))
#   
#   return(summary)
# }
# 
# # Create summary tables
# summary_acc <- create_summary_table(all_results, "accuracy")
# 
# summary_prob <- create_summary_table(all_results, "prob_accuracy")
# 
# # write.csv(summary_acc, "summary_statistics_accuracy_er.csv", row.names = FALSE)
# # write.csv(summary_prob, "summary_statistics_prob_accuracy_er.csv", row.names = FALSE)

# p1 <- ggplot(all_results_filtered, aes(x = method, y = accuracy)) +
#   geom_violin(alpha = 0.5, fill = "lightblue") +
#   scale_x_discrete(limits = c("ACE", "PastML","SIMMAP", "TreeTime_1.0", 
#                               "SAASI", "SAASI_pars"),
#                    labels = c("ace", "PastML", "simmap","TreeTime",
#                               "SAASI", "SAASI*")) +
#   labs(title = "a", x = "", y = "Accuracy") +
#   theme_minimal() +
#   theme(text = element_text(size = 15, family = "serif"),
#         plot.title = element_text(size = 15),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title.y = element_text(hjust = 0.5),  # Center the y-axis label
#         legend.title = element_text(size = 15, face = "bold")) +
#   ylim(0.5, 1)
# 
# p2 <- ggplot(all_results_filtered, aes(x = method, y = prob_accuracy)) +
#   geom_violin(alpha = 0.5, fill = "lightcoral") +
#   scale_x_discrete(limits = c("ACE", "PastML", "SIMMAP",
#                               "TreeTime_1.0",
#                               "SAASI", "SAASI_pars"),
#                    labels = c("ace", "PastML","simmap",
#                               "TreeTime",
#                               "SAASI", "SAASI*")) +
#   labs(title = "b",
#        x = "", y = "") +
#   theme_minimal() +
#   theme(text = element_text(size = 15, family = "serif", hjust = 1),
#         plot.title = element_text(size = 15),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.title = element_text(size = 15, face = "bold")) +
#   ylim(0.5, 1)
# 
# p3 <- ggarrange(p1, p2, nrow = 1)
# 
# 
# all_timing <- all_timing %>% filter(time_seconds > 0)
# 
# new_timing <- all_timing %>% 
#   filter(!(method %in% c("SAASI", "SAASI_pars", "SIMMAP") & time_seconds < 10))
# 
# timing_results_filtered <- new_timing  %>% 
#   filter(!method %in% c("TreeTime_2.0", "TreeTime_5.0", "TreeTime_10.0"))
# 
# timing_results_filtered$log_time <- log10(timing_results_filtered$time_seconds)
# p4 <- ggplot(timing_results_filtered, aes(x = method, y = time_seconds)) +
#   geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.3, size = 2) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#   scale_color_distiller(palette = "Spectral", direction = 1)+
#   scale_x_discrete(limits = c("ACE", "PastML", "SIMMAP",
#                               "TreeTime_1.0",
#                               "SAASI", "SAASI_pars"),
#                    labels = c("ace", "PastML", 
#                               "simmap","TreeTime", 
#                               "SAASI", "SAASI*")) +
#   labs(title = "",
#        x = "", y = "Time (seconds)",
#        color = "Tree size") +
#   theme_minimal() + 
#   theme(text = element_text(size = 15, family = "serif"),
#         plot.title = element_text(size = 15),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title.y = element_text(hjust = 0.5),
#         axis.title.x = element_text(hjust = 0.5),
#         legend.title = element_text(size = 15, face = "bold"))
# 
# 
# 
# overall_summary <- all_results_filtered %>%
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
# overall_summary
# # Save summary tables
# 
# # saveRDS(timing_results_filtered,"1000_sims_er_time.rds")
# # saveRDS(all_results_filtered,"1000_sims_er.rds")
# 
# 


#######
#RDS
d <- readRDS("Code/Results/1000_sims_er.rds")
dt <- readRDS("Code/Results/1000_sims_er_time.rds")


p1 <- ggplot(d, aes(x = method, y = accuracy)) +
  geom_violin(alpha = 0.5, fill = "lightblue") +
  scale_x_discrete(limits = c("ACE", "PastML","SIMMAP", "TreeTime_1.0",
                              "SAASI", "SAASI_pars"),
                   labels = c("ace", "PastML", "simmap","TreeTime",
                              "SAASI", "SAASI*")) +
  labs(title = "a", x = "", y = "Accuracy") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "serif"),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(hjust = 0.5),  # Center the y-axis label
        legend.title = element_text(size = 15, face = "bold")) +
  ylim(0.5, 1)

p2 <- ggplot(d, aes(x = method, y = prob_accuracy)) +
  geom_violin(alpha = 0.5, fill = "lightcoral") +
  scale_x_discrete(limits = c("ACE", "PastML", "SIMMAP",
                              "TreeTime_1.0",
                              "SAASI", "SAASI_pars"),
                   labels = c("ace", "PastML","simmap",
                              "TreeTime",
                              "SAASI", "SAASI*")) +
  labs(title = "b",
       x = "", y = "") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "serif", hjust = 1),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 15, face = "bold")) +
  ylim(0.5, 1)

p3 <- ggarrange(p1, p2, nrow = 1)

p4 <- ggplot(dt, aes(x = method, y = time_seconds)) +
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.3, size = 2) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  scale_color_distiller(palette = "Spectral", direction = 1)+
  scale_x_discrete(limits = c("ACE", "PastML", "SIMMAP",
                              "TreeTime_1.0",
                              "SAASI", "SAASI_pars"),
                   labels = c("ace", "PastML",
                              "simmap","TreeTime",
                              "SAASI", "SAASI*")) +
  labs(title = "",
       x = "", y = "Time (seconds)",
       color = "Tree size") +
  theme_minimal() +
  theme(text = element_text(size = 15, family = "serif"),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(hjust = 0.5),
        axis.title.x = element_text(hjust = 0.5),
        legend.title = element_text(size = 15, face = "bold"))
