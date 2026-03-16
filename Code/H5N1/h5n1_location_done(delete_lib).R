# Generate Figure 5

# library("ape")
# library("deSolve")
# library("strap")
# library("stringr")
# library("ips")
# require("treeio")
# library("ggplot2")
# require("ggtree")
# library("reshape2")
# library("ggstance")
# library("diversitree")
# library("ggimage")
# library("ggpubr")
# library("dplyr")
# library("patchwork")
# library("tidyverse")
# library("janitor")
# library("stringr")
# library("skimr")
# library("phytools")
# library(dplyr)
# library(rstatix)
# library(jsonlite)

source(file.path("SAASI/ode_solve.R"))
source(file.path("SAASI/saasi.R"))
source(file.path("Code/H5N1/h5n1_phy_data_cleaning.R"))

####################
# unmodified phylogeny
unmod <- ggtree(tree)

# modified phylogeny
tree_modified <- di2multi(tree, tol = 1e-08)

drop_tip <- c(1,104)
tree_modified <- drop.tip(tree_modified,drop_tip)
tree_info <- tree_info[-c(1,104),]


# modify the tree as well by deleting human tip, only one, and
# deleting the root ace estimates complex transition rates.

tips_to_drop <- grep("human", tree_modified$tip.label, value = TRUE, ignore.case = TRUE)
tree_modified <- drop.tip(tree_modified, tips_to_drop)
tree_info <- tree_info %>% 
  filter(species != "human")

# plot the phylogeny and color the tips based on species
phy_species <- ggtree(tree_modified)
phy_species <- phy_species %<+% tree_info + geom_tippoint(aes(color=species))

# plot the phylogeny and color the tips based on locations
phy_locations <- ggtree(tree_modified)
phy_locations <- phy_locations %<+% tree_info + geom_tippoint(aes(color=location))

# plot the phylogeny and color all the nodes based on beast phy
phy_beast <- ggtree(tree_modified)
phy_beast <- phy_beast %<+% tree_withnode + geom_point(aes(color=nodelabels))


####################
# Analysis (Ancestral State Reconstruction ASR)
####################

####################
# ASR using ace
ace_phy <- multi2di(tree_modified)
ace_phy$state <- tree_info$location

ans_ace<-ace(ace_phy$state, ace_phy,type = "discrete",method = "ML", model="ER")

# ace estimate the transition rate
qij_rate <- ans_ace$rates


####################
# saasi analysis

replace_matrix_with_vector <- function(matrix, vector) {
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- vector[matrix[i,j]]
    }
  }
  return(matrix)
}

qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)
mm<-qij_matrix

pars <- c(21.1,  21.1, 21.1, 21.1, 21.1, 21.1,21.1,21.1,21.1,21.1,21.1,21.1,21.1,    
          6.7, 6.7, 6.7, 6.7, 6.7, 6.7,6.7,6.7,6.7,6.7,6.7,6.7,6.7, 
          5,  5,  5,  5,  5,  5, 5,5,5,5,5,5,5) 

vect <- rep(0.1,13*13-13)
pars <- c(pars,vect)

our_phy <- multi2di(tree_modified)
our_phy$tip.label <- tree_info$location
lab <- our_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab
our_phy$tip.state <- numeric_vector

asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,10,11,12,13,2,3,4,5,6,7,8,9)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:13,color=c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
  "#FFD92F", "#E5C494", "#B3B3B3"
))

our_species <- ggtree(ace_phy,mrsd = "2024-04-06",options(ignore.negative.edge=TRUE))
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=location),size=0.9)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
    "#FFD92F", "#E5C494", "#B3B3B3"
  ),labels = c("New_Mexico" = "New Mexico", "North_Carolina" = "North Carolina"))+
  theme_tree2()+ggtitle("a") +labs(color = "Location")+
  geom_hilight(node=120, fill="purple",alpha=0.05) +
  geom_vline(xintercept = 2023.97, color = "red", linetype = "dashed", size = 0.5)+
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15), 
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))

p1 <- inset(our_species, our_pie,width = 0.03,height = 0.03,hjust=0.005)


############################
# now plot the alluvial plot
# Extract node states

asr_our <- saasi(pars,our_phy,mm)

node_states <- apply(asr_our, 1, which.max)
node_states <- factor(node_states, levels = 1:13, labels = c(
  "California", "Indiana", "Kansas","Maryland","Michigan","Minnesota","Montana",
  "New_Mexico","North_Carolina","Ohio","Oklahoma","Texas","Wisconsin"))

all_states <- c(ace_phy$state, node_states)
# 
# Create a data frame for plotting
plot_data <- data.frame(
  node = c(1:length(ace_phy$state), 1:ace_phy$Nnode + length(ace_phy$state)),
  descendant = all_states,
  is_tip = c(rep(TRUE, length(ace_phy$state)), rep(FALSE, ace_phy$Nnode))
)
# 
# Add parent information
plot_data$parent <- sapply(plot_data$node, function(n) {
  if (n == length(ace_phy$state) + 1) return(NA)  # root node
  parent <- ace_phy$edge[ace_phy$edge[,2] == n, 1]
  if (length(parent) == 0) return(NA)  # tip
  return(parent)
})


node_to_node <- plot_data[!plot_data$is_tip & !is.na(plot_data$parent),]
node_to_node$ancestor <- plot_data$descendant[match(node_to_node$parent, plot_data$node)]

node_to_tip <- plot_data[plot_data$is_tip & !is.na(plot_data$parent),]
node_to_tip$ancestor <- plot_data$descendant[match(node_to_tip$parent, plot_data$node)]

alluvial_data <- rbind(node_to_node, node_to_tip)

alluvial_long <- pivot_longer(alluvial_data, cols = c(descendant, ancestor),
                              names_to = "generation", values_to = "location")

unique_location <- unique(alluvial_long$location)
lo <- c(
  "California", "Indiana", "Kansas", "Maryland", "Michigan", "Minnesota", "Montana",
  "New_Mexico", "North_Carolina", "Ohio", "Oklahoma", "Texas", "Wisconsin")

# incorrect order
alluvial_long$location[alluvial_long$location == lo[1]] <- 1
alluvial_long$location[alluvial_long$location == lo[2]] <- 2
alluvial_long$location[alluvial_long$location == lo[3]] <- 3
alluvial_long$location[alluvial_long$location == lo[4]] <- 4
alluvial_long$location[alluvial_long$location == lo[5]] <- 5
alluvial_long$location[alluvial_long$location == lo[6]] <- 6
alluvial_long$location[alluvial_long$location == lo[7]] <- 7
alluvial_long$location[alluvial_long$location == lo[8]] <- 8
alluvial_long$location[alluvial_long$location == lo[9]] <- 9
alluvial_long$location[alluvial_long$location == lo[10]] <- 10
alluvial_long$location[alluvial_long$location == lo[11]] <- 11
alluvial_long$location[alluvial_long$location == lo[12]] <- 12
alluvial_long$location[alluvial_long$location == lo[13]] <- 13


# Create the alluvial plot
alluvial_p1 <- ggplot(alluvial_long,
                       aes(x = generation, stratum = location, alluvium = node, 
                           fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Location") +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  ggtitle("c") +
  scale_fill_manual(values = c(
    "#E41A1C", "#66C2A5", "#FC8D62", "#8DA0CB", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#FFFF33", "#F781BF", "#999999","#1B9E77" ,"#D95F02"),
    labels=c(lo[1], lo[10], lo[11], lo[12],lo[13], lo[2], lo[3], lo[4], lo[5], lo[6],lo[7], lo[8], lo[9])
  ) +
  theme(text = element_text(size = 15, family = "serif"), plot.title = element_text(size=15),legend.position = "none")


#######################
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)
mm<-qij_matrix

pars <- c(21.1,  21.1, 21.1, 21.1, 21.1, 21.1,21.1,21.1,21.1,21.1,21.1,21.1,21.1,    
          6.7, 6.7, 6.7, 6.7, 6.7, 6.7,6.7,6.7,6.7,6.7,6.7,6.7,6.7, 
          2,  2,  2,  2,  2,  2, 2,2,2,2,2,10,2) 
vect <- rep(0.1,13*13-13)
pars <- c(pars,vect)
our_phy <- multi2di(tree_modified)
our_phy$tip.label <- tree_info$location
lab <- our_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab
our_phy$tip.state <- numeric_vector
asr_our <- saasi(pars,our_phy,mm)

# internal_node_ids <- root_node:nnode

colnames(asr_our) = c(1,10,11,12,13,2,3,4,5,6,7,8,9)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:13,color=c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
  "#FFD92F", "#E5C494", "#B3B3B3"
))

our_species <- ggtree(ace_phy,mrsd = "2024-04-06",options(ignore.negative.edge=TRUE))
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=location),size=0.9)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
    "#FFD92F", "#E5C494", "#B3B3B3"
  ),labels = c("New_Mexico" = "New Mexico", "North_Carolina" = "North Carolina"))+
  theme_tree2()+ggtitle("b") +labs(color = "Location")+
  geom_hilight(node=120, fill="purple",alpha=0.05) +
  geom_vline(xintercept = 2023.97, color = "red", linetype = "dashed", size = 0.5)+
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15), 
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))


p2 <- inset(our_species, our_pie,width = 0.03,height = 0.03,hjust=0.005)

############################
asr_our <- saasi(pars,our_phy,mm)

node_states <- apply(asr_our, 1, which.max)
node_states <- factor(node_states, levels = 1:13, labels = c(
  "California", "Indiana", "Kansas", "Maryland", "Michigan", "Minnesota", "Montana",
  "New_Mexico", "North_Carolina", "Ohio", "Oklahoma", "Texas", "Wisconsin"))
all_states <- c(ace_phy$state, node_states)

# Create a data frame for plotting
plot_data <- data.frame(
  node = c(1:length(ace_phy$state), 1:ace_phy$Nnode + length(ace_phy$state)),
  descendant = all_states,
  is_tip = c(rep(TRUE, length(ace_phy$state)), rep(FALSE, ace_phy$Nnode))
)

# Add parent information
plot_data$parent <- sapply(plot_data$node, function(n) {
  if (n == length(ace_phy$state) + 1) return(NA)  # root node
  parent <- ace_phy$edge[ace_phy$edge[,2] == n, 1]
  if (length(parent) == 0) return(NA)  # tip
  return(parent)
})

node_to_node <- plot_data[!plot_data$is_tip & !is.na(plot_data$parent),]
node_to_node$ancestor <- plot_data$descendant[match(node_to_node$parent, plot_data$node)]

node_to_tip <- plot_data[plot_data$is_tip & !is.na(plot_data$parent),]
node_to_tip$ancestor <- plot_data$descendant[match(node_to_tip$parent, plot_data$node)]

alluvial_data <- rbind(node_to_node, node_to_tip)

alluvial_long <- pivot_longer(alluvial_data, cols = c(descendant, ancestor),
                              names_to = "generation", values_to = "location")

unique_location <- unique(alluvial_long$location)
lo <- c(
  "California", "Indiana", "Kansas", "Maryland", "Michigan", "Minnesota", "Montana",
  "New_Mexico", "North_Carolina", "Ohio", "Oklahoma", "Texas", "Wisconsin")

# incorrect order
alluvial_long$location[alluvial_long$location == lo[1]] <- 1
alluvial_long$location[alluvial_long$location == lo[2]] <- 2
alluvial_long$location[alluvial_long$location == lo[3]] <- 3
alluvial_long$location[alluvial_long$location == lo[4]] <- 4
alluvial_long$location[alluvial_long$location == lo[5]] <- 5
alluvial_long$location[alluvial_long$location == lo[6]] <- 6
alluvial_long$location[alluvial_long$location == lo[7]] <- 7
alluvial_long$location[alluvial_long$location == lo[8]] <- 8
alluvial_long$location[alluvial_long$location == lo[9]] <- 9
alluvial_long$location[alluvial_long$location == lo[10]] <- 10
alluvial_long$location[alluvial_long$location == lo[11]] <- 11
alluvial_long$location[alluvial_long$location == lo[12]] <- 12
alluvial_long$location[alluvial_long$location == lo[13]] <- 13


# Create the alluvial plot
alluvial_p2 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, 
                            fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Location") +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  ggtitle("d") +
  scale_fill_manual(values = c(
    "#E41A1C", "#66C2A5", "#FC8D62", "#8DA0CB", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#FFFF33", "#F781BF", "#999999","#1B9E77" ,"#D95F02"),
    labels=c(lo[1], lo[10], lo[11], lo[12],lo[13], lo[2], lo[3], lo[4], lo[5], lo[6],lo[7], lo[8], lo[9])
  ) +
  theme(text = element_text(size = 15, family = "serif"), plot.title = element_text(size= 15),legend.position = "none")


p_comb <- ggarrange(p1,p2,alluvial_p1,alluvial_p2,heights = c(1.75,1), nrow=2,ncol = 2, common.legend = TRUE,legend = "bottom")
p_comb

