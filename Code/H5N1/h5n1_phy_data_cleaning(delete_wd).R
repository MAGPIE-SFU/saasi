library(tidyr)
library(stringr)
require("treeio")

setwd("/Users/yexuan-magpie/Desktop/SAASI_NC_codes/")

options(warn=-1)


tree <- ape::read.nexus('Data/H5N1.tre')
tree$tip.label <- gsub("'", "", tree$tip.label)

tree_info_raw <- tree$tip.label

tree_info <- data.frame()

###################################
# id - the entire name
id <- c()
for(i in 1:range(length(tree_info_raw))){
  id[i] <- tree_info_raw[i]
}

id_df <- data.frame(id)
tree_info <- rbind(tree_info,id_df)

###################################
# location - the second last element before "|"
location <- c()
for(i in 1:range(length(tree_info_raw))){
  location[i] <- tail(str_split(tree_info_raw[i], "\\|")[[1]],2)[1]
}

location_df <- data.frame(location)
tree_info <- cbind(tree_info,location_df)

###################################
# species 
species <- c()
for(i in 1:range(length(tree_info_raw))){
  species[i] <- tail(str_split(tree_info_raw[i], "\\|")[[1]],3)[1]
}

species_df <- data.frame(species)
tree_info <- cbind(tree_info,species_df)

###################################
# segment
segment <- c()
for(i in 1:range(length(tree_info_raw))){
  segment[i] <- "HA"
}

segment_df <- data.frame(segment)
tree_info <- cbind(tree_info,segment_df)

###################################
# strains
strains <- c()
for(i in 1:range(length(tree_info_raw))){
  strains[i] <- str_split(tree_info_raw[i], "\\|")[[1]][5]
}

strains_df <- data.frame(strains)
tree_info <- cbind(tree_info,strains_df)

###################################
# date
date <- c()
for(i in 1:range(length(tree_info_raw))){
  date[i] <- tail(str_split(tree_info_raw[i], "\\|")[[1]],1)
}

date_df <- data.frame(date)
tree_info <- cbind(tree_info,date_df)



###################################
# adding info into the phylogeny 
beast_info <- read.beast('Data/H5N1.tre') #this is a treedata object
beast_df <- as_tibble(beast_info)

internal_info <- data.frame(beast_df$label,beast_df$Host)

tree_withnode <- as_tibble(tree)
tree_withnode$nodelabels = internal_info$beast_df.Host



