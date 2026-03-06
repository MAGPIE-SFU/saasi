#' Plot a phylogenetic tree with SAASI tip annotations
#'
#' `plot_saasi()` plots a phylogenetic tree annoted with [saasi()] output. Tips are annotated with their states.
#' Internal nodes are annotated with a pie chart showing the marginal distribution over states for each node.
#' The [ape::plot.phylo()] function is used for plotting.
#'
#' @param tree An object of class `phylo` with `tip.state` attribute. 
#' @param saasi_result Output of `saasi()` containing internal node annotations.
#' @param colours Optional character vector specifying the node colour for each state. By default, colours are automatically
#'   generated based on the number of states.
#' @param tip_cex Size of the tips. Default value is 0.5.
#' @param node_cex Size of node pie charts. Default value is 0.2.
#' @param save_file Optional string specifying the file path where the plot will be saved.
#'   By default the plot is displayed on the current device.
#' @param width Width of the plot in pixels. Only relevant if saving to file. Default value is 3000 pixels.
#' @param height Height of the plot in pixels. Only relevant if saving to file. Default value is 3000 pixels.
#' @param res Plot resolution in dpi. Only relevant if saving to file. Default value is 300 dpi.
#' @return If a save file is specified, the plot will be saved as directed. 
#' If not, then the plot will be displayed on the users device.
#' The function returns the output of `plot.phylo()`, a list containing plot specifications.
#' @export
#' @examples
#' # Run SAASI
#' saasi_res <- saasi(demo_tree_prepared, demo_Q, demo_pars)
#' 
#' # Plot the results using the default settings
#' plot_saasi(demo_tree_prepared, saasi_res)
#'
#' \dontrun{
#' # Use custom colours and save the result to tree.png
#' plot_saasi(tree, saasi_result,
#'            colours = viridis::viridis(3, begin=0.1, end=0.9),
#'            node_cex = 0.3,
#'            save_file = "tree.png")
#' }
#' 
plot_saasi <- function(tree,
                       saasi_result,
                       colours = NULL,
                       tip_cex = 0.5,
                       node_cex = 0.2,
                       save_file = NULL,
                       width = 3000, 
                       height = 3000,
                       res = 300) {
  n_states <- length(unique(tree$tip.state))
  
  if(is.null(colours)){
    colours <- viridis::cividis(n_states, begin=0.1, end=0.9)
  }
  
  if(length(colours) < n_states){
    stop("You must specify at least as many colours as there are states.")
  }
  
  # Map tip states to colours
  state_factor <- factor(tree$tip.state)
  tip_colours <- colours[as.numeric(state_factor)]
  pch_value <- 21
  
  # Open file device if saving
  if(!is.null(save_file)){
    grDevices::png(save_file, width = width, height = height, res = res)
  }
  
  # Plotting 
  plot_phy <- ape::ladderize(tree, right = FALSE)
  res <- plot(plot_phy, show.tip.label = FALSE)
  
  ape::tiplabels(bg = tip_colours, cex = tip_cex, adj = 0.5, pch = pch_value)
  
  if(!is.null(saasi_result)){
    ape::nodelabels(pie = saasi_result,
                    piecol = colours[1:ncol(saasi_result)],
                    cex = node_cex)
  }
  
  graphics::legend("topleft",
                   legend = levels(state_factor),
                   col = colours[1:n_states],
                   pch = 19,
                   pt.cex = 2,
                   bty = "n",
                   cex = 0.8)
  
  if(!is.null(save_file)){
    grDevices::dev.off()
  }
  return(invisible(res))
}