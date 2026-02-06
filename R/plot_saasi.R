#' Plot phylogenetic tree with estimated ancestral states from SAASI
#'
#' @param tree A phylogenetic tree of class `phylo` prepared using `prepare_tree_for_saasi`. This is the same tree as in the `saasi` input. 
#' @param saasi_result Output from `saasi`.
#' @param colors Character vector of colors. If NULL, colors are automatically
#'   generated based on the number of states (default: NULL)
#' @param tip_cex Size of tip circles (default: 0.5)
#' @param node_cex Size of node pie charts (default: 0.2)
#' @param save_file Character. File path to save the plot (e.g. "tree.png").
#'   If NULL, plot is drawn to current device (default: NULL)
#' @param width Plot width in pixels when saving (default: 3000)
#' @param height Plot height in pixels when saving (default: 3000)
#' @param res Plot resolution in dpi when saving (default: 300)
#' @return Invisibly returns the color vector used
#' @export
#' @examples
#' # Basic plot
#' plot_saasi(tree, saasi_result)
#'
#' # Custom colors and save to file
#' plot_saasi(tree, saasi_result,
#'            colors = c("red", "blue", "green"),
#'            node_cex = 0.3,
#'            save_file = "tree.png")
#' 
plot_saasi <- function(tree,
                       saasi_result,
                       colors = NULL,
                       tip_cex = 0.5,
                       node_cex = 0.2,
                       save_file = NULL,
                       width = 3000, 
                       height = 3000,
                       res = 300) {
  n_states <- length(unique(tree$tip.state))
  
  if(is.null(colors)){
    colors <- grDevices::rainbow(n_states)
  }
  
  if(length(colors) < n_states){
    stop("Number of colors must be >= number of states")
  }
  
  # Map tip states to colors
  state_factor <- factor(tree$tip.state)
  tip_colors <- colors[as.numeric(state_factor)]
  pch_value <- 21
  
  # Open file device if saving
  if(!is.null(save_file)){
	  grDevices::png(save_file, width = width, height = height, res = res)
  }
  
  # Plotting 
  plot_phy <- ape::ladderize(tree, right = FALSE)
  plot(plot_phy, show.tip.label = FALSE)
  
  ape::tiplabels(bg = tip_colors, cex = tip_cex, adj = 0.5, pch = pch_value)
  
  ape::nodelabels(pie = saasi_result,
                  piecol = colors[1:ncol(saasi_result)],
                  cex = node_cex)
  
  graphics::legend("topleft",
         legend = levels(state_factor),
         col = colors[1:n_states],
         pch = 19,
         pt.cex = 2,
         bty = "n",
         cex = 0.8)
  
  if(!is.null(save_file)){
	  grDevices::dev.off()
  }
  invisible(colors)
}

