#' Copy/pasted from private `diversitree::make.tree.bisse` function, with some
#' edits to incorporate sampling events and allow three or more states.
#'
#' Adding sampling events throws a cyclomatic complexity error, but I will
#' ignore that because I do not want to touch the original `diversitree` code.
#' Similar thing for lack of snake_case variable names.
#'
#' I did introduce less startling changes to the original code re: other linting
#' errors.
#'
#' @noRd
make.tree.bisse_modified <- function(pars, k, max.taxa = Inf, max.t = Inf, x0) {
  # Matrix representation of state changes to consider
  to <- matrix(unlist(lapply(1 : k, function(i) (1 : k)[-i])), k, k - 1, TRUE)

  extinct <- FALSE
  split   <- FALSE
  # Adding Sampling
  sam <- FALSE

  parent <- 0
  n.i <- rep(0, k)
  r.i <- rowSums(pars)
  len <- 0
  t <- 0
  hist <- list()

  # Always single lineage in our code
  states <- x0
  # We use 1-based states
  n.taxa <- lineages <- n.i[x0] <- 1
  start <- 0

  while (n.taxa <= max.taxa && n.taxa > 0) {
    ## When does an event happen?
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    if (t > max.t) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }

    len[lineages] <- len[lineages] + dt

    # ADDING A NEW EVENT TYPE: SAMPLING
    # type: 1: speciation, 2: extinction, 3: sampling, >3: char change
    # sampling assumes the sample is tested and hence becomes a taxa
    # in the phylogeny, this allows the phylogeny be non-ultrametric
    state <- sample(k, 1, FALSE, r.n / r.tot)

    ## Pick a lineage for that state:
    j <- sample(n.i[state], 1)
    lineage <- lineages[states[lineages] == state][j]
    type <- sample(3 + (k - 1), 1, FALSE, pars[state, ])

    if (type == 1) {
      ## Speciating:
      if (n.taxa == max.taxa)
        ## Don't add this one.
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      sam[new.i] <- FALSE
      states[new.i] <- state
      parent[new.i] <- lineage
      start[new.i] <- t
      len[new.i] <- 0

      n.i[state] <- n.i[state] + 1
      n.taxa <- n.taxa + 1

      lineages <- which(!split & !extinct & !sam)
    } else if (type == 2) {
      ## Extinct
      extinct[lineage] <- TRUE

      lineages <- which(!split & !extinct & !sam)

      n.i[state] <- n.i[state] - 1
      n.taxa <- n.taxa - 1
    } else if (type == 3) {
      ## Sampling; this is new
      sam[lineage] <- TRUE
      n.taxa <- n.taxa - 1
      n.i[state] <- n.i[state] - 1
      lineages <- which(!split & !extinct & !sam)
    } else {
      ## Character switch:
      # Must account for having more than one possible state
      states[lineage] <- state.new <- to[state, type - 3]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1, -1)
      hist[[length(hist) + 1]] <- c(lineage, t, state, state.new)
    }
  }

  info <- data.frame(idx = seq_along(extinct), len = len, parent = parent,
                     start = start, state = states, extinct = extinct,
                     split = split, sam = sam)

  hist <- as.data.frame(do.call(rbind, hist))
  if (nrow(hist) == 0)
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0

  attr(info, "t") <- t
  attr(info, "hist") <- hist
  info
}

#' Copy/pasted from private `diversitree::me.to.ape.bisse` function.
#'
#' I did introduce minor changes to the original code to fix linting issues.
#' However, I did not change variable names to snake_case.
#'
#' @noRd
me.to.ape.bisse <- function(x, root.state) {
  if (nrow(x) == 0)
    return(NULL)
  Nnode <- sum(!x$split) - 1
  n.tips <- sum(!x$split)

  x$idx2 <- NA
  x$idx2[!x$split] <- 1:n.tips
  x$idx2[x$split] <- order(x$idx[x$split]) + n.tips + 1

  i <- match(x$parent, x$idx)
  x$parent2 <- x$idx2[i]
  x$parent2[is.na(x$parent2)] <- n.tips + 1

  tip.label <- ifelse(subset(x, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)

  x$name <- NA
  x$name[!x$split] <- tip.label
  ## More useful, but I don't want to clobber anything...
  x$name2 <- c(tip.label, node.label)[x$idx2]

  tip.state <- x$state[match(1:n.tips, x$idx2)]
  names(tip.state) <- tip.label

  node.state <- x$state[match(1:Nnode + n.tips, x$idx2)]
  names(node.state) <- node.label
  node.state["nd1"] <- root.state

  hist <- attr(x, "hist")
  if (!is.null(hist)) {
    hist$idx2 <- x$idx2[match(hist$idx, x$idx)]
    hist$name2 <- x$name2[match(hist$idx, x$idx)]
    if (nrow(hist) > 0)
      hist <- hist[order(hist$idx2), ]
  }

  phy <- ape::reorder.phylo(structure(
    list(
      edge = cbind(x$parent2, x$idx2),
      Nnode = Nnode,
      tip.label = tip.label,
      tip.state = tip.state,
      node.label = node.label,
      node.state = node.state,
      edge.length = x$len,
      orig = x,
      hist = hist
    ),
    class = "phylo"
  ))
  phy$edge.state <- x$state[match(phy$edge[, 2], x$idx2)]
  phy
}

#' Simulate a birth/death/sampling tree
#'
#' This tree can be passed to [asrproject::saasi], and will include speciation,
#' extinction, sampling, and mutation events.
#'
#' @param x0 Natural number used as root state in the returned tree. Must be a
#' state declared in `params_df`.
#' @param max_taxa Maximum number of nodes allowed in the returned tree.
#' @param max_t Maximum depth allowed in the returned tree.
#' @param include_extinct Boolean declaring whether extinct taxa are included in
#' the returned tree.
#' @return A `phylo` phylogenetic tree (`ape` format).
#' @inheritParams saasi
#' @export
sim_bds_tree <- function(params_df, q_matrix, x0, max_taxa = 100, max_t = 100,
                         include_extinct = FALSE) {
  # Format similar to matrix generated in original diversitree fn
  k <- nrow(params_df)
  pars <- cbind(
    data.matrix(params_df[3:5]),
    matrix(t(q_matrix)[col(q_matrix) != row(q_matrix)], k, byrow = TRUE)
  )

  phy <- NULL
  while (is.null(phy)) {
    info <- make.tree.bisse_modified(pars, k, max_taxa, max_t, x0)
    phy <- me.to.ape.bisse(info[-1, ], info$state[1])
    if (!is.null(phy) && !include_extinct) {
      phy <- diversitree::prune(phy)
    }
  }
  return(phy)
}
