plot_lr_wiring <- function(ligand_exprs, receptor_exprs, cell_labels, thresh, title_txt = NULL) {
  # coerce & sanitize
  ligand_exprs   <- suppressWarnings(as.numeric(ligand_exprs))
  receptor_exprs <- suppressWarnings(as.numeric(receptor_exprs))
  ligand_exprs[!is.finite(ligand_exprs)]     <- 0
  receptor_exprs[!is.finite(receptor_exprs)] <- 0
  
  n_cell   <- length(cell_labels)
  ligand_exprs[ligand_exprs < thresh | is.na(ligand_exprs)]     <- 0
  receptor_exprs[receptor_exprs < thresh | is.na(receptor_exprs)] <- 0
  
  # NxN product; replace non-finite with 0
  norm_mat <- outer(ligand_exprs, receptor_exprs, `*`)
  norm_mat[!is.finite(norm_mat)] <- 0
  
  total <- sum(norm_mat, na.rm = TRUE)
  if (total > 0) {
    final_mat <- norm_mat / total
    final_mat[final_mat < 0.1] <- 0
  } else {
    final_mat <- norm_mat
  }
  
  dimnames(final_mat) <- list(cell_labels, cell_labels)
  
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  net <- igraph::graph_from_adjacency_matrix(final_mat, mode = "directed", weighted = TRUE)
  
  # drop zero/NA/neg edges
  if (igraph::ecount(net)) {
    w <- igraph::E(net)$weight
    kill <- which(!is.finite(w) | w <= 0)
    if (length(kill)) net <- igraph::delete_edges(net, kill)
  }
  
  # scale widths
  if (igraph::ecount(net)) {
    w <- igraph::E(net)$weight
    igraph::E(net)$width <- 4 * (w / max(w))
  }
  
  co <- igraph::layout_in_circle(net)
  
  op <- graphics::par(mar = c(0.3, 0.3, 2.5, 0.3)); on.exit(graphics::par(op), add = TRUE)
  
  if (igraph::ecount(net) > 0) {
    graphics::plot(net,
                   edge.label         = round(igraph::E(net)$weight, 2),
                   edge.curved        = TRUE,
                   vertex.label.font  = 2,
                   edge.arrow.size    = 0.6,
                   layout             = co,
                   margin             = c(0.4, 0.4, 0.4, 0.4),
                   vertex.size        = 25,
                   vertex.label.cex   = 2,
                   vertex.label.dist  = 4,
                   vertex.label.degree= -pi/2,
                   edge.color         = "black",
                   vertex.color       = "white",
                   edge.label.color   = "brown",
                   edge.label.cex     = 1.5,
                   edge.label.font    = 2
    )
  } else {
    graphics::plot(net,
                   lty = 0, vertex.label.font = 2, layout = co,
                   margin = c(0.4, 0.4, 0.4, 0.4), vertex.size = 25,
                   vertex.label.cex = 2, vertex.label.dist = 4, vertex.label.degree = -pi/2,
                   vertex.color = "white"
    )
  }
  
  if (!is.null(title_txt) && nzchar(title_txt)) {
    graphics::title(main = title_txt, cex.main = 1.4, line = 0.5)
  }
  
  invisible(net)
}
