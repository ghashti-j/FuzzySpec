#' Plot method for fuzzy clustering results
#'
#' Creates a visualization of fuzzy clustering results with optional
#' uncertainty representation and cluster coloring.
#'
#' @param data A fuzzy clustering object or list containing X (data), U (membership matrix),
#'   y (cluster assignments), and k (number of clusters)
#' @param plotFuzzy Logical; if TRUE, point sizes represent membership uncertainty
#' @param colorCluster Logical; if TRUE, points are colored by cluster assignment
#' @param ... Additional arguments (unused)
#'
#' @return A ggplot2 object
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous scale_color_brewer
#' @importFrom ggplot2 theme_classic theme element_rect element_text labs coord_equal
#' @importFrom viridisLite viridis
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- gen.fuzzy(n = 300, dataset = "spirals", noise = 0.18)
#' plot.fuzzy(d)
#'
plot.fuzzy <- function(x, plotFuzzy = TRUE, colorCluster = TRUE, ...) {
  X <- x$X
  U <- x$U
  y <- x$y
  k <- x$k
  plotDF <- data.frame(
    x1 = X[, 1],
    x2 = X[, 2],
    cluster = factor(y),
    uncertainty = 1 - apply(U, 1, max)
  )
  p <- ggplot2::ggplot(plotDF, ggplot2::aes(x = x1, y = x2))
  if (colorCluster && plotFuzzy) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = cluster, size = uncertainty),
      alpha = 0.7
    )
  } else if (colorCluster && !plotFuzzy) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = cluster),
      size = 2,
      alpha = 0.7
    )
  } else if (!colorCluster && plotFuzzy) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(size = uncertainty),
      color = "black",
      alpha = 0.7
    )
  } else {
    p <- p + ggplot2::geom_point(
      color = "black",
      size = 2,
      alpha = 0.7
    )
  }
  if (plotFuzzy) {
    p <- p + ggplot2::scale_size_continuous(
      range = c(0.5, 5),
      name = "Uncertainty",
      breaks = c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.99),
      labels = c("0.01", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90", "0.99")
    )
  }
  if (colorCluster) {
    if (k <= 8) {
      p <- p + ggplot2::scale_color_brewer(palette = "Set1", name = "Cluster")
    } else {
      # Use viridisLite for the color scale instead of scales
      p <- p + ggplot2::scale_color_manual(
        values = viridisLite::viridis(k),
        name = "Cluster"
      )
    }
  }
  p <- p +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12),
      legend.position = if (colorCluster || plotFuzzy) "bottom" else "none"
    ) +
    ggplot2::labs(
      x = expression(X[1]),
      y = expression(X[2])
    ) +
    ggplot2::coord_equal()

  return(p)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("x1", "x2", "cluster", "uncertainty"))
}
