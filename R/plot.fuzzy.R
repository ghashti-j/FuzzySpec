if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("x1", "x2", "cluster", "uncertainty"))
}

plot.fuzzy <- function(data, plotFuzzy = TRUE, colorCluster = TRUE) {
  X <- data$X; U <- data$U; y <- data$y; k <- data$k
  plotDF <- data.frame(
    x1 = X[,1],
    x2 = X[,2],
    cluster = factor(y),
    uncertainty = 1 - apply(U, 1, max)
  )

  p <- ggplot2::ggplot(plotDF, aes(x = x1, y = x2))

  if (colorCluster && plotFuzzy) p <- p + ggplot2::geom_point(aes(color = cluster, size = uncertainty), alpha = 0.7)
  else if (colorCluster && !plotFuzzy) p <- p + ggplot2::geom_point(aes(color = cluster), size = 2, alpha = 0.7)
  else if (!colorCluster && plotFuzzy) p <- p + ggplot2::geom_point(aes(size = uncertainty), color = "black", alpha = 0.7)
  else p <- p + ggplot2::geom_point(color = "black", size = 2, alpha = 0.7)

  if (plotFuzzy) {
    p <- p + ggplot2::scale_size_continuous(
      range = c(0.5, 5),
      name = "Uncertainty",
      breaks = c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.99),
      labels =  c("0.01", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90", "0.99")
    )
  }

  if (colorCluster) if (k <= 8) p <- p + ggplot2::scale_color_brewer(palette = "Set1", name = "Cluster") else  p <- p + scale_color_viridis_d(name = "Cluster")
  p <- p +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = element_rect(NA, "black", 1),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = if(colorCluster || plotFuzzy) "bottom" else "none"
    ) +
    ggplot2::labs(
      x = expression(X[1]),
      y = expression(X[2])
    ) +
    ggplot2::coord_equal()
  return(p)
}

