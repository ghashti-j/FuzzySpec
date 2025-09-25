compute.sigma <- function(distance, r = NULL) {
  if(is.null(distance)) stop("Argument 'distance' required.")
  if(nrow(distance) != ncol(distance)) stop("Argument 'distance' must be an n x n matrix.")
  if (is.null(r)) {
    warning("Argument 'r' not provided, determining 'r' via find.radius()")
    r <- find.radius(distance)
  }
  n <- nrow(distance)
  sigma <- numeric(n)
  for (i in 1:n) sigma[i] <- sort(distance[i, ])[r + 1]
  return(list(sigma = sigma, radius = r))
}
