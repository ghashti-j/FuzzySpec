compute.SNN <- function(similarity = NULL, r = NULL){
  if(is.null(similarity)) stop("Argument 'similarity' required.")
  if(nrow(similarity) != ncol(similarity)) stop("Argument 'similarity' must be an n x n matrix.")
  if (is.null(r)) stop("Argument 'r' missing with no default")
  n <- nrow(similarity)
  SNN.S <- matrix(0, n, n)
  for(i in 1:(n-1)) {
    neighbours.i <- order(similarity[i,], decreasing = T)[-1][1:r]
    for(j in (i+1):n) {
      neighbours.j <- order(similarity[j,], decreasing = T)[-1][1:r]
      shared <- length(intersect(neighbours.i, neighbours.j))
      min.k <- r
      SNN.S[i,j] <- SNN.S[j,i] <- shared/min.k
    }
  }
  diag(SNN.S) <- 1
  return(SNN.S)
}
