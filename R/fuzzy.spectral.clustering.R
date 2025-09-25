fuzzy.spectral.clustering <- function(W = NULL,
                                      k = NULL,
                                      m = NULL,
                                      method = "CM",
                                      nstart = 10,
                                      max.iter = 1000) {
  if(is.null(W)) stop("Adjacency (similarity) matrix 'W' required.")
  if(is.null(k)) stop("Number of clusters 'k' required.")
  if(is.null(m)){
    warning("Argument 'm', missing, defaulting to m=2")
    m <- 2
  }
  diag(W) <- 0
  d <- 1/sqrt(rowSums(W))
  l <- d*W%*%diag(d)
  XI <- eigen(l)
  xi <- XI$vectors[,1:k]
  yi <- xi/sqrt(rowSums(xi^2))
  if(method == "CM"){
    km <- fclust::FKM(X = yi, k = k, m = m, RS = nstart, maxit = max.iter)
    return(list(cluster = km$clus[,1],
                u = km$U,
                evecs   = yi,
                centers = km$H))
  }
  if(method == "GMM"){
    km <- mclust::Mclust(yi, G=k)
    return(list(cluster = km$classification,
                u = km$z,
                evecs   = yi,
                centers = km$parameters$mean))
  }
}
