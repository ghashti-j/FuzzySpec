######## the following script allows us to the run the SNNs algorithm
######## Xiucai Ye; T. Sakurai

# Xiucai Ye and T. Sakurai, 
# "Spectral clustering using robust similarity measure based on closeness of shared Nearest Neighbors," 
# 2015 International Joint Conference on Neural Networks (IJCNN), Killarney, Ireland, 
# 2015, pp. 1-8, doi: 10.1109/IJCNN.2015.7280495.

scsnn <- function(X, k = 10) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (k >= n) stop("k must be < n")
  D  <- as.matrix(dist(X))                   
  knn.idx <- t(apply(D, 1, function(row) {   
    order(row, decreasing = FALSE)[-1L][1:k] 
  }))                                         
  pos.mat <- matrix(0L, n, n)
  for (i in 1:n) pos.mat[i, knn.idx[i, ]] <- seq_len(k)
  w <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      shared <- intersect(knn.idx[i, ], knn.idx[j, ])
      if (length(shared)) {
        for (r in shared) {
          pi <- pos.mat[i, r]         
          pj <- pos.mat[j, r]         
          w[i, j] <- w[i, j]+(k-pi+1)*(k-pj+1)
        }
      }
      if (pos.mat[i, j] > 0L && pos.mat[j, i] > 0L) {
        pi <- pos.mat[i, j]           
        pj <- pos.mat[j, i]           
        w[i, j] <- w[i, j]+(k-pi+1)*(k-pj+1)
      }
      w[j, i] <- w[i, j]            
    }
  }
  max.w <- max(w)      
  S <- w/max.w
  diag(S) <- 1
  return(S)
}
