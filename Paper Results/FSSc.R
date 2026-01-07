######## the following script allows us to the run the FSSc algorithm
######## C. Borgelt

# Borgelt, C. (2009). Fuzzy Subspace Clustering. 
# In: Fink, A., Lausen, B., Seidel, W., Ultsch, A. (eds) 
# Advances in Data Analysis, Data Handling and Business Intelligence. 
# Studies in Classification, Data Analysis, and Knowledge Organization. 
# Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-01044-6_8


FSC <- function(data, K, label.old = NULL, tao = 2, sigm = 0.01, m = 2, RS = 5, maxit = 500) {
  data.num <- nrow(data)
  data.dim <- ncol(data)
  best.result <- NULL
  best.obj <- Inf
  for (rs in 1:RS) {
    eps <- 1e-6  
    fitness <- numeric(maxit)
    U <- matrix(0, data.num, K) 
    w.temp <- matrix(0, K, data.dim)
    distants <- matrix(0, data.num, K)
    para.miu <- matrix(0, K, data.dim)
    para.w <- matrix(1/data.dim, K, data.dim) 
    if (rs == 1) {
      # initialize w/ k-means++
      centers.idx <- sample(1:data.num, K)
      for (i in 1:data.num) {
        dists <- numeric(K)
        for (k in 1:K) {
          dists[k] <- sum((data[i,]-data[centers.idx[k],])^2)
        }
        # dist to fuzz
        for (k in 1:K) {
          if (min(dists) == 0) {
            U[i,] <- 0
            U[i, which.min(dists)] <- 1
          } else {
            U[i,k] <- 1/sum((dists[k]/dists)^(2/(m-1)))
          }
        }
      }
    } else {
      # subsequent starts
      U <- matrix(runif(data.num*K), data.num, K)
      U <- U/rowSums(U)  
    }
    
    converged <- FALSE
    iter <- 0
    for (t in 1:maxit) {
      iter <- t
      # centers update
      U.m <- U^m  
      miu.up <- t(U.m) %*% data.matrix(data)
      U.sum <- colSums(U.m)
      U.sum[U.sum < eps] <- eps  
      para.miu <- miu.up/(U.sum %o% rep(1, data.dim))
      # weight update
      for (k in 1:K) {
        diff.squared <- sweep(data, 2, para.miu[k,], "-")^2
        weighted.diff <- sweep(diff.squared, 1, U.m[,k], "*")
        w.temp[k,] <- colSums(weighted.diff)+sigm
      }
      # features update
      if (tao > 1) {
        w.up <- w.temp^(-1/(tao-1))
        w.sum <- rowSums(w.up)
        w.sum[w.sum == 0] <- 1 
        para.w <- w.up/(w.sum %o% rep(1, data.dim))
      }
      for (k in 1:K) {
        diff <- sweep(data, 2, para.miu[k,], "-")^2
        weighted.diff <- sweep(diff, 2, para.w[k,]^tao, "*")
        distants[,k] <- rowSums(weighted.diff)
      }
      # fuzzy matrix update
      U.new <- matrix(0, data.num, K)
      for (i in 1:data.num) {
        zero.dist <- which(distants[i,] < eps)
        if (length(zero.dist) > 0) {
          U.new[i,] <- 0
          U.new[i, zero.dist] <- 1/length(zero.dist)
        } else {
          for (k in 1:K) {
            U.new[i,k] <- 1/sum((distants[i,k]/distants[i,])^(2/(m-1)))
          }
        }
      }
      fitness[t] <- sum(U.new^m*distants)+sigm*sum(para.w^tao)
      # convergence check
      if (t > 1) {
        if (abs(fitness[t]-fitness[t-1]) < eps || 
            max(abs(U.new-U)) < eps) {
          converged <- TRUE
          U <- U.new
          break
        }
      }
      U <- U.new
    }
    label <- apply(U, 1, which.max)
    # update all
    if (fitness[iter] < best.obj) {
      best.obj <- fitness[iter]
      best.result <- list(
        U = U,             
        H = para.miu,     
        clus = cbind(label, U),  
        iter = iter,
        func.val = fitness[iter],
        W = para.w,       
        converged = converged,
        m = m             
      )
    }
  }
  class(best.result) <- c("fclust", "list")
  return(best.result)
}