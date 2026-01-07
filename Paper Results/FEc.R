######## the following script allows us to the run the FEc algorithm
######## D. Tran; M. Wagner


# D. Tran and M. Wagner, "Fuzzy entropy clustering," 
# Ninth IEEE International Conference on Fuzzy Systems. FUZZ- IEEE 2000 (Cat. No.00CH37063), 
# San Antonio, TX, USA, 2000, pp. 152-157 vol.1, doi: 10.1109/FUZZY.2000.838650.

FEC <- function(data, K, label.old = NULL, gamma = 1, m = 2, RS = 5, maxit = 500) {
  data.num <- nrow(data)
  data.dim <- ncol(data)
  best.result <- NULL
  best.obj <- Inf
  for (rs in 1:RS) {
    eps <- 1e-6 
    esp <- 1e-6 
    fitness <- numeric(maxit)
    distant <- matrix(0, data.num, K)
    responsivity <- matrix(0, data.num, K)
    para.miu <- matrix(0, K, data.dim)
    
    # center init
    if (is.null(label.old) || rs > 1) {
      if (rs == 1) {
        centers.idx <- numeric(K)
        centers.idx[1] <- sample(1:data.num, 1)
        for (k in 2:K) {
          dists <- numeric(data.num)
          for (i in 1:data.num) {
            min.dist <- Inf
            for (j in 1:(k-1)) {
              d <- sum((data[i,]-data[centers.idx[j],])^2)
              if (d < min.dist) min.dist <- d
            }
            dists[i] <- min.dist
          }
          centers.idx[k] <- sample(1:data.num, 1, prob = dists)
        }
        for (k in 1:K) {
          para.miu[k,] <- data.matrix(data[centers.idx[k],])
        }
      } else {
        for (d in 1:data.dim) {
          para.miu[,d] <- runif(K, min(data[,d]), max(data[,d]))
        }
      }
    } else {
      for (k in 1:K) {
        X.k <- data[label.old == k, , drop = FALSE]
        if (nrow(X.k) > 0) {
          para.miu[k,] <- colMeans(X.k)
        } else {
          para.miu[k,] <- data[sample(1:data.num, 1),]
        }
      }
    }
    converged <- FALSE
    iter <- 0
    for (t in 1:maxit) {
      iter <- t
      for (k in 1:K) {
        diff <- sweep(data, 2, para.miu[k,], "-")
        distant[,k] <- rowSums(diff^2)
      }
      # membership update
      R.up <- exp(-distant/gamma)
      R.sum <- rowSums(R.up)
      R.sum[R.sum == 0] <- eps  
      responsivity <- R.up/R.sum
      # centers update
      miu.up <- t(responsivity) %*% data.matrix(data)
      resp.sum <- colSums(responsivity)
      resp.sum[resp.sum == 0] <- eps 
      para.miu <- miu.up/(resp.sum %o% rep(1, data.dim))
      # objective fcn
      entropy.term <- responsivity*log(responsivity+eps)
      fitness[t] <- sum(responsivity*distant)+gamma*sum(entropy.term)
      # convergence check
      if (t > 1) {
        if (abs(fitness[t]-fitness[t-1]) < esp) {
          converged <- TRUE
          break
        }
      }
    }
    label <- apply(responsivity, 1, which.max)
    entropy <- -sum(responsivity*log(responsivity+eps))/data.num
    # update all
    if (fitness[iter] < best.obj) {
      best.obj <- fitness[iter]
      best.result <- list(
        U = responsivity, 
        H = para.miu,     
        clus = cbind(label, responsivity),  
        iter = iter,
        func.val = fitness[iter],
        entropy = entropy,
        gamma = gamma,
        converged = converged
      )
    }
  }
  class(best.result) <- c("fclust", "list")
  return(best.result)
}