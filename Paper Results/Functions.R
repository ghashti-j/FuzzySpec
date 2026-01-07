pkgs <- c(
  "ppclust","mclust","MASS","mvtnorm","np","ggplot2","dplyr","tidyr",
  "Rfast","Spectrum","pracma","Matrix","grid","e1071","abind","lmtest",
  "ggforce","cluster","dbscan","kernlab","fclust"
)

missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing)) install.packages(missing)
invisible(lapply(pkgs, library, character.only = TRUE))

determine.u <- function(y_true) {
  ydisc <- factor(y_true)     
  n <- length(ydisc)
  k <- nlevels(ydisc)
  U <- matrix(0L, nrow = n, ncol = k,
              dimnames = list(NULL, levels(ydisc)))
  U[cbind(seq_len(n), as.integer(ydisc))] <- 1L
  U
}

fari <-
  function(a, b){
    #A, B bonding matrices
    n <- nrow(a)
    A <- a %*% t(a)
    B <- b %*% t(b)
    j <- matrix(1, n, n)
    Na <- (sum(A*j)/sum(A*A))*A
    Nb <- (sum(B*j)/sum(B*B))*B
    ri <- (sum(Na*Nb)+sum((j-Na)*(j-Nb))-n)/(2*choose(n,2))
    M <- j/n
    R <- diag(n)-M
    Eri <- ((2*sum(A*j)*sum(B*j)/(sum(A*A)*sum(B*B)))*(sum(M*A)*sum(M*B)+((1/(n-1))*sum(R*A)*sum(R*B)))-(sum(A*j)^2)/sum(A*A)-(sum(B*j)^2)/sum(B*B)+n^2-n)/(2*choose(n,2))
    ari <- (ri-Eri)/(1-Eri)
    store <- list(fri=ri, fari=ari)
    return(store)
  }


clustering_accuracy <- function(trueLabs, predLabs) {
  confMat <- table(trueLabs, predLabs)
  nTrue <- nrow(confMat)
  nPred <- ncol(confMat)
  correct <- 0
  if (nTrue <= nPred) {
    for (i in 1:nTrue) {
      correct <- correct+max(confMat[i,])
    }
  } else {
    for (j in 1:nPred) {
      correct <- correct+max(confMat[,j])
    }
  }
  accuracy <- correct/length(trueLabs)
  return(accuracy)
}

fuzzy.spectral.clustering <- function(W,k,m,method = NULL) {
  if(is.null(method)) stop("Input GMM, CM, or BK")
  diag(W) <- 0
  d <- 1/sqrt(rowSums(W))
  l <- d*W%*%diag(d)
  XI <- eigen(l)
  eg <- XI$values[k]-XI$values[k+1]
  xi <- XI$vectors[,1:k]
  yi <- xi/sqrt(rowSums(xi^2))
  if(method == "CM"){
    km <- e1071::cmeans(yi, centers = k,iter.max = 1000, m = m)
    return(list(cluster = km$cluster,
                u = km$membership,
                evecs   = yi,
                centers = km$centers,
                eigengap = eg))
  }
  if(method == "GMM"){
    km <- mclust::Mclust(yi, G=k)
    return(list(cluster = km$classification,
                u = km$z,
                evecs   = yi,
                centers = km$parameters$mean,
                eigengap = eg))
  }
  if(method == "BK"){
    km <- bootk(yi, groups = k, itermax = 1000)
    return(list(cluster = km$hardoob,
                u = km$fuzzoob,
                evecs   = yi,
                eigengap = eg))
  }
}

kNN.dist <- function(D, k) {
  n <- nrow(D)
  id <- t(apply(D, 1, function(d) order(d)[-1][1:k]))
  list(id = id)
}


find.radius <- function(D_w_dist) {
  n <- nrow(D_w_dist)
  r <- 1
  nn <- ceiling(sqrt(n))
  while(TRUE) {
    if(r > nn) nn <- r+10
    dist.obj <- kNN.dist(D_w_dist, nn)
    NaN_Num <- tabulate(dist.obj$id[, 1:r], nbins = n)
    NaN_Num_0 <- length(NaN_Num[NaN_Num == 0])
    if(r == 1) Nan_Num_0_Upd <- NaN_Num_0
    if(r > 1 & Nan_Num_0_Upd == NaN_Num_0) break
    Nan_Num_0_Upd <- length(NaN_Num[NaN_Num == 0])
    r <- r+1
  }
  return(r)
}

compute.sigma <- function(distance, method = NULL){
  r <- find.radius(distance)
  n <- nrow(distance)
  if(is.null(method)) stop ("Euclidean or DKSS")
  if(method == "vw") {
    sigma.kss <- numeric(n)
    for(i in 1:n) {
      sigma.kss[i] <- sort(distance[i,])[r+1]
    }
    return(list(sigma = sigma.kss, radius = r))
  }
  if(method == "eu"){
    sigma.eu <- numeric(n)
    for(i in 1:n) {
      sigma.eu[i] <- sort(distance[i,])[r+1]
    }
    return(list(sigma = sigma.eu, radius = r))
  }
}

compute.SNN <- function(S.kss, r){
  n <- nrow(S.kss)
  SNN.S <- matrix(0, n, n)
  for(i in 1:(n-1)) {
    neighbors_i <- order(S.kss[i,], decreasing = T)[-1][1:r]
    for(j in (i+1):n) {
      neighbors_j <- order(S.kss[j,], decreasing = T)[-1][1:r]
      shared <- length(intersect(neighbors_i, neighbors_j))
      min_k <- r
      SNN.S[i,j] <- SNN.S[j,i] <- shared/min_k
    }
  }
  diag(SNN.S) <- 1
  return(SNN.S)
}

make.similarity <- function(data, method = NULL, n, 
                            sigma.eu, sigma.ks, 
                            rho.eu, rho.ks, 
                            snn.eu, snn.ks, 
                            S.kss, S.eu, 
                            D.kss.sq, D.eu.sq, D.eu){
  if(method == "Weu-id")  return(W = S.eu)
  if(method == "Wvw-id")  return(W = S.kss)
  if(method == "Weula-id")  return(W = exp(-D.eu.sq/outer(sigma.eu, sigma.eu)))
  if(method == "Wvwla-id")  return(W = exp(-D.kss.sq/outer(sigma.ks, sigma.ks)))
  if(method == "Weu-sim")  return(W = S.eu*(1+sqrt(outer(rho.eu, rho.eu)))/2)
  if(method == "Wvw-sim")  return(W = S.kss*(1+sqrt(outer(rho.ks, rho.ks)))/2)
  if(method == "Weu-simsnns")  return(W = S.eu*ifelse(snn.eu > 0, (1+snn.eu)/2, 0)*(1+sqrt(outer(rho.eu, rho.eu)))/2)
  if(method == "Wvw-simsnns")  return(W = S.kss*ifelse(snn.ks > 0, (1+snn.ks)/2, 0)*(1+sqrt(outer(rho.ks, rho.ks)))/2)
  if(method == "Weu-snns")  return(W = S.eu*ifelse(snn.eu > 0, (1+snn.eu)/2, 0))
  if(method == "Wvw-snns")  return(W = S.kss*ifelse(snn.ks > 0, (1+snn.ks)/2, 0))
  if(method == "Weu-snn")  return(W = S.eu*(1+snn.eu)/2)
  if(method == "Wvw-snn")  return(W = S.kss*(1+snn.ks)/2)
  if(method == "Weu-simsnn")  return(W = S.eu*((1+snn.eu)/2)*((1+sqrt(outer(rho.eu, rho.eu)))/2))
  if(method == "Wvw-simsnn")  return(W = S.kss*((1+snn.ks)/2)*((1+sqrt(outer(rho.ks, rho.ks)))/2))
  if(method == "Weula-snns")  return(W = exp(-(D.eu.sq/outer(sigma.eu, sigma.eu)))*ifelse(snn.eu > 0, (snn.eu+1)/2, 0))
  if(method == "Wvwla-snns")  return(W = exp(-(D.kss.sq/outer(sigma.ks, sigma.ks)))*ifelse(snn.ks > 0, (snn.ks+1)/2, 0))
  if(method == "Weula-snn")  return(W = exp(-(D.eu.sq/outer(sigma.eu, sigma.eu)))*(snn.eu+1)/2)
  if(method == "Wvwla-snn")  return(W = exp(-(D.kss.sq/outer(sigma.ks, sigma.ks)))*(snn.ks+1)/2)
  if(method == "snneu") return(W = snn.eu)
  if(method == "snnvw") return(W = snn.ks)
  
  if (method == "LDAs") return(ldasm(data))
  if (method == "SNNs") return(scsnn(data, k = find.radius(D.eu)))
  if (method == "SCLSM") return(sclsm(data, k = find.radius(D.eu)))
}

