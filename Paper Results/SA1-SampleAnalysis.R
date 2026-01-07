
library(ggplot2)
library(mvtnorm)
library(np)
library(e1071)
library(mclust)
library(ppclust)
library(cluster)
library(dbscan)
library(kernlab)
library(lmtest)
library(abind)
library(fclust)

source("VWKFC.R")
source("FCMM.R")
source("FSSC.R")
source("FEc.R")
source("SNNs.R")
source("LDAs.R")
source("Functions.R")
source("SA1-Datasets.R")


# PARAMETERS
numIters <- 1000
numPts <- 100
m <- 2
saveResultsInterval <- 50


# METHOD NAMES
methodNames <- c(
  "Weu-id", "Wvw-id", "Weula-id", "Wvwla-id", 
  "Weu-sim", "Wvw-sim", 
  "Weu-simsnn", "Wvw-simsnn", "Weu-simsnns", "Wvw-simsnns", 
  "Weu-snn", "Wvw-snn", "Weu-snns", "Wvw-snns", 
  "Weula-snn", "Wvwla-snn", "Weula-snns", "Wvwla-snns",
  
  # COMPETING
  "LDAs", "SNNs", "FCMn", "GFCm", "MFCm", 
  "PFCm", "FSSc", "FEc", "GMM", "VWKFC", "FCMM"
)

numMethods <- length(methodNames)


# INITIALIZE RESULT MATRICES
accuracyMat <- matrix(NA, nrow = numMethods, ncol = numIters,
                          dimnames = list(methodNames, paste0("iter_", 1:numIters)))
fariMat <- matrix(NA, nrow = numMethods, ncol = numIters,
                      dimnames = list(methodNames, paste0("iter_", 1:numIters)))
runtimeMat <- matrix(NA, nrow = numMethods, ncol = numIters,
                         dimnames = list(methodNames, paste0("iter_", 1:numIters)))
entropyMat <- matrix(NA, nrow = numMethods, ncol = numIters,
                    dimnames = list(methodNames, paste0("iter_", 1:numIters)))
xbMat <- matrix(NA, nrow = numMethods, ncol = numIters,
                    dimnames = list(methodNames, paste0("iter_", 1:numIters)))

rholist <- list()
sigmalist <- list()
radiusMat <- matrix(0, 2, numIters)
snntimes <- numeric(numIters)
sigma.times <- numeric(numIters)
bw.times <- numeric(numIters)
cluster.sizes <- matrix(NA, nrow = 20, ncol = numIters)
rownames(radiusMat) <- c("eu", "vw")

# toggle off if desired
verbose = TRUE


# MAIN SIMULATION LOOP
####### NOTE ######
####### Here, any of the desired datasets from SA1 can be used instead by
####### simply switching the function call four lines below to the desired
####### data generating function
####### we have used the spiral dataset

#### options include 'genHypeData', 'genGaussianData', 'genWormData', 'genRingData'
#### 'genSpiralData, or 'genPieData'.
for (iter in 1:numIters) {
  if(verbose) print(sprintf("Iteration %d/%d", iter, numIters))
  if(verbose) if (iter %% 10 == 0) cat(sprintf(" [%.1f%%]", 100 * iter / numIters))
  
  data <- genSpiralData(nPts = numPts, seed = 44869410 + iter)
  X <- data$X
  uTrue <- data$uTrue
  yTrue <- data$yTrue
  currk <- data$k

  n <- nrow(X)
  p <- ncol(X)
  
  ################## NOTE ###################################
  # the code
  # bw <- npudensbw(X, bwmethod = "normal-reference", nmulti = 2)
  # is included to increase computation times
  # change to bwmethod = "cv.ls" to obtain results as in the paper
  tryCatch({
    t0 <- Sys.time()
    bw <- npudensbw(X, bwmethod = "normal-reference", nmulti = 2)
    bw.times[iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    h <- bw$bw
  }, error = function(e) {
    h <- apply(X, 2, sd) * n^(-1/(4+p))
  })
  
  w <- 1/h^2
  Xsc <- sweep(X, 2, sqrt(w), "*")
  D.kss <- as.matrix(dist(Xsc))
  D.kss.sq <- D.kss^2
  D.kss.sc <- D.kss/max(D.kss)
  S.kss <- exp(-D.kss.sq/2)
  S.kss.sc <- exp(-D.kss.sc)
  S.kss.sq <- S.kss^2
  D.eu <- as.matrix(dist(X))
  D.eu.sq <- D.eu^2
  S.eu <- 1/(1 + D.eu)
  
  t0 <- Sys.time()
  sigmaks <- compute.sigma(D.kss, "vw")
  sigma.times[iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  r.ks <- sigmaks$radius
  sigma.ks <- sigmaks$sigma
  
  t0 <- Sys.time()
  snn.ks <- compute.SNN(S.kss, r.ks)
  snntimes[iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  sigmaeu <- compute.sigma(D.eu, "eu")
  r.eu <- sigmaeu$radius
  sigma.eu <- sigmaeu$sigma
  snn.eu <- compute.SNN(S.eu, r.eu)
  
  rho.eu <- numeric(n)
  for(i in 1:n) {
    ss <- sort(S.eu[i,], decreasing = TRUE)
    rho.eu[i] <- ss[r.eu + 1]
  }
  
  rho.ks <- numeric(n)
  for(i in 1:n) {
    ss <- sort(S.kss[i,], decreasing = TRUE)
    rho.ks[i] <- ss[r.ks + 1]
  }
  
  rholist[[iter]] <- list(eu = rho.eu, ks = rho.ks)
  sigmalist[[iter]] <- list(eu = sigma.eu, ks = sigma.ks)
  radiusMat[,iter] <- c(r.eu, r.ks)
  
  spectralMethods <- methodNames[!methodNames %in% c("FCMn", "GFCm", "MFCm", 
                                                     "PFCm", "FSSc", "FEc", "GMM", "VWKFC", "FCMM")]
  
  adjacencyMatrix <- list()
  for(method in spectralMethods){
    tryCatch({
      adjacencyMatrix[[method]] <- make.similarity(data = X, method = method, n = n, 
                                         sigma.eu = sigma.eu, sigma.ks = sigma.ks, 
                                         rho.eu = rho.eu, rho.ks = rho.ks, 
                                         snn.eu = snn.eu, snn.ks = snn.ks, 
                                         S.kss = S.kss, S.eu = S.eu, 
                                         D.kss.sq = D.kss.sq, 
                                         D.eu.sq = D.eu.sq, D.eu = D.eu)
    }, error = function(e) { adjacencyMatrix[[method]] <- NULL })
  }
  
  # Run spectral methods
  for (method in names(adjacencyMatrix)) {
    if (!is.null(adjacencyMatrix[[method]])) {
      tryCatch({
        t0 <- Sys.time()
        res <- fuzzy.spectral.clustering(adjacencyMatrix[[method]], k = currk, m = m, method = "CM")
        runtimeMat[method, iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        if (!is.null(res)) {
          fariMat[method, iter] <- fari(uTrue, res$u)$fari
          accuracyMat[method, iter] <- clustering_accuracy(yTrue, res$cluster)
          entropyMat[method, iter] <- fclust::PE(res$u)
          xbMat[method, iter] <- fclust::XB(res$evecs, res$u, res$centers, m)
        }
      }, error = function(e) {
        message(sprintf("\nError in method %s: %s", method, e$message))
      })
    }
  }
  
  # FCMn
  tryCatch({
    t0 <- Sys.time()
    res <- fclust::FKM(X, k = currk, m = m, RS = 5, maxit = 500)
    runtimeMat["FCMn", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
      fariMat["FCMn", iter] <- fari(uTrue, res$U)$fari
      accuracyMat["FCMn", iter] <- clustering_accuracy(yTrue, res$clus[,1])
      entropyMat["FCMn", iter] <- fclust::PE(res$U)
      xbMat["FCMn", iter] <- fclust::XB(X, res$U, res$H, m)
    }
  }, error = function(e) {
    message(sprintf("\nError in FCMean: %s", e$message))
  })
  
  # GFCm - Gustafson-Kessel
  tryCatch({
    t0 <- Sys.time()
    res <- fclust::FKM.gk(X, k = currk, m = m, RS = 5, maxit = 500)
    runtimeMat["GFCm", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
      fariMat["GFCm", iter] <- fari(uTrue, res$U)$fari
      accuracyMat["GFCm", iter] <- clustering_accuracy(yTrue, res$clus[,1])
      entropyMat["GFCm", iter] <- fclust::PE(res$U)
      xbMat["GFCm", iter] <- fclust::XB(X, res$U, res$H, m)
    }
  }, error = function(e) {
    message(sprintf("\nError in GK: %s", e$message))
  })
  
  # MFCm - Fuzzy C-Medoids
  tryCatch({
    t0 <- Sys.time()
    res <- fclust::FKM.med(X, k = currk, m = m, RS = 5, maxit = 500)
    runtimeMat["MFCm", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
      fariMat["MFCm", iter] <- fari(uTrue, res$U)$fari
      accuracyMat["MFCm", iter] <- clustering_accuracy(yTrue, res$clus[,1])
      entropyMat["MFCm", iter] <- fclust::PE(res$U)
      xbMat["MFCm", iter] <- fclust::XB(X, res$U, res$H, m)
    }
  }, error = function(e) {
    message(sprintf("\nError in FCMed: %s", e$message))
  })
  
  # PFCm - Fuzzy C-Means with polynomial fuzzifier
  tryCatch({
    t0 <- Sys.time()
    res <- fclust::FKM.pf(X, k = currk, b = 0.5, RS = 5, maxit = 500)
    runtimeMat["PFCm", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
      fariMat["PFCm", iter] <- fari(uTrue, res$U)$fari
      accuracyMat["PFCm", iter] <- clustering_accuracy(yTrue, res$clus[,1])
      entropyMat["PFCm", iter] <- fclust::PE(res$U)
      xbMat["PFCm", iter] <- fclust::XB(X, res$U, res$H, m)
    }
  }, error = function(e) {
    message(sprintf("\nError in FCPoly: %s", e$message))
  })
  
  # FSSc fuzzy subspace clustering
  tryCatch({
    t0 <- Sys.time()
    res <- FSC(X, m = m, K = currk)
    runtimeMat["FSSc", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
      fariMat["FSSc", iter] <- fari(uTrue, res$U)$fari
      accuracyMat["FSSc", iter] <- clustering_accuracy(yTrue, res$clus[,1])
      entropyMat["FSSc", iter] <- fclust::PE(res$U)
      xbMat["FSSc", iter] <- fclust::XB(X, res$U, res$H, m)
    }
  }, error = function(e) {
    message(sprintf("\nError in FSC: %s", e$message))
  })
  
  # fuzzy entropy clustering
  tryCatch({
    t0 <- Sys.time()
    res <- FEC(X, m = m, K = currk)
    runtimeMat["FEc", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
      fariMat["FEc", iter] <- fari(uTrue, res$U)$fari
      accuracyMat["FEc", iter] <- clustering_accuracy(yTrue, res$clus[,1])
      entropyMat["FEc", iter] <- fclust::PE(res$U)
      xbMat["FEc", iter] <- fclust::XB(X, res$U, res$H, m)
    }
  }, error = function(e) {
    message(sprintf("\nError in FEC: %s", e$message))
  })
  
  # GMM
  tryCatch({
    t0 <- Sys.time()
    res <- mclust::Mclust(data = X, G = currk)
    runtimeMat["GMM", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
      fariMat["GMM", iter] <- fari(uTrue, res$z)$fari
      accuracyMat["GMM", iter] <- clustering_accuracy(yTrue, res$classification)
      entropyMat["GMM", iter] <- fclust::PE(res$z)
      xbMat["GMM", iter] <- fclust::XB(X, res$z, t(res$parameters$mean), 2)
    }
  }, error = function(e) {
    message(sprintf("\nError in GMM: %s", e$message))
  })
  
  # VWKFC
  tryCatch({
    t0 <- Sys.time()
    res <- vwkfc(X, C = currk, m = m, sigma = 1, gamma = 1, eps = 1e-6, iM = 100, verbose = FALSE)
    runtimeMat["VWKFC", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
    fariMat["VWKFC", iter] <- fari(uTrue, t(res$U))$fari
    accuracyMat["VWKFC", iter] <- clustering_accuracy(yTrue, apply(t(res$U), 1, which.max))
    entropyMat["VWKFC", iter] <- fclust::PE(t(res$U))
    xbMat["VWKFC", iter] <- fclust::XB(X, t(res$U), res$W, 2)
    }
  }, error = function(e) {
    message(sprintf("\nError in VWKFC: %s", e$message))
  })
  
  # FCMM
  tryCatch({
    t0 <- Sys.time()
    res <- fcmm(X, c = currk, q = 20, r = 2, alpha = 1,
                max_iter = 200, tol = 1e-6, init = "kmeans",
                verbose = FALSE)
    runtimeMat["FCMM", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(res)) {
    fariMat["FCMM", iter] <- fari(uTrue, res$L)$fari
    accuracyMat["FCMM", iter] <- clustering_accuracy(yTrue, apply(res$L, 1, which.max))
    entropyMat["FCMM", iter] <- fclust::PE(res$L)
    xbMat["FCMM", iter] <- fclust::XB(X, res$L, res$V, 2)
    }
  }, error = function(e) {
    message(sprintf("\nError in VWKFC: %s", e$message))
  })
  
  if (iter %% saveResultsInterval == 0) {
    resultsList <- list(
      accuracy = accuracyMat[, 1:iter],
      fari = fariMat[, 1:iter],
      runtime = runtimeMat[, 1:iter],
      PE = entropyMat[, 1:iter],
      XB = xbMat[, 1:iter],
      rho = rholist[1:iter],
      sigma = sigmalist[1:iter],
      radii = radiusMat[,1:iter],
      snntimes = snntimes[1:iter],
      sigma.times = sigma.times[1:iter],
      bw.times = bw.times[1:iter],
      clussize = cluster.sizes[,1:iter],
      parameters = list(
        numIters_completed = iter,
        numPts = numPts,
        k = currk,
        m = m,
        methodNames = methodNames
      )
    )
    
    saveRDS(resultsList, file = sprintf("spirals_%d.rds", iter))
    if(verbose) cat(sprintf("\n[Saved checkpoint at iteration %d]\n", iter))
  }
}


finResults <- list(
  accuracy = accuracyMat,
  fari = fariMat,
  runtime = runtimeMat,
  PE = entropyMat,
  XB = xbMat,
  rho = rholist,
  sigma = sigmalist,
  radii = radiusMat,
  snntimes = snntimes,
  sigma.times = sigma.times,
  clussize = cluster.sizes,
  bw.times = bw.times,
  parameters = list(
    numPts = numPts,
    k = currk,
    m = m,
    methodNames = methodNames
  )
)

# saveRDS(finResults, file = "finalSpiralResults.rds")
