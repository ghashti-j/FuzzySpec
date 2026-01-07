# SA3 generation, visualization and results from paper

# required functions
source("VWKFC.R")
source("FCMM.R")
source("FSSC.R")
source("FEc.R")
source("SNNs.R")
source("LDAs.R")
source("Functions.R")

# data generation function 
genNoiseVarData <- function(nNoiseDims = 0, nPointsPerCluster = 250, seed = NULL,
                            clusterSeparation = 4,clusterSigma = 0.8, 
                            noiseRange = c(-5, 5)) {
  if (!is.null(seed)) set.seed(seed)
  center1 <- c(-clusterSeparation/2, 0)
  center2 <- c(clusterSeparation/2, 0)

  cluster1Data <- mvtnorm::rmvnorm(nPointsPerCluster,
                                   mean = center1,
                                   sigma = diag(2) * clusterSigma^2)
  cluster2Data <- mvtnorm::rmvnorm(nPointsPerCluster,
                                   mean = center2,
                                   sigma = diag(2) * clusterSigma^2)
  XBase <- rbind(cluster1Data, cluster2Data)
  nTotal <- nrow(XBase)
  if (nNoiseDims > 0) {
    noiseDims <- matrix(runif(nTotal * nNoiseDims,
                              min = noiseRange[1],
                              max = noiseRange[2]),
                        nrow = nTotal,
                        ncol = nNoiseDims)
    X <- cbind(XBase, noiseDims)
  } else X <- XBase

  yTrue <- c(rep(1, nPointsPerCluster),
             rep(2, nPointsPerCluster))
  uFuzzy <- matrix(0, nTotal, 2)
  
  for (i in 1:nTotal) {
    d1 <- mvtnorm::dmvnorm(XBase[i,], mean = center1, sigma = diag(2) * clusterSigma^2)
    d2 <- mvtnorm::dmvnorm(XBase[i,], mean = center2, sigma = diag(2) * clusterSigma^2)
    totalD <- d1 + d2
    if (totalD > 0) {
      uFuzzy[i, 1] <- d1 / totalD
      uFuzzy[i, 2] <- d2 / totalD
    } else {
      dist1 <- sqrt(sum((XBase[i,] - center1)^2))
      dist2 <- sqrt(sum((XBase[i,] - center2)^2))
      if (dist1 < dist2) uFuzzy[i, ] <- c(0.9, 0.1) else uFuzzy[i, ] <- c(0.1, 0.9)
    }
  }
  
  # Bhattacharyya coefficient overlap in 2D space
  d <- sqrt(sum((center1 - center2)^2))
  overlap <- exp(-0.25 * d^2 / (2 * clusterSigma^2))
  centers <- rbind(colMeans(X[which(yTrue == 1),]), 
                   colMeans(X[which(yTrue == 2),]))
  
  return(list(
    X = X,
    XBase = XBase, 
    uTrue = uFuzzy,
    yTrue = yTrue,
    k = 2,
    nNoiseDims = nNoiseDims,
    totalDims = ncol(X),
    centers = centers,
    clusterSigma = clusterSigma,
    clusterOverlap = overlap,
    noiseRange = noiseRange
  ))
}

viz2DClust <- function(data, title = NULL, dimA = 1, dimB = 2) {
  df <- data.frame(
    x = data$X[, dimA],
    y = data$X[, dimB],
    cluster = as.factor(data$yTrue)
  )
  df$uMax <- apply(data$uTrue, 1, max)
  centersDF <- data.frame(
    x = data$centers[, dimA],
    y = data$centers[, dimB],
    cluster = as.factor(1:2)
  )
  colours <- c("#1F78B4", "#33A02C")
  
  p <- ggplot(df, aes(x = x, y = y, color = cluster)) +
    geom_point(aes(size = uMax), alpha = 0.8) +
    geom_point(data = centersDF, aes(x = x, y = y),
               color = "black", size = 4, shape = 4) +
    scale_color_manual(values = colours, guide = "none") +
    scale_size_continuous(range = c(5, 1.5)) +
    coord_cartesian(xlim = c(-5.1, 5.1), ylim = c(-5.1, 5.1)) +
    labs(x = bquote(X[.(dimA)]),
         y = bquote(X[.(dimB)]),
         size = expression(U[i*","*C[j]])) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      panel.background = element_rect(NA, "black", 1),
      legend.position = "bottom",
      legend.key = element_blank()
    )
  if (dimA == 1 && dimB == 2) {
    for (i in 1:2) {
      circleData <- data.frame(
        x = data$centers[i, 1],
        y = data$centers[i, 2],
        r = 2 * data$clusterSigma
      )
      p <- p + ggforce::geom_circle(
        data = circleData,
        aes(x0 = x, y0 = y, r = r),
        color = colours[i],
        fill = NA, linetype = "dashed", inherit.aes = FALSE
      )
    }
  }
  p
}

# PARAMETERS
numItersPerNoise <- 1000
maxNoiseDims <- 20
minNoiseDims <- 0
nPointsPerCluster <- 50
m <- 2
saveInterval <- 20
baseSeed <- 44869410
clusterSeparation <- 3
clusterSigma <- 0.8
noiseRange <- c(-5, 5)

verbose <- TRUE

# METHOD NAMES
methodNames <- c(
  "Weu-id", "Wvw-id", "Weula-id", "Wvwla-id",
  "Weu-sim", "Wvw-sim",
  "Weu-simsnn", "Wvw-simsnn", "Weu-simsnns", "Wvw-simsnns",
  "Weu-snn", "Wvw-snn", "Weu-snns", "Wvw-snns",
  "Weula-snn", "Wvwla-snn", "Weula-snns", "Wvwla-snns",
  
  # COMPETING METHODS
  "LDAs", "SNNs", "FCMn", "GFCm", "MFCm",
  "PFCm", "FSSc", "FEc", "GMM", "VWKFC", "FCMM"
)
numMethods <- length(methodNames)

# INITIALIZE STORAGE
totalIterations <- numItersPerNoise * (maxNoiseDims - minNoiseDims + 1)
resultsByNoise <- list()
for (nNoise in minNoiseDims:maxNoiseDims) {
  resultsByNoise[[paste0("noise", nNoise)]] <- list(
    accuracy = matrix(NA, nrow = numMethods, ncol = numItersPerNoise,
                      dimnames = list(methodNames, paste0(1:numItersPerNoise))),
    fari = matrix(NA, nrow = numMethods, ncol = numItersPerNoise,
                  dimnames = list(methodNames, paste0(1:numItersPerNoise))),
    runtime = matrix(NA, nrow = numMethods, ncol = numItersPerNoise,
                     dimnames = list(methodNames, paste0(1:numItersPerNoise))),
    entropy = matrix(NA, nrow = numMethods, ncol = numItersPerNoise,
                     dimnames = list(methodNames, paste0(1:numItersPerNoise))),
    xb = matrix(NA, nrow = numMethods, ncol = numItersPerNoise,
                dimnames = list(methodNames, paste0(1:numItersPerNoise))),
    rholist = list(),
    sigmalist = list(),
    radiusMat = matrix(0, 2, numItersPerNoise),
    snntimes = numeric(numItersPerNoise),
    sigma.times = numeric(numItersPerNoise),
    bw.times = numeric(numItersPerNoise)
  )
}
rownames(resultsByNoise[[1]]$radiusMat) <- c("eu", "vw")

# MAIN EXPERIMENT LOOP
overallIter <- 0
for (currentNoise in minNoiseDims:maxNoiseDims) {
  if(verbose) cat(sprintf("\n  Start %d noise dimensions \n", currentNoise))
  for (iter in 1:numItersPerNoise) {
    overallIter <- overallIter + 1
    data <- genNoiseVarData(
      nNoiseDims = currentNoise,
      nPointsPerCluster = nPointsPerCluster,
      seed = 1000 * currentNoise + iter,
      clusterSeparation = clusterSeparation,
      clusterSigma = clusterSigma,
      noiseRange = noiseRange)
    
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
      resultsByNoise[[paste0("noise", currentNoise)]]$bw.times[iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
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
    
    # sigma values
    t0 <- Sys.time()
    sigmaks <- compute.sigma(D.kss, "vw")
    resultsByNoise[[paste0("noise", currentNoise)]]$sigma.times[iter] <- 
      as.numeric(difftime(Sys.time(), t0, units = "secs"))
    r.ks <- sigmaks$radius
    sigma.ks <- sigmaks$sigma
    
    t0 <- Sys.time()
    snn.ks <- compute.SNN(S.kss, r.ks)
    resultsByNoise[[paste0("noise", currentNoise)]]$snntimes[iter] <- 
      as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    sigmaeu <- compute.sigma(D.eu, "eu")
    r.eu <- sigmaeu$radius
    sigma.eu <- sigmaeu$sigma
    snn.eu <- compute.SNN(S.eu, r.eu)
    
    # rho values
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
    
    resultsByNoise[[paste0("noise", currentNoise)]]$rholist[[iter]] <- list(eu = rho.eu, ks = rho.ks)
    resultsByNoise[[paste0("noise", currentNoise)]]$sigmalist[[iter]] <- list(eu = sigma.eu, ks = sigma.ks)
    resultsByNoise[[paste0("noise", currentNoise)]]$radiusMat[,iter] <- c(r.eu, r.ks)
    
    # spectral methods
    spectralMethods <- methodNames[!methodNames %in% c("FCMn", "GFCm", "MFCm",
                                                       "PFCm", "FSSc", "FEc", "GMM", "VWKFC", "FCMM")]
    
    adjacencyMatrix <- list()
    for(method in spectralMethods) {
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
    
    for (method in names(adjacencyMatrix)) {
      if (!is.null(adjacencyMatrix[[method]])) {
        tryCatch({
          t0 <- Sys.time()
          res <- fuzzy.spectral.clustering(adjacencyMatrix[[method]], k = currk, m = m, method = "CM")
          resultsByNoise[[paste0("noise", currentNoise)]]$runtime[method, iter] <- 
            as.numeric(difftime(Sys.time(), t0, units = "secs"))
          if (!is.null(res)) {
            resultsByNoise[[paste0("noise", currentNoise)]]$fari[method, iter] <- fari(uTrue, res$u)$fari
            resultsByNoise[[paste0("noise", currentNoise)]]$accuracy[method, iter] <- 
              clustering_accuracy(yTrue, res$cluster)
            resultsByNoise[[paste0("noise", currentNoise)]]$entropy[method, iter] <- fclust::PE(res$u)
            resultsByNoise[[paste0("noise", currentNoise)]]$xb[method, iter] <- 
              fclust::XB(res$evecs, res$u, res$centers, m)
          }
        }, error = function(e) {
          if(verbose) message(sprintf("\nError in method %s: %s", method, e$message))
        })
      }
    }
    
    # FCMn
    tryCatch({
      t0 <- Sys.time()
      res <- fclust::FKM(X, k = currk, m = m, RS = 5, maxit = 500)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["FCMn", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["FCMn", iter] <- fari(uTrue, res$U)$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["FCMn", iter] <- 
          clustering_accuracy(yTrue, res$clus[,1])
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["FCMn", iter] <- fclust::PE(res$U)
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["FCMn", iter] <- 
          fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FCMean: %s", e$message))
    })
    
    # GFCm - Gustafson-Kessel
    tryCatch({
      t0 <- Sys.time()
      res <- fclust::FKM.gk(X, k = currk, m = m, RS = 5, maxit = 500)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["GFCm", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["GFCm", iter] <- fari(uTrue, res$U)$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["GFCm", iter] <- 
          clustering_accuracy(yTrue, res$clus[,1])
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["GFCm", iter] <- fclust::PE(res$U)
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["GFCm", iter] <- 
          fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in GK: %s", e$message))
    })
    
    # MFCm - Fuzzy C-Medoids
    tryCatch({
      t0 <- Sys.time()
      res <- fclust::FKM.med(X, k = currk, m = m, RS = 5, maxit = 500)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["MFCm", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["MFCm", iter] <- fari(uTrue, res$U)$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["MFCm", iter] <- 
          clustering_accuracy(yTrue, res$clus[,1])
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["MFCm", iter] <- fclust::PE(res$U)
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["MFCm", iter] <- 
          fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FCMed: %s", e$message))
    })
    
    # PFCm - Fuzzy C-Means with polynomial fuzzifier
    tryCatch({
      t0 <- Sys.time()
      res <- fclust::FKM.pf(X, k = currk, b = 0.5, RS = 5, maxit = 500)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["PFCm", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["PFCm", iter] <- fari(uTrue, res$U)$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["PFCm", iter] <- 
          clustering_accuracy(yTrue, res$clus[,1])
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["PFCm", iter] <- fclust::PE(res$U)
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["PFCm", iter] <- 
          fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FCPoly: %s", e$message))
    })
    
    # FSSc - Fuzzy Subspace Clustering
    tryCatch({
      t0 <- Sys.time()
      res <- FSC(X, m = m, K = currk)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["FSSc", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["FSSc", iter] <- fari(uTrue, res$U)$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["FSSc", iter] <- 
          clustering_accuracy(yTrue, res$clus[,1])
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["FSSc", iter] <- fclust::PE(res$U)
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["FSSc", iter] <- 
          fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FSC: %s", e$message))
    })
    
    # FEc - Fuzzy Entropy Clustering
    tryCatch({
      t0 <- Sys.time()
      res <- FEC(X, m = m, K = currk)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["FEc", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["FEc", iter] <- fari(uTrue, res$U)$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["FEc", iter] <- 
          clustering_accuracy(yTrue, res$clus[,1])
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["FEc", iter] <- fclust::PE(res$U)
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["FEc", iter] <- 
          fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FEC: %s", e$message))
    })
    
    # GMM
    tryCatch({
      t0 <- Sys.time()
      res <- mclust::Mclust(data = X, G = currk)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["GMM", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["GMM", iter] <- fari(uTrue, res$z)$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["GMM", iter] <- 
          clustering_accuracy(yTrue, res$classification)
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["GMM", iter] <- fclust::PE(res$z)
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["GMM", iter] <- 
          fclust::XB(X, res$z, t(res$parameters$mean), 2)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in GMM: %s", e$message))
    })
    
    # VWKFC
    tryCatch({
      t0 <- Sys.time()
      res <- vwkfc(X, C = currk, m = m, sigma = 1, gamma = 1, eps = 1e-6, iM = 100, verbose = FALSE)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["VWKFC", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["VWKFC", iter] <- fari(uTrue, t(res$U))$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["VWKFC", iter] <- 
          clustering_accuracy(yTrue, apply(t(res$U), 1, which.max))
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["VWKFC", iter] <- fclust::PE(t(res$U))
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["VWKFC", iter] <- 
          fclust::XB(X, t(res$U), res$W, 2)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in VWKFC: %s", e$message))
    })
    
    # FCMM
    tryCatch({
      t0 <- Sys.time()
      res <- fcmm(X, c = currk, q = 20, r = 2, alpha = 1,
                  max_iter = 200, tol = 1e-6, init = "kmeans",
                  verbose = FALSE)
      resultsByNoise[[paste0("noise", currentNoise)]]$runtime["FCMM", iter] <- 
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsByNoise[[paste0("noise", currentNoise)]]$fari["FCMM", iter] <- fari(uTrue, res$L)$fari
        resultsByNoise[[paste0("noise", currentNoise)]]$accuracy["FCMM", iter] <- 
          clustering_accuracy(yTrue, apply(res$L, 1, which.max))
        resultsByNoise[[paste0("noise", currentNoise)]]$entropy["FCMM", iter] <- fclust::PE(res$L)
        resultsByNoise[[paste0("noise", currentNoise)]]$xb["FCMM", iter] <- 
          fclust::XB(X, res$L, res$V, 2)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FCMM: %s", e$message))
    })
    
    if (iter %% saveInterval == 0) {
      saveRDS(resultsByNoise,
              file = sprintf("progressiveNoiseDims_n%d_iter%d.rds", currentNoise, iter))
      if(verbose) cat(sprintf("\n Save checkpoint at noise_dims=%d, iter=%d \n", currentNoise, iter))
    }
  }

  saveRDS(resultsByNoise,
          file = sprintf("progressiveNoiseDimsComplete_n%d.rds", currentNoise))
  if(verbose) cat(sprintf("\n Completed noise_dims=%d \n", currentNoise))
}


# FINAL RESULTS
finalResults <- list(
  resultsByNoise = resultsByNoise,
  parameters = list(
    numItersPerNoise = numItersPerNoise,
    minNoiseDims = minNoiseDims,
    maxNoiseDims = maxNoiseDims,
    nPointsPerCluster = nPointsPerCluster,
    clusterSeparation = clusterSeparation,
    clusterSigma = clusterSigma,
    noiseRange = noiseRange,
    m = m,
    methodNames = methodNames,
    baseSeed = baseSeed,
    date = Sys.Date()
  )
)

saveRDS(finalResults, file = "finalGaussianNoiseResults.rds")