# SA2 generation, visualization and results from paper

# required functions
source("VWKFC.R")
source("FCMM.R")
source("FSSC.R")
source("FEc.R")
source("SNNs.R")
source("LDAs.R")
source("Functions.R")

# generation function 
genOverlapData <- function(nPoints = 500, k = 5, separationFactor = 4.0,
                          seed = NULL,
                          baseSigma = 0.5,
                          arrangement = "circle") {
  if (!is.null(seed)) set.seed(seed)
  nPerCluster <- floor(nPoints / k)
  nRemaining <- nPoints - k * nPerCluster
  clusterSizes <- rep(nPerCluster, k)
  if (nRemaining > 0) {
    clusterSizes[1:nRemaining] <- clusterSizes[1:nRemaining] + 1
  }

  if (arrangement == "circle") {
    angles <- seq(0, 2*pi, length.out = k + 1)[1:k]
    baseCenters <- cbind(cos(angles), sin(angles))
  } else if (arrangement == "line") {
    positions <- seq(-1, 1, length.out = k)
    baseCenters <- cbind(positions, rep(0, k))
  } else if (arrangement == "grid") {
    gridSize <- ceiling(sqrt(k))
    xPos <- rep(seq(-1, 1, length.out = gridSize), gridSize)[1:k]
    yPos <- rep(seq(-1, 1, length.out = gridSize), each = gridSize)[1:k]
    baseCenters <- cbind(xPos, yPos)
  } else stop("arrangement must be 'circle', 'line', or 'grid'")
  
  
  centers <- baseCenters * separationFactor
  covMatrix <- diag(2) * baseSigma^2
  datlist <- list()
  for (i in 1:k) {
    datlist[[i]] <- mvtnorm::rmvnorm(clusterSizes[i],
                                     mean = centers[i,],
                                     sigma = covMatrix)
  }

  X <- do.call("rbind", datlist)
  nTotal <- nrow(X)
  yTrue <- unlist(mapply(rep, 1:k, clusterSizes))
  uFuzzy <- matrix(0, nTotal, k)
  piVec <- clusterSizes / nTotal
  densities <- matrix(0, nTotal, k)
  for (j in 1:k) densities[, j] <- mvtnorm::dmvnorm(X, mean = centers[j,], sigma = covMatrix)
  weightedDensities <- sweep(densities, 2, piVec, "*")
  totalDensity <- rowSums(weightedDensities)
  zeroDensity <- totalDensity == 0
  if (any(zeroDensity)) {
    for (i in which(zeroDensity)) {
      distances <- sqrt(rowSums((centers - matrix(X[i,], nrow = k, ncol = 2, byrow = TRUE))^2))
      nearest <- which.min(distances)
      uFuzzy[i, nearest] <- 1
    }
  }

  uFuzzy[!zeroDensity, ] <- weightedDensities[!zeroDensity, ] / totalDensity[!zeroDensity]
  minDistance <- Inf
  avgDistance <- 0
  nPairs <- 0
  
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      d <- sqrt(sum((centers[i,] - centers[j,])^2))
      minDistance <- min(minDistance, d)
      avgDistance <- avgDistance + d
      nPairs <- nPairs + 1
    }
  }
  avgDistance <- avgDistance / nPairs
  
  # overlap metric
  overlapMetric <- mean(apply(uFuzzy, 1, function(row) {
    trueCluster <- which.max(row)
    max(row[-trueCluster])
  }))
  
  # effective overlap 
  overlapThreshold <- 4 * baseSigma
  overlapPairs <- sum(as.matrix(dist(centers)) < overlapThreshold &
                        as.matrix(dist(centers)) > 0) / 2
  overlapProportion <- overlapPairs / (k * (k-1) / 2)
  
  return(list(
    X = X,
    uTrue = uFuzzy,
    yTrue = yTrue,
    k = k,
    centers = centers,
    clusterSizes = clusterSizes,
    sigma = baseSigma,
    separationFactor = separationFactor,
    minClusterDistance = minDistance,
    avgClusterDistance = avgDistance,
    overlapMetric = overlapMetric,
    overlapProportion = overlapProportion,
    covMatrix = covMatrix
  ))
}

# visualization function
vizSepDecay <- function(data, title = NULL) {
  df <- data.frame(
    x = data$X[,1],
    y = data$X[,2],
    cluster = as.factor(data$yTrue),
    u = apply(data$uTrue, 1, max)
  )
  centersDF <- data.frame(
    x = data$centers[,1],
    y = data$centers[,2],
    cluster = as.factor(1:data$k)
  )
  xRange <- range(data$X[,1])
  yRange <- range(data$X[,2])
  xExpand <- diff(xRange) * 0.1
  yExpand <- diff(yRange) * 0.1
  colours <- c("#1F78B4", "#33A02C", "#FB9A99", "#E31A1C", "#c497cf", "magenta4", "orange4")

  p <- ggplot(df, aes(x = x, y = y, color = cluster)) +
    geom_point(aes(size = u), alpha = 0.6) +
    geom_point(data = centersDF, aes(x = x, y = y),
               color = "black", size = 4, shape = 4) +
    scale_color_manual(values = colours, guide = "none") +
    scale_size_continuous(range = c(5, 1.5)) +
    coord_equal() +
    xlim(xRange[1] - xExpand, xRange[2] + xExpand) +
    ylim(yRange[1] - yExpand, yRange[2] + yExpand) +
    labs(x = expression(X[1]),
         y = expression(X[2]),
         color = expression(C[j]),
         size = expression(U[i*","*C[j]])) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          panel.background = element_rect(NA, "black", 1),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.box.margin = margin(t = -4, r = 0, b = 0, l = 0, unit = "mm"))

  for (i in 1:data$k) {
    circleData <- data.frame(
      x = data$centers[i,1],
      y = data$centers[i,2],
      r = 2 * data$sigma)
    p <- p + geom_circle(data = circleData,
                         aes(x0 = x, y0 = y, r = r),
                         color = colours[i],
                         fill = NA, linetype = "dashed",
                         alpha = 0.5, inherit.aes = FALSE)
  }
  return(p)
}

# PARAMETERS
numItersPerSep <- 1000
sepSequence <- seq(2.5, 0.5, -0.25)  # 5 separation levels: 2.5, 2.0, 1.5, 1.0, 0.5
numPts <- 500
currk <- 5  
m <- 2
saveInterval <- 20
baseSeed <- 44869410

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
resultsBySep <- list()
for (sep in sepSequence) {
  sepKey <- paste0("sep", gsub("\\.", "_", as.character(sep)))
  resultsBySep[[sepKey]] <- list(
    accuracy = matrix(NA, nrow = numMethods, ncol = numItersPerSep,
                      dimnames = list(methodNames, paste0(1:numItersPerSep))),
    fari = matrix(NA, nrow = numMethods, ncol = numItersPerSep,
                  dimnames = list(methodNames, paste0(1:numItersPerSep))),
    runtime = matrix(NA, nrow = numMethods, ncol = numItersPerSep,
                     dimnames = list(methodNames, paste0(1:numItersPerSep))),
    entropy = matrix(NA, nrow = numMethods, ncol = numItersPerSep,
                     dimnames = list(methodNames, paste0(1:numItersPerSep))),
    xb = matrix(NA, nrow = numMethods, ncol = numItersPerSep,
                dimnames = list(methodNames, paste0(1:numItersPerSep))),
    rholist = list(),
    sigmalist = list(),
    radiusMat = matrix(0, 2, numItersPerSep),
    snntimes = numeric(numItersPerSep),
    sigma.times = numeric(numItersPerSep),
    bw.times = numeric(numItersPerSep),
    overlapMetrics = numeric(numItersPerSep),
    minDistances = numeric(numItersPerSep)
  )
}
rownames(resultsBySep[[1]]$radiusMat) <- c("eu", "vw")


# MAIN EXPERIMENT
overallIter <- 0
totalIterations <- numItersPerSep * length(sepSequence)

for (sepIdx in seq_along(sepSequence)) {
  currentSep <- sepSequence[sepIdx]
  sepKey <- paste0("sep", gsub("\\.", "_", as.character(currentSep)))
  
  if(verbose) cat(sprintf("\n Separation = %.1f \n", currentSep))
  
  for (iter in 1:numItersPerSep) {
    overallIter <- overallIter + 1
    iterSeed <- baseSeed + iter
    data <- genOverlapData(
      nPoints = numPts,
      k = currk,
      separationFactor = currentSep,
      seed = iterSeed,
      baseSigma = 0.5,
      arrangement = "circle"
    )
    
    X <- data$X
    uTrue <- data$uTrue
    yTrue <- data$yTrue
    n <- nrow(X)
    p <- ncol(X)
    
    resultsBySep[[sepKey]]$overlapMetrics[iter] <- data$overlapMetric
    resultsBySep[[sepKey]]$minDistances[iter] <- data$minClusterDistance
    
    ################## NOTE ###################################
    # the code
    # bw <- npudensbw(X, bwmethod = "normal-reference", nmulti = 2)
    # is included to increase computation times
    # change to bwmethod = "cv.ls" to obtain results as in the paper
    tryCatch({
      t0 <- Sys.time()
      bw <- npudensbw(X, bwmethod = "normal-reference", nmulti = 2)
      resultsBySep[[sepKey]]$bw.times[iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
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
    resultsBySep[[sepKey]]$sigma.times[iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    r.ks <- sigmaks$radius
    sigma.ks <- sigmaks$sigma
    
    t0 <- Sys.time()
    snn.ks <- compute.SNN(S.kss, r.ks)
    resultsBySep[[sepKey]]$snntimes[iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
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
    
    resultsBySep[[sepKey]]$rholist[[iter]] <- list(eu = rho.eu, ks = rho.ks)
    resultsBySep[[sepKey]]$sigmalist[[iter]] <- list(eu = sigma.eu, ks = sigma.ks)
    resultsBySep[[sepKey]]$radiusMat[,iter] <- c(r.eu, r.ks)
    
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
          resultsBySep[[sepKey]]$runtime[method, iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
          if (!is.null(res)) {
            resultsBySep[[sepKey]]$fari[method, iter] <- fari(uTrue, res$u)$fari
            resultsBySep[[sepKey]]$accuracy[method, iter] <- clustering_accuracy(yTrue, res$cluster)
            resultsBySep[[sepKey]]$entropy[method, iter] <- fclust::PE(res$u)
            resultsBySep[[sepKey]]$xb[method, iter] <- fclust::XB(res$evecs, res$u, res$centers, m)
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
      resultsBySep[[sepKey]]$runtime["FCMn", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["FCMn", iter] <- fari(uTrue, res$U)$fari
        resultsBySep[[sepKey]]$accuracy["FCMn", iter] <- clustering_accuracy(yTrue, res$clus[,1])
        resultsBySep[[sepKey]]$entropy["FCMn", iter] <- fclust::PE(res$U)
        resultsBySep[[sepKey]]$xb["FCMn", iter] <- fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FCMean: %s", e$message))
    })
    
    # GFCm - Gustafson-Kessel
    tryCatch({
      t0 <- Sys.time()
      res <- fclust::FKM.gk(X, k = currk, m = m, RS = 5, maxit = 500)
      resultsBySep[[sepKey]]$runtime["GFCm", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["GFCm", iter] <- fari(uTrue, res$U)$fari
        resultsBySep[[sepKey]]$accuracy["GFCm", iter] <- clustering_accuracy(yTrue, res$clus[,1])
        resultsBySep[[sepKey]]$entropy["GFCm", iter] <- fclust::PE(res$U)
        resultsBySep[[sepKey]]$xb["GFCm", iter] <- fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in GK: %s", e$message))
    })
    
    # MFCm - Fuzzy C-Medoids
    tryCatch({
      t0 <- Sys.time()
      res <- fclust::FKM.med(X, k = currk, m = m, RS = 5, maxit = 500)
      resultsBySep[[sepKey]]$runtime["MFCm", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["MFCm", iter] <- fari(uTrue, res$U)$fari
        resultsBySep[[sepKey]]$accuracy["MFCm", iter] <- clustering_accuracy(yTrue, res$clus[,1])
        resultsBySep[[sepKey]]$entropy["MFCm", iter] <- fclust::PE(res$U)
        resultsBySep[[sepKey]]$xb["MFCm", iter] <- fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FCMed: %s", e$message))
    })
    
    # PFCm - Fuzzy C-Means with polynomial fuzzifier
    tryCatch({
      t0 <- Sys.time()
      res <- fclust::FKM.pf(X, k = currk, b = 0.5, RS = 5, maxit = 500)
      resultsBySep[[sepKey]]$runtime["PFCm", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["PFCm", iter] <- fari(uTrue, res$U)$fari
        resultsBySep[[sepKey]]$accuracy["PFCm", iter] <- clustering_accuracy(yTrue, res$clus[,1])
        resultsBySep[[sepKey]]$entropy["PFCm", iter] <- fclust::PE(res$U)
        resultsBySep[[sepKey]]$xb["PFCm", iter] <- fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FCPoly: %s", e$message))
    })
    
    # FSSc - Fuzzy Subspace Clustering
    tryCatch({
      t0 <- Sys.time()
      res <- FSC(X, m = m, K = currk)
      resultsBySep[[sepKey]]$runtime["FSSc", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["FSSc", iter] <- fari(uTrue, res$U)$fari
        resultsBySep[[sepKey]]$accuracy["FSSc", iter] <- clustering_accuracy(yTrue, res$clus[,1])
        resultsBySep[[sepKey]]$entropy["FSSc", iter] <- fclust::PE(res$U)
        resultsBySep[[sepKey]]$xb["FSSc", iter] <- fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FSC: %s", e$message))
    })
    
    # FEc - Fuzzy Entropy Clustering
    tryCatch({
      t0 <- Sys.time()
      res <- FEC(X, m = m, K = currk)
      resultsBySep[[sepKey]]$runtime["FEc", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["FEc", iter] <- fari(uTrue, res$U)$fari
        resultsBySep[[sepKey]]$accuracy["FEc", iter] <- clustering_accuracy(yTrue, res$clus[,1])
        resultsBySep[[sepKey]]$entropy["FEc", iter] <- fclust::PE(res$U)
        resultsBySep[[sepKey]]$xb["FEc", iter] <- fclust::XB(X, res$U, res$H, m)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FEC: %s", e$message))
    })
    
    # GMM
    tryCatch({
      t0 <- Sys.time()
      res <- mclust::Mclust(data = X, G = currk)
      resultsBySep[[sepKey]]$runtime["GMM", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["GMM", iter] <- fari(uTrue, res$z)$fari
        resultsBySep[[sepKey]]$accuracy["GMM", iter] <- clustering_accuracy(yTrue, res$classification)
        resultsBySep[[sepKey]]$entropy["GMM", iter] <- fclust::PE(res$z)
        resultsBySep[[sepKey]]$xb["GMM", iter] <- fclust::XB(X, res$z, t(res$parameters$mean), 2)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in GMM: %s", e$message))
    })
    
    # VWKFC
    tryCatch({
      t0 <- Sys.time()
      res <- vwkfc(X, C = currk, m = m, sigma = 1, gamma = 1, eps = 1e-6, iM = 100, verbose = FALSE)
      resultsBySep[[sepKey]]$runtime["VWKFC", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["VWKFC", iter] <- fari(uTrue, t(res$U))$fari
        resultsBySep[[sepKey]]$accuracy["VWKFC", iter] <- clustering_accuracy(yTrue, apply(t(res$U), 1, which.max))
        resultsBySep[[sepKey]]$entropy["VWKFC", iter] <- fclust::PE(t(res$U))
        resultsBySep[[sepKey]]$xb["VWKFC", iter] <- fclust::XB(X, t(res$U), res$W, 2)
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
      resultsBySep[[sepKey]]$runtime["FCMM", iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (!is.null(res)) {
        resultsBySep[[sepKey]]$fari["FCMM", iter] <- fari(uTrue, res$L)$fari
        resultsBySep[[sepKey]]$accuracy["FCMM", iter] <- clustering_accuracy(yTrue, apply(res$L, 1, which.max))
        resultsBySep[[sepKey]]$entropy["FCMM", iter] <- fclust::PE(res$L)
        resultsBySep[[sepKey]]$xb["FCMM", iter] <- fclust::XB(X, res$L, res$V, 2)
      }
    }, error = function(e) {
      if(verbose) message(sprintf("\nError in FCMM: %s", e$message))
    })
    
    if (iter %% saveInterval == 0) {
      saveRDS(resultsBySep,
              file = sprintf("separationDecay_sep%.1f_iter%d.rds", currentSep, iter))
      if(verbose) cat(sprintf("\n Save at sep=%.1f, iter=%d \n", currentSep, iter))
    }
  }

  saveRDS(resultsBySep,
          file = sprintf("separationDecayComplete_sep%.1f.rds", currentSep))
  if(verbose) cat(sprintf("\n Completed sep=%.1f \n", currentSep))
}


# FINAL RESULTS
finalResults <- list(
  resultsBySep = resultsBySep,
  parameters = list(
    numItersPerSep = numItersPerSep,
    sepSequence = sepSequence,
    numPts = numPts,
    k = currk,
    baseSigma = 0.5,
    arrangement = "circle",
    m = m,
    methodNames = methodNames,
    baseSeed = baseSeed,
    date = Sys.Date()
  )
)

saveRDS(finalResults, file = "finalGaussianOverlapResults.rds")