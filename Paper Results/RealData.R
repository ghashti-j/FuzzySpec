################################### Processing Real Data
#### all data is stored in the saved RDS file, 
#### the commented lines generated the saved datasets datasets.rds
#### upon running the below code, all datasets will be processed, with
#### all results saved in the list "resList" with additional rds files 
#### saved with the results.


##### NOTE #####
# at line 134, the code
# bw <- npudensbw(X, bwmethod = "normal-reference", nmulti = 2)
# is included to increase computation times
# change to bwmethod = "cv.ls" to obtain results as in the paper
source("VWKFC.R")
source("FCMM.R")
source("FSSC.R")
source("FEc.R")
source("SNNs.R")
source("LDAs.R")
source("Functions.R")
datasets <- readRDS("datasets.rds")

# iono <- read.csv("Data/ionosphere.data", header = F)
# iono <- iono[,-2]
# iono.cl <- iono[,34]
# iono.u <- determine.u(iono.cl)
# iono.r <- data.matrix(iono[,-34])
# iono.s <- data.matrix(scale(iono[,-34]))
# 
# iris <- read.csv("Data/iris.data", header = F)
# iris.cl <- iris[,5]
# iris.u <- determine.u(iris.cl)
# iris.r <- data.matrix(iris[,-5])
# iris.s <- data.matrix(scale(iris[,-5]))
# 
# isol <- read.csv("Data/isolet.data", header = F)
# isol.cl <- isol[,618]
# isol.u <- determine.u(isol.cl)
# isol.r <- data.matrix(isol[,-618])
# isol.s <- data.matrix(scale(isol[,-618]))
# 
# lett <- read.csv("Data/letter.data", header = F)
# lett <- lett[lett$V1 %in% c("A", "E", "I", "O", "U", "Y"), ]
# lett.cl <- lett$V1
# lett.u <- determine.u(lett.cl)
# lett.r <- data.matrix(lett[,-1])
# lett.s <- data.matrix(scale(lett[,-1]))
# 
# seed <- read.table("Data/seed.txt", header = F)
# seed.cl <- seed$V8
# seed.u <- determine.u(seed.cl)
# seed.r <- data.matrix(seed[,-8])
# seed.s <- data.matrix(scale(seed[,-8]))
# 
# segm <- read.csv("Data/segmentation.data", header = F)
# segm.cl <- segm$V1
# segm.u <- determine.u(segm.cl)
# segm.r <- data.matrix(segm[,-c(1,3,4,5)])
# segm.s <- data.matrix(scale(segm[,-c(1,3,4,5)]))
# 
# spec <- read.csv("Data/spect.data", header = F)
# spec.cl <- spec$V1
# spec.u <- determine.u(spec.cl)
# spec.r <- data.matrix(spec[,-1])
# spec.s <- data.matrix(scale(spec[,-1]))
# 
# wine <- read.csv("Data/wine.data", header = F)
# wine.cl <- wine$V1
# wine.u <- determine.u(wine.cl)
# wine.r <- data.matrix(wine[,-1])
# wine.s <- data.matrix(scale(wine[,-1]))
# 
# yeas <- read.table("Data/yeast.data", header = F)
# yeas.cl <- yeas$V10
# yeas.u <- determine.u(yeas.cl)
# yeas.r <- data.matrix(yeas[,-c(1,10)])
# yeas.s <- data.matrix(scale(yeas[,-c(1,10)]))
# 
# wdbc <- read.table("Data/wdbc.data", header = F, sep = ",")
# wdbc.cl <- wdbc$V2
# wdbc.u <- determine.u(wdbc.cl)
# wdbc.r <- data.matrix(wdbc[,-c(1,2)])
# wdbc.s <- data.matrix(scale(wdbc[,-c(1,2)]))
# 
# bank <- read.table("Data/bank.txt", header = F, sep = ",")
# bank$V5 <- bank$V5 + 1
# bank.cl <- bank$V5
# bank.u <- determine.u(bank.cl)
# bank.r <- data.matrix(bank[,-5])
# bank.s <- data.matrix(scale(bank[,-5]))
# 
# datasets <- list(
#   list(name = "bank", X = bank.r, X_scaled = bank.s, U = bank.u, cl = bank.cl, k = ncol(bank.u)),
#   list(name = "iono", X = iono.r, X_scaled = iono.s, U = iono.u, cl = iono.cl, k = ncol(iono.u)),
#   list(name = "iris", X = iris.r, X_scaled = iris.s, U = iris.u, cl = iris.cl, k = ncol(iris.u)),
#   list(name = "lett", X = lett.r, X_scaled = lett.s, U = lett.u, cl = lett.cl, k = ncol(lett.u)),
#   list(name = "seed", X = seed.r, X_scaled = seed.s, U = seed.u, cl = seed.cl, k = ncol(seed.u)),
#   list(name = "segm", X = segm.r, X_scaled = segm.s, U = segm.u, cl = segm.cl, k = ncol(segm.u)),
#   list(name = "spec", X = spec.r, X_scaled = spec.s, U = spec.u, cl = spec.cl, k = ncol(spec.u)),
#   list(name = "wdbc", X = wdbc.r, X_scaled = wdbc.s, U = wdbc.u, cl = wdbc.cl, k = ncol(wdbc.u)),
#   list(name = "wine", X = wine.r, X_scaled = wine.s, U = wine.u, cl = wine.cl, k = ncol(wine.u)),
#   list(name = "yeas", X = yeas.r, X_scaled = yeas.s, U = yeas.u, cl = yeas.cl, k = ncol(yeas.u))
#  )

# saveRDS(datasets, "datasets.rds")
# datasets <- readRDS("datasets.rds")

# define method names
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

# process a single dataset
processData <- function(X, data.u, data.cl, k, m = 2, dataset_name, verbose = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  dataResults <- list(
    fari = numeric(length(methodNames)),
    PE = numeric(length(methodNames)),
    XB = numeric(length(methodNames)),
    accuracy = numeric(length(methodNames))
  )
  names(dataResults$fari) <- methodNames
  names(dataResults$PE) <- methodNames
  names(dataResults$XB) <- methodNames
  names(dataResults$accuracy) <- methodNames
  
  # Compute bandwidth
  tryCatch({
    bw <- npudensbw(X, bwmethod = "normal-reference", nmulti = 2)
    h <- bw$bw
    w <- 1/h^2
  }, error = function(e) {
    h <- apply(X, 2, sd) * n^(-1/(4+p))
    w <- 1/h^2
  })
  Xsc <- sweep(X, 2, sqrt(w), "*")
  
  # distance and similarity matrices
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
  sigmaks <- compute.sigma(D.kss, "vw")
  r.ks <- sigmaks$radius
  sigma.ks <- sigmaks$sigma
  snn.ks <- compute.SNN(S.kss, r.ks)
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
  
  # Process spectral methods 
  spectralMethods <- methodNames[!methodNames %in% c("FCMn", "GFCm", "MFCm", "PFCm", "FSSc", "FEc")]
  adjacencyMatrix <- list()
  for(method in spectralMethods){
    if(verbose) paste0(method)
    tryCatch({
      adjacencyMatrix[[method]] <- make.similarity(data = X, method = method, n = n, 
                                         sigma.eu = sigma.eu, sigma.ks = sigma.ks, 
                                         rho.eu = rho.eu, rho.ks = rho.ks, 
                                         snn.eu = snn.eu, snn.ks = snn.ks, 
                                         S.kss = S.kss, S.eu = S.eu, 
                                         D.kss.sq = D.kss.sq, 
                                         D.eu.sq = D.eu.sq, D.eu = D.eu)
    }, error = function(e) { 
      adjacencyMatrix[[method]] <- NULL 
      message(sprintf("Error creating similarity matrix for %s: %s", method, e$message))
    })
  }
  
  # Run spectral methods
  for (method in names(adjacencyMatrix)) {
    if(verbose) paste0(method)
    if (!is.null(adjacencyMatrix[[method]])) {
      tryCatch({
        res <- fuzzy.spectral.clustering(adjacencyMatrix[[method]], k = k, m = m, method = "CM")
        if (!is.null(res)) {
          dataResults$fari[method] <- fari(data.u, res$u)$fari
          dataResults$accuracy[method] <- clustering_accuracy(data.cl, res$cluster)
          dataResults$PE[method] <- fclust::PE(res$u)
          dataResults$XB[method] <- fclust::XB(res$evecs, res$u, res$centers, m)
        }
      }, error = function(e) {
        message(sprintf("Error in spectral method %s: %s", method, e$message))
        dataResults$fari[method] <- NA
        dataResults$accuracy[method] <- NA
        dataResults$PE[method] <- NA
        dataResults$XB[method] <- NA
      })
    } else {
      dataResults$fari[method] <- NA
      dataResults$accuracy[method] <- NA
      dataResults$PE[method] <- NA
      dataResults$XB[method] <- NA
    }
  }
  
  # FCMn - Fuzzy C-Means
  tryCatch({
    res <- fclust::FKM(X, k = k, m = m, RS = 5, maxit = 500)
    dataResults$fari["FCMn"] <- fari(data.u, res$U)$fari
    dataResults$accuracy["FCMn"] <- clustering_accuracy(data.cl, res$clus[,1])
    dataResults$PE["FCMn"] <- fclust::PE(res$U)
    dataResults$XB["FCMn"] <- fclust::XB(X, res$U, res$H, m)
  }, error = function(e) {
    message(sprintf("Error in FCMn: %s", e$message))
    dataResults$fari["FCMn"] <- NA
    dataResults$accuracy["FCMn"] <- NA
    dataResults$PE["FCMn"] <- NA
    dataResults$XB["FCMn"] <- NA
  })
  
  # GFCm - Gustafson-Kessel
  tryCatch({
    res <- fclust::FKM.gk(X, k = k, m = m, RS = 5, maxit = 500)
    dataResults$fari["GFCm"] <- fari(data.u, res$U)$fari
    dataResults$accuracy["GFCm"] <- clustering_accuracy(data.cl, res$clus[,1])
    dataResults$PE["GFCm"] <- fclust::PE(res$U)
    dataResults$XB["GFCm"] <- fclust::XB(X, res$U, res$H, m)
  }, error = function(e) {
    message(sprintf("Error in GFCm: %s", e$message))
    dataResults$fari["GFCm"] <- NA
    dataResults$accuracy["GFCm"] <- NA
    dataResults$PE["GFCm"] <- NA
    dataResults$XB["GFCm"] <- NA
  })
  
  # MFCm - Fuzzy C-Medoids
  tryCatch({
    res <- fclust::FKM.med(X, k = k, m = m, RS = 5, maxit = 500)
    dataResults$fari["MFCm"] <- fari(data.u, res$U)$fari
    dataResults$accuracy["MFCm"] <- clustering_accuracy(data.cl, res$clus[,1])
    dataResults$PE["MFCm"] <- fclust::PE(res$U)
    dataResults$XB["MFCm"] <- fclust::XB(X, res$U, res$H, m)
  }, error = function(e) {
    message(sprintf("Error in MFCm: %s", e$message))
    dataResults$fari["MFCm"] <- NA
    dataResults$accuracy["MFCm"] <- NA
    dataResults$PE["MFCm"] <- NA
    dataResults$XB["MFCm"] <- NA
  })
  
  # PFCm - Fuzzy C-Means with polynomial fuzzifier
  tryCatch({
    res <- fclust::FKM.pf(X, k = k, b = 0.5, RS = 5, maxit = 500)
    dataResults$fari["PFCm"] <- fari(data.u, res$U)$fari
    dataResults$accuracy["PFCm"] <- clustering_accuracy(data.cl, res$clus[,1])
    dataResults$PE["PFCm"] <- fclust::PE(res$U)
    dataResults$XB["PFCm"] <- fclust::XB(X, res$U, res$H, m)
  }, error = function(e) {
    message(sprintf("Error in PFCm: %s", e$message))
    dataResults$fari["PFCm"] <- NA
    dataResults$accuracy["PFCm"] <- NA
    dataResults$PE["PFCm"] <- NA
    dataResults$XB["PFCm"] <- NA
  })
  
  # FSSc fuzzy subspace clustering
  tryCatch({
    res <- FSC(X, m = m, K = k)
    dataResults$fari["FSSc"] <- fari(data.u, res$U)$fari
    dataResults$accuracy["FSSc"] <- clustering_accuracy(data.cl, res$clus[,1])
    dataResults$PE["FSSc"] <- fclust::PE(res$U)
    dataResults$XB["FSSc"] <- fclust::XB(X, res$U, res$H, m)
  }, error = function(e) {
    message(sprintf("Error in FSSc: %s", e$message))
    dataResults$fari["FSSc"] <- NA
    dataResults$accuracy["FSSc"] <- NA
    dataResults$PE["FSSc"] <- NA
    dataResults$XB["FSSc"] <- NA
  })
  
  # fuzzy entropy clustering
  tryCatch({
    res <- FEC(X, m = m, K = k)
    dataResults$fari["FEc"] <- fari(data.u, res$U)$fari
    dataResults$accuracy["FEc"] <- clustering_accuracy(data.cl, res$clus[,1])
    dataResults$PE["FEc"] <- fclust::PE(res$U)
    dataResults$XB["FEc"] <- fclust::XB(X, res$U, res$H, m)
  }, error = function(e) {
    message(sprintf("Error in FEc: %s", e$message))
    dataResults$fari["FEc"] <- NA
    dataResults$accuracy["FEc"] <- NA
    dataResults$PE["FEc"] <- NA
    dataResults$XB["FEc"] <- NA
  })
  
  
  # GMM
  tryCatch({
    res <- mclust::Mclust(X, G = k)
    dataResults$fari["GMM"] <- fari(data.u, res$z)$fari
    dataResults$accuracy["GMM"] <- clustering_accuracy(data.cl, res$classification)
    dataResults$PE["GMM"] <- fclust::PE(res$z)
    dataResults$XB["GMM"] <- fclust::XB(X, res$z, t(res$parameters$mean), 2)
  }, error = function(e) {
    message(sprintf("Error in GMM: %s", e$message))
    dataResults$fari["GMM"] <- NA
    dataResults$accuracy["GMM"] <- NA
    dataResults$PE["GMM"] <- NA
    dataResults$XB["GMM"] <- NA
  })
  
  
  # VWKFC
  tryCatch({
    res <- vwkfc(X, C = k, m = m, sigma = 1, gamma = 1, eps = 1e-6, iM = 100, verbose = FALSE)
    dataResults$fari["VWKFC"] <- fari(data.u, t(res$U))$fari
    dataResults$accuracy["VWKFC"] <- clustering_accuracy(data.cl, apply(t(res$U), 1, which.max))
    dataResults$PE["VWKFC"] <- fclust::PE(t(res$U))
    dataResults$XB["VWKFC"] <- fclust::XB(X, t(res$U), res$W, 2)
  }, error = function(e) {
    message(sprintf("Error in VWKFC: %s", e$message))
    dataResults$fari["VWKFC"] <- NA
    dataResults$accuracy["VWKFC"] <- NA
    dataResults$PE["VWKFC"] <- NA
    dataResults$XB["VWKFC"] <- NA
  })
  
  # FCMM
  tryCatch({
    res <- fcmm(X, c = k, q = 20, r = 2, alpha = 1,
                max_iter = 200, tol = 1e-6, init = "kmeans",
                verbose = FALSE)
    dataResults$fari["FCMM"] <- fari(data.u, res$L)$fari
    dataResults$accuracy["FCMM"] <- clustering_accuracy(data.cl, apply(res$L, 1, which.max))
    dataResults$PE["FCMM"] <- fclust::PE(res$L)
    dataResults$XB["FCMM"] <- fclust::XB(X, res$L, res$V, 2)
  }, error = function(e) {
    message(sprintf("Error in FCMM: %s", e$message))
    dataResults$fari["FCMM"] <- NA
    dataResults$accuracy["FCMM"] <- NA
    dataResults$PE["FCMM"] <- NA
    dataResults$XB["FCMM"] <- NA
  })
  
  return(dataResults)
}

# Main processing loop
resList <- list()
itercount <- 1
for (dataset in datasets[3]) {
  print(paste0("dataset", itercount,"(unscaled version) of ", length(datasets)))
  resList[[paste0(dataset$name, ".r")]] <- processData(
    X = dataset$X,
    data.u = dataset$U,
    data.cl = dataset$cl,
    k = dataset$k,
    dataset_name = paste0(dataset$name, ".r")
  )
  saveRDS(resList, paste0("realres_raw", itercount, ".rds"))
  print(paste0("dataset", itercount,"(scaled version) of ", length(datasets)))
  resList[[paste0(dataset$name, ".s")]] <- processData(
    X = dataset$X_scaled,
    data.u = dataset$U,
    data.cl = dataset$cl,
    k = dataset$k,
    dataset_name = paste0(dataset$name, ".s")
  )
  saveRDS(resList, paste0("realres_scaled", itercount, ".rds"))
  itercount <- itercount + 1
}


