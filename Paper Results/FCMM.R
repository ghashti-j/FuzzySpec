######## the following script allows us to the run the FCMM algorithm
######## Xiaojun Yang; Mingjun Zhu; Bo Sun; Zheng Wang; Feiping Nie

# X. Yang, M. Zhu, B. Sun, Z. Wang and F. Nie, 
# "Fuzzy C-Multiple-Means Clustering for Hyperspectral Image," 
# in IEEE Geoscience and Remote Sensing Letters, 
# vol. 20, pp. 1-5, 2023, Art no. 5503205, doi: 10.1109/LGRS.2023.3246633.

fcmm <- function(X, c, q, r = 2, alpha = 1,
                 max_iter = 200, tol = 1e-6,
                 init = c("random", "kmeans"),
                 seed = NULL, verbose = FALSE) {
  init <- match.arg(init)
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  if (c < 2) stop("c must be >= 2.")
  if (q < 2) stop("q must be >= 2.")
  if (r <= 1) stop("r must be > 1.")
  if (alpha < 0) stop("alpha must be >= 0.")
  eps <- 1e-12
  # initialization of M and V then Y and Z 
  if (init == "kmeans") {
    kmQ <- stats::kmeans(X, centers = q, nstart = 5)
    M <- kmQ$centers
    kmCentres <- stats::kmeans(M, centers = c, nstart = 5)
    V <- kmCentres$centers
  } else {
    idxM <- sample.int(n, q, replace = FALSE)
    M <- X[idxM, , drop = FALSE]
    idxV <- sample.int(q, c, replace = FALSE)
    V <- M[idxV, , drop = FALSE]
  }
  # compute squared distances between rows of A and B
  sqdist <- function(A, B) {
    nA <- nrow(A); nB <- nrow(B)
    out <- matrix(0, nrow = nA, ncol = nB)
    for (j in 1:nB) {
      diff <- sweep(A, 2, B[j, ], "-")
      out[, j] <- rowSums(diff^2)
    }
    pmax(out, eps)
  }
  # Eq 7, Y update given M
  updateY <- function(X, M, r) {
    D <- sqdist(X, M)  
    power <- 1/(r-1)
    Y <- matrix(0, nrow = n, ncol = q)
    for (i in 1:n) {
      ratios <- outer(D[i, ], D[i, ], "/") 
      denom <- rowSums(ratios^power)
      Y[i, ] <- 1/denom
    }
    Y
  }
  # Eq 9, Z update given M,V
  updateZ <- function(M, V, r) {
    D <- sqdist(M, V)  # q x c
    power <- 1/(r-1)
    Z <- matrix(0, nrow = q, ncol = c)
    for (k in 1:q) {
      ratios <- outer(D[k, ], D[k, ], "/") # c x c
      denom <- rowSums(ratios^power)
      Z[k, ] <- 1/denom
    }
    Z
  }
  # Eq 11, M update given Y,Z,V
  updateM <- function(X, Y, Z, V, r, alpha) {
    Yr <- Y^r          
    Zr <- Z^r          
    newM <- matrix(0, nrow = q, ncol = d)
    # numerator1k
    num1 <- t(Yr)%*%X          
    den1 <- colSums(Yr)            
    # numerator2k 
    num2 <- alpha*(Zr%*%V)     
    den2 <- alpha*rowSums(Zr)   
    den <- den1+den2
    den[den == 0] <- eps
    newM <- sweep(num1+num2, 1, den, "/")
    newM
  }
  # Eq 12, V update given M,Z
  updateV <- function(M, Z, r) {
    Zr <- Z^r            
    num <- t(Zr)%*%M   
    den <- rowSums(t(Zr))
    den[den == 0] <- eps
    sweep(num, 1, den, "/")
  }
  # Objective (Eq 5)
  objective <- function(X, M, Y, V, Z, r, alpha) {
    D_xm <- sqdist(X, M)    
    D_mv <- sqdist(M, V)   
    J1 <- sum((Y^r)*D_xm)
    J2 <- sum((Z^r)*D_mv)
    J1+alpha*J2
  }
  # Initialize Y and Z once 
  Y <- updateY(X, M, r)
  Z <- updateZ(M, V, r)
  
  prevJ <- objective(X, M, Y, V, Z, r, alpha)
  
  for (it in 1:max_iter) {
    M <- updateM(X, Y, Z, V, r, alpha)
    V <- updateV(M, Z, r)
    Y <- updateY(X, M, r)
    Z <- updateZ(M, V, r)
    J <- objective(X, M, Y, V, Z, r, alpha)
    if (verbose) message(sprintf("iter=%d  J=%.3f  |Jprev-Jnew|=%.3e", it, J, abs(prevJ-J)))
    if (abs(prevJ-J) < tol) break
    prevJ <- J
  }
  Lmat <- Y%*%Z  
  list(L = Lmat, Y = Y, Z = Z, M = M, V = V, objective = prevJ)
}
