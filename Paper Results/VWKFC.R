######## the following script allows us to the run the VWKFC algorithm
######## proposed by Yiming Tang; Zhifu Pan; Witold Pedrycz; Fuji Ren; Xiaocheng Song

# Y. Tang, Z. Pan, W. Pedrycz, F. Ren and X. Song, 
# "Viewpoint-Based Kernel Fuzzy Clustering With Weight Information Granules," 
# in IEEE Transactions on Emerging Topics in Computational Intelligence, 
# vol. 7, no. 2, pp. 342-356, April 2023, doi: 10.1109/TETCI.2022.3201620.


minmax01 <- function(X) {
  # Eq 5 --- min-max scaling 
  X <- as.matrix(X)
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  denom <- maxs - mins
  denom[denom == 0] <- 1
  Xs <- sweep(sweep(X, 2, mins, "-"), 2, denom, "/")
  list(X = Xs, mins = mins, maxs = maxs)
}

# RBF squared distance
rbfPhiDist <- function(a, b, sigma) (2 - 2 * exp(- (a - b)^2 / (2 * sigma^2)))

# Frobenius norm
frNorm <- function(M) sqrt(sum(M * M))


# KHDI Eqs 5,6,8,9,10,11
khdiInit <- function(X, C, sigma = 1, a = 2) {
  X <- as.matrix(X)
  N <- nrow(X); L <- ncol(X)
  
  # Eq 5
  mm <- minmax01(X)
  Xn <- mm$X
  D <- as.matrix(dist(Xn, method = "euclidean"))
  diag(D) <- 0
  
  # Eq 8
  dvals <- D[upper.tri(D)]
  r <- (max(dvals) - min(dvals)) / C
  
  # Eq 6
  rhoPrime <- rowSums(D < r)  
  rhoMin <- min(rhoPrime); rhoMax <- max(rhoPrime)
  rhoDenom <- (rhoMax - rhoMin); if (rhoDenom == 0) rhoDenom <- 1
  rho <- (rhoPrime - rhoMin) / rhoDenom
  
  # Eq 9
  deltaPrime <- rep(NA_real_, N)
  k <- which.max(rho)
  
  for (i in seq_len(N)) {
    higher <- which(rho > rho[i])
    if (length(higher) == 0) {
      deltaPrime[i] <- max(D[i, ])
    } else {
      deltaPrime[i] <- min(D[i, higher])
    }
  }
  deltaMin <- min(deltaPrime); deltaMax <- max(deltaPrime)
  deltaDenom <- (deltaMax - deltaMin); if (deltaDenom == 0) deltaDenom <- 1
  delta <- (deltaPrime - deltaMin) / deltaDenom
  
  # Eq 10
  tau <- rho * delta
  ord <- order(tau, decreasing = TRUE)
  Xs <- Xn[ord, , drop = FALSE]
  
  # Eq 11
  dc <- max(dvals) / (a * C)
  V <- matrix(NA_real_, nrow = 0, ncol = L)
  v1 <- Xs[1, , drop = FALSE]
  V <- rbind(V, v1)
  xd <- as.numeric(v1)
  num <- 1
  index <- 2
  
  while (num < C) {
    if (index > N) break 
    xk <- Xs[index, , drop = FALSE]
    distCentre <- apply(V, 1, function(v) sqrt(sum((xk - v)^2)))
    if (min(distCentre) < dc) {
      index <- index + 1
      next
    }
    V <- rbind(V, xk)
    num <- num + 1
    index <- index + 1
  }
  if (nrow(V) < C) warning("KHDI selected fewer than C centers.")
  list(V = V, xd = xd, Xn = Xn, mm = mm, r = r, dc = dc, tau = tau, ord = ord)
}


# VWKFC updates (Eqs 17,18,19)
compDij <- function(X, G, W, sigma) {
  N <- nrow(X); L <- ncol(X); C <- nrow(G)
  Dij <- matrix(0, nrow = C, ncol = N)
  for (i in seq_len(C)) {
    acc <- rep(0, N)
    for (l in seq_len(L)) {
      d2 <- rbfPhiDist(X[, l], G[i, l], sigma) # length N
      acc <- acc + W[i, l] * d2
    }
    Dij[i, ] <- acc
  }
  Dij
}

updateU <- function(X, G, W, m, sigma, tiny = 1e-12) {
  # Eq 17
  Dij <- compDij(X, G, W, sigma)
  Dij <- pmax(Dij, tiny)
  pow <- -1 / (m - 1)
  A <- Dij ^ pow  
  denom <- colSums(A)
  U <- sweep(A, 2, denom, "/")
  U
}

updateW <- function(X, U, G, m, sigma, gamma, alpha = NULL) {
  # Eq 18
  X <- as.matrix(X)
  N <- nrow(X); L <- ncol(X); C <- nrow(G)
  if (is.null(alpha)) alpha <- rep(1, N)
  if (length(alpha) != N) stop("alpha must have length N.")
  Um <- U ^ m  
  W <- matrix(0, nrow = C, ncol = L)
  for (i in seq_len(C)) {
    scores <- numeric(L)
    wj <- alpha * Um[i, ] 
    for (l in seq_len(L)) {
      d2 <- rbfPhiDist(X[, l], G[i, l], sigma)
      scores[l] <- sum(wj * d2)
    }
    z <- -gamma * scores
    z <- z - max(z)
    ex <- exp(z)
    W[i, ] <- ex / sum(ex)
  }
  W
}

updateG <- function(X, U, Gprev, W, m, sigma, alpha = NULL, xd, q = 1) {
  # Eq 19
  X <- as.matrix(X)
  N <- nrow(X); L <- ncol(X); C <- nrow(Gprev)
  if (is.null(alpha)) alpha <- rep(1, N)
  if (length(alpha) != N) stop("alpha must have length N.")
  if (length(xd) != L) stop("xd must have length L.")
  if (q < 1 || q > C) stop("q must be between 1 and C.")
  Um <- U ^ m
  Gnew <- Gprev
  for (i in seq_len(C)) {
    if (i == q) {
      Gnew[i, ] <- xd
      next
    }
    wjBase <- alpha * Um[i, ]  
    for (l in seq_len(L)) {
      ker <- exp(- (X[, l] - Gprev[i, l])^2 / (2 * sigma^2))
      wil <- W[i, l]
      wj <- wjBase * wil * ker
      denom <- sum(wj)
      if (denom <= 0 || !is.finite(denom)) {
        Gnew[i, l] <- Gprev[i, l]
      } else {
        Gnew[i, l] <- sum(wj * X[, l]) / denom
      }
    }
  }
  Gnew
}


# Main VWKFC 
vwkfc <- function(X, C, m, sigma = 1, gamma = 1, eps = 1e-5, iM = 100,
                  alpha = NULL, a_khdi = 2, q = 1, verbose = FALSE) {
  
  X <- as.matrix(X)
  N <- nrow(X); L <- ncol(X)
  if (m <= 1) stop("m must be > 1.")
  if (C < 2) stop("C must be >= 2.")
  if (eps <= 0) stop("eps must be > 0.")
  if (iM < 1) stop("iM must be >= 1.")
  if (is.null(alpha)) alpha <- rep(1, N)
  if (length(alpha) != N) stop("alpha must have length N.")
  
  # KHDI init
  init <- khdiInit(X, C = C, sigma = sigma, a = a_khdi)
  Xn <- init$Xn
  G <- init$V
  xd <- init$xd
  
  W <- matrix(1 / L, nrow = C, ncol = L) # init weights W
  U <- matrix(0, nrow = C, ncol = N)  # init memberships 
  iter <- 0
  Gprev <- G
  
  repeat {
    iter <- iter + 1
    # Eq 17
    U <- updateU(Xn, G, W, m = m, sigma = sigma)
    # Eq 19
    Gnew <- updateG(Xn, U, Gprev = G, W = W,
                           m = m, sigma = sigma,
                           alpha = alpha, xd = xd, q = q)
    # Eq 18
    Wnew <- updateW(Xn, U, Gnew,
                           m = m, sigma = sigma,
                           gamma = gamma, alpha = alpha)
    diff <- frNorm(Gnew - G)
    if (verbose) message("iter=", iter, "  ||Gcurr-Gprev||=", signif(diff, 6))
    Gprev <- G
    G <- Gnew
    W <- Wnew
    if (diff < eps || iter >= iM) break
  }
  list(U = U, G = G, W = W, xd = xd, iter = iter, khdi = init)
}

  