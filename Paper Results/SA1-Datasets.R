library(ggplot2)
library(mvtnorm)
library(viridis)

#########################################
# spirals data
genSpiralData <- function(nPts = 600, seed = NULL, noise = 0.1, fuzzy = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  n_spiral <- floor(nPts / 3)
  nPtsRemain <- nPts - 3 * n_spiral
  n1 <- n_spiral + (nPtsRemain > 0)
  n2 <- n_spiral + (nPtsRemain > 1)
  n3 <- n_spiral
  t_vals <- list()
  
  # Spiral 1
  t1 <- seq(0, pi, length.out = n1)
  t_vals[[1]] <- t1
  r1 <- 0.5 + 0.8 * t1
  x1_clean <- r1 * cos(t1)
  y1_clean <- r1 * sin(t1)
  x1 <- x1_clean + rnorm(n1, 0, noise)
  y1 <- y1_clean + rnorm(n1, 0, noise)
  spiral1 <- cbind(x1, y1)
  
  # Spiral 2 
  t2 <- seq(0, pi, length.out = n2)
  t_vals[[2]] <- t2
  r2 <- 0.5 + 0.8 * t2
  angle_offset2 <- 2 * pi / 3
  x2_clean <- r2 * cos(t2 + angle_offset2)
  y2_clean <- r2 * sin(t2 + angle_offset2)
  x2 <- x2_clean + rnorm(n2, 0, noise)
  y2 <- y2_clean + rnorm(n2, 0, noise)
  spiral2 <- cbind(x2, y2)
  
  # Spiral 3 
  t3 <- seq(0, pi, length.out = n3)
  t_vals[[3]] <- t3
  r3 <- 0.5 + 0.8 * t3
  angle_offset3 <- 4 * pi / 3
  x3_clean <- r3 * cos(t3 + angle_offset3)
  y3_clean <- r3 * sin(t3 + angle_offset3)
  x3 <- x3_clean + rnorm(n3, 0, noise)
  y3 <- y3_clean + rnorm(n3, 0, noise)
  spiral3 <- cbind(x3, y3)
  
  X <- rbind(spiral1, spiral2, spiral3)
  yTrue <- c(rep(1, n1), rep(2, n2), rep(3, n3))
  
  if (!fuzzy) {
    uTrue <- matrix(0, nrow = nrow(X), ncol = 3)
    for (i in 1:nrow(X)) {
      uTrue[i, yTrue[i]] <- 1
    }
  } else {
    nTotal <- nrow(X)
    uFuzzy <- matrix(0, nTotal, 3)
    for (i in 1:nTotal) {
      point <- X[i, ]
      distances <- numeric(3)
      for (s in 1:3) {
        tDense <- seq(0, pi, length.out = 100)
        rDense <- 0.5 + 0.8 * tDense
        if (s == 1) {
          xSpiral <- rDense * cos(tDense)
          ySpiral <- rDense * sin(tDense)
        } else if (s == 2) {
          xSpiral <- rDense * cos(tDense + 2*pi/3)
          ySpiral <- rDense * sin(tDense + 2*pi/3)
        } else {
          xSpiral <- rDense * cos(tDense + 4*pi/3)
          ySpiral <- rDense * sin(tDense + 4*pi/3)
        }
        
        dists <- sqrt((xSpiral - point[1])^2 + (ySpiral - point[2])^2)
        distances[s] <- min(dists)
      }
      
      bandwidth <- 0.5
      for (s in 1:3) uFuzzy[i, s] <- exp(-distances[s]^2 / (2 * bandwidth^2))
      uFuzzy[i, yTrue[i]] <- uFuzzy[i, yTrue[i]] * 2
      distCentre <- sqrt(point[1]^2 + point[2]^2)
      if (distCentre < 1) {
        ctrWgt <- exp(-distCentre)
        uFuzzy[i, ] <- uFuzzy[i, ] * (1 - 0.5 * ctrWgt) + 
          0.5 * ctrWgt * rep(1/3, 3)
      }
      total <- sum(uFuzzy[i, ])
      if (total > 0) uFuzzy[i, ] <- uFuzzy[i, ] / total else uFuzzy[i, yTrue[i]] <- 1
    }
    uTrue <- uFuzzy
  }
  
  return(list(X = X, uTrue = uTrue, yTrue = yTrue, k = 3))
}

#####################################################
# gaussian data
genGaussianData <- function(nPts = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  centres <- matrix(c(-2, 0, 2, 0, 0, 3), ncol = 2, byrow = TRUE)
  cov1 <- matrix(c(1, 0.3, 0.3, 1), 2)
  cov2 <- matrix(c(1, -0.3, -0.3, 1), 2)
  cov3 <- matrix(c(0.8, 0, 0, 0.8), 2)

  n1 <- n2 <- n3 <- floor(nPts/3)
  X1 <- mvtnorm::rmvnorm(n1, centres[1,], cov1)
  X2 <- mvtnorm::rmvnorm(n2, centres[2,], cov2)
  X3 <- mvtnorm::rmvnorm(n3, centres[3,], cov3)
  
  X <- rbind(X1, X2, X3)
  d1 <- mvtnorm::dmvnorm(X, centres[1,], cov1)
  d2 <- mvtnorm::dmvnorm(X, centres[2,], cov2)
  d3 <- mvtnorm::dmvnorm(X, centres[3,], cov3)
  
  PIvec <- c(n1, n2, n3) / nPts
  wd1 <- PIvec[1] * d1
  wd2 <- PIvec[2] * d2
  wd3 <- PIvec[3] * d3
  
  densTotal <- wd1 + wd2 + wd3
  uTrue <- cbind(wd1/densTotal, wd2/densTotal, wd3/densTotal)
  yTrue <- apply(uTrue, 1, which.max)
  return(list(X = X, uTrue = uTrue, yTrue = yTrue, k = 3))
}

##############################################
# pie wedge dataset
genPieData <- function(nPts = 800, seed = NULL, innerRadius = 1, 
                                      outerRadius = 4, gapAngle = 0.15) {
  if (!is.null(seed)) set.seed(seed)
  nWedge <- 8
  nPtsWedge <- floor(nPts / nWedge)
  nPtsRemain <- nPts - nWedge * nPtsWedge
  gapTotal <- gapAngle * nWedge
  wedgeAngle <- (2 * pi - gapTotal) / nWedge
  X <- matrix(0, nrow = 0, ncol = 2)
  yTrue <- c()

  for (i in 1:nWedge) {
    nCurr <- nPtsWedge + (i <= nPtsRemain)
    startAngle <- (i - 1) * (wedgeAngle + gapAngle)
    endAngle <- startAngle + wedgeAngle
    angles <- runif(nCurr, startAngle, endAngle)
    radii <- runif(nCurr, innerRadius, outerRadius)
    x <- radii * cos(angles) + rnorm(nCurr, 0, 0.02)
    y <- radii * sin(angles) + rnorm(nCurr, 0, 0.02)
    wedgePts <- cbind(x, y)
    X <- rbind(X, wedgePts)
    yTrue <- c(yTrue, rep(i, nCurr))
  }
  
  # fuzzy memberships based on angular distance
  nTotal <- nrow(X)
  uFuzzy <- matrix(0, nTotal, nWedge)
  for (i in 1:nTotal) {
    point <- X[i, ]
    r <- sqrt(point[1]^2 + point[2]^2)
    theta <- atan2(point[2], point[1])
    
    if (theta < 0) theta <- theta + 2 * pi
    
    for (j in 1:nWedge) {
      ctrAngle <- (j - 1) * (wedgeAngle + gapAngle) + wedgeAngle / 2
      distAngle <- abs(theta - ctrAngle)
      distAngle <- min(distAngle, 2 * pi - distAngle)
      # membership based on angular distance
      membAngle <- exp(-distAngle^2 / (0.3^2))
      radPenalty <- 1
      if (r < innerRadius + 0.5) {
        radPenalty <- (r - innerRadius) / 0.5
      } else if (r > outerRadius - 0.5) {
        radPenalty <- (outerRadius - r) / 0.5
      }
      radPenalty <- max(0, min(1, radPenalty))
      
      uFuzzy[i, j] <- membAngle * radPenalty
    }
    
    uFuzzy[i, yTrue[i]] <- uFuzzy[i, yTrue[i]] * 3
    total <- sum(uFuzzy[i, ])
    
    if (total > 0) uFuzzy[i, ] <- uFuzzy[i, ] / total else uFuzzy[i, yTrue[i]] <- 1
  }
  return(list(X = X, uTrue = uFuzzy, yTrue = yTrue, k = nWedge))
}

#############################################
# Rings Dataset
genRingData <- function(nPts = 600, seed = NULL, noise = 0.05) {
  if (!is.null(seed)) set.seed(seed)
  nRings <- 3
  radii <- c(1, 2.5, 4)
  widths <- c(0.3, 0.4, 0.5)
  nPtsRings <- floor(nPts / nRings)
  nPtsRemain <- nPts - nRings * nPtsRings
  X <- matrix(0, nrow = 0, ncol = 2)
  yTrue <- c()

  for (i in 1:nRings) {
    nCurr <- nPtsRings + (i <= nPtsRemain)
    angles <- runif(nCurr, 0, 2*pi)
    rBase <- radii[i]
    rVar <- runif(nCurr, -widths[i]/2, widths[i]/2)
    r <- rBase + rVar
    x <- r * cos(angles) + rnorm(nCurr, 0, noise)
    y <- r * sin(angles) + rnorm(nCurr, 0, noise)
    
    ptsRings <- cbind(x, y)
    X <- rbind(X, ptsRings)
    yTrue <- c(yTrue, rep(i, nCurr))
  }
  
  # fuzzy memberships based on radial distance
  nTotal <- nrow(X)
  uFuzzy <- matrix(0, nTotal, nRings)
  
  for (i in 1:nTotal) {
    point <- X[i, ]
    r <- sqrt(point[1]^2 + point[2]^2)
    for (j in 1:nRings) {
      distRings <- abs(r - radii[j])
      uFuzzy[i, j] <- exp(-distRings^2 / (2 * widths[j]^2))
    }
    total <- sum(uFuzzy[i, ])
    if (total > 0) uFuzzy[i, ] <- uFuzzy[i, ] / total else uFuzzy[i, yTrue[i]] <- 1
  }
  return(list(X = X, uTrue = uFuzzy, yTrue = yTrue, k = nRings))
}

######################################################
# Worms Dataset
genWormData <- function(nPts = 600, seed = NULL, noise = 0.1) {
  if (!is.null(seed)) set.seed(seed)
  
  nWorms <- 4
  nPtsWorms <- floor(nPts / nWorms)
  nPtsRemain <- nPts - nWorms * nPtsWorms
  X <- matrix(0, nrow = 0, ncol = 2)
  yTrue <- c()

  paramWorms <- list(
    list(amp = 1.5, freq = 1, phase = 0, yOff = 4),
    list(amp = 1.2, freq = 0.8, phase = pi/2, yOff = 2),
    list(amp = 1.8, freq = 1.2, phase = pi, yOff = -1),
    list(amp = 1.3, freq = 0.9, phase = 3*pi/2, yOff = -4)
  )
  
  for (i in 1:nWorms) {
    nCurr <- nPtsWorms + (i <= nPtsRemain)
    t <- seq(0, 2*pi, length.out = nCurr)
    params <- paramWorms[[i]]
    xraw <- t - pi  
    x <- 2*xraw
    y <- params$amp * sin(params$freq * t + params$phase) + params$yOff
    width <- 0.3
    perpOff <- runif(nCurr, -width/2, width/2)
    dx <- c(diff(x), x[nCurr] - x[nCurr-1])
    dy <- c(diff(y), y[nCurr] - y[nCurr-1])
    norm <- sqrt(dx^2 + dy^2)
    perpX <- -dy / norm * perpOff
    perpY <- dx / norm * perpOff
    x <- x + perpX + rnorm(nCurr, 0, noise)
    y <- y + perpY + rnorm(nCurr, 0, noise)
    
    ptsWorms <- cbind(x, y)
    X <- rbind(X, ptsWorms)
    yTrue <- c(yTrue, rep(i, nCurr))
  }
  
  # fuzzy memberships based on distance to worm centerlines
  nTotal <- nrow(X)
  uFuzzy <- matrix(0, nTotal, nWorms)
  
  for (i in 1:nTotal) {
    point <- X[i, ]
    for (j in 1:nWorms) {
      params <- paramWorms[[j]]
      ttest <- seq(0, 2*pi, length.out = 100)
      xtest <- ttest - pi
      ytest <- params$amp * sin(params$freq * ttest + params$phase) + params$yOff
      distances <- sqrt((xtest - point[1])^2 + (ytest - point[2])^2)
      mdist <- min(distances)
      uFuzzy[i, j] <- exp(-mdist^2 / (2 * 0.5^2))
    }
    
    total <- sum(uFuzzy[i, ])
    if (total > 0) uFuzzy[i, ] <- uFuzzy[i, ] / total else uFuzzy[i, yTrue[i]] <- 1
  }
  return(list(X = X, uTrue = uFuzzy, yTrue = yTrue, k = nWorms))
}


###################################################
# hyperbolas dataset

genHypeDataHard <- function(nPts = 500, seed = NULL, noise = 0.1) {
  if (!is.null(seed)) set.seed(seed)
  nballWeight <- floor(nPts * 0.2)
  n_hyp <- floor(nPts * 0.2)
  nPtsRemain <- nPts - nballWeight - 4 * n_hyp
  nballWeight <- nballWeight + nPtsRemain
  ball <- mvtnorm::rmvnorm(nballWeight, mean = c(0, 0), sigma = 0.2 * diag(2))
  a <- 0.5
  b <- 0.5
  offset <- 1.2
  # Hyperbola 1
  t1 <- seq(-2, 2, length.out = n_hyp)
  x1 <- offset + a * cosh(t1)
  y1 <- b * sinh(t1)
  hyp1 <- cbind(x1, y1) + matrix(rnorm(n_hyp * 2, 0, noise), ncol = 2)
  # Hyperbola 2
  t2 <- seq(-2, 2, length.out = n_hyp)
  x2 <- -offset - a * cosh(t2)
  y2 <- b * sinh(t2)
  hyp2 <- cbind(x2, y2) + matrix(rnorm(n_hyp * 2, 0, noise), ncol = 2)
  # Hyperbola 3
  t3 <- seq(-2, 2, length.out = n_hyp)
  x3 <- a * sinh(t3)
  y3 <- offset + b * cosh(t3)
  hyp3 <- cbind(x3, y3) + matrix(rnorm(n_hyp * 2, 0, noise), ncol = 2)
  # Hyperbola 4 
  t4 <- seq(-2, 2, length.out = n_hyp)
  x4 <- a * sinh(t4)
  y4 <- -offset - b * cosh(t4)
  hyp4 <- cbind(x4, y4) + matrix(rnorm(n_hyp * 2, 0, noise), ncol = 2)
  X <- rbind(ball, hyp1, hyp2, hyp3, hyp4)
  yTrue <- c(rep(1, nballWeight), rep(2, n_hyp), rep(3, n_hyp), rep(4, n_hyp), rep(5, n_hyp))
  uTrue <- matrix(0, nrow = nrow(X), ncol = 5)
  for (i in 1:nrow(X)) uTrue[i, yTrue[i]] <- 1
  return(list(X = X, uTrue = uTrue, yTrue = yTrue, k = 5))
}


genHypeData <- function(nPts = 500, seed = NULL, noise = 0.1) {
  if (!is.null(seed)) set.seed(seed)
  data <- genHypeDataHard(nPts, seed = seed, noise = noise)
  X <- data$X
  yTrue <- data$yTrue
  nTotal <- nrow(X)
  uFuzzy <- matrix(0, nTotal, 5)
  ballCentre <- c(0, 0)
  ballCov <- 0.2 * diag(2)
  a <- 0.5
  b <- 0.5
  offset <- 1.2
  
  for (i in 1:nTotal) {
    point <- X[i, ]
    ballDist <- mvtnorm::dmvnorm(matrix(point, nrow = 1), mean = ballCentre, sigma = ballCov)
    tRange <- seq(-2, 2, length.out = 50)
    
    curve1 <- cbind(offset + a * cosh(tRange), b * sinh(tRange))
    d1 <- min(sqrt(rowSums((curve1 - matrix(point, nrow = 50, ncol = 2, byrow = TRUE))^2)))
    curve2 <- cbind(-offset - a * cosh(tRange), b * sinh(tRange))
    d2 <- min(sqrt(rowSums((curve2 - matrix(point, nrow = 50, ncol = 2, byrow = TRUE))^2)))
    curve3 <- cbind(a * sinh(tRange), offset + b * cosh(tRange))
    d3 <- min(sqrt(rowSums((curve3 - matrix(point, nrow = 50, ncol = 2, byrow = TRUE))^2)))
    curve4 <- cbind(a * sinh(tRange), -offset - b * cosh(tRange))
    d4 <- min(sqrt(rowSums((curve4 - matrix(point, nrow = 50, ncol = 2, byrow = TRUE))^2)))
    
    bandwidth <- 0.5
    wballWeight <- ballDist * 50 
    w1 <- exp(-d1^2 / (2 * bandwidth^2))
    w2 <- exp(-d2^2 / (2 * bandwidth^2))
    w3 <- exp(-d3^2 / (2 * bandwidth^2))
    w4 <- exp(-d4^2 / (2 * bandwidth^2))

    weights <- c(wballWeight, w1, w2, w3, w4)
    weights[yTrue[i]] <- weights[yTrue[i]] * 3
    
    totalWeight <- sum(weights)
    if (totalWeight > 0) uFuzzy[i, ] <- weights / totalWeight else uFuzzy[i, yTrue[i]] <- 1
  }
  return(list(X = X, uTrue = uFuzzy, yTrue = yTrue, k = 5))
}
