gen.fuzzy <- function(n = 500,
                      dataset = c("gaussian", "hyperbolas", "spirals",
                                          "wedges", "rings", "worms", "random"),
                      k = NULL,
                      noise = 0.1,
                      covType = c("spherical", "diagonal", "rotated", "correlated"),
                      seed = NULL) {

  dataset <- match.arg(dataset)
  covType <- match.arg(covType)
  if (!is.null(seed)) set.seed(seed)
  genGaussian <- function(n, noise) {
    centres <- matrix(c(-2, 0, 2, 0, 0, 3), ncol = 2, byrow = TRUE)
    cov1 <- matrix(c(1, 0.3, 0.3, 1), 2)
    cov2 <- matrix(c(1, -0.3, -0.3, 1), 2)
    cov3 <- matrix(c(0.8, 0, 0, 0.8), 2)
    covs <- list(cov1, cov2, cov3)
    nPerCl <- floor(n/3); nRem <- n-3*nPerCl; sizes <- rep(nPerCl, 3)
    if(nRem > 0) sizes[1:nRem] <- sizes[1:nRem]+1
    X <- matrix(0, n, 2); yTrue <- numeric(n); idx <- 1
    for(i in 1:3) {
      pts <- mvtnorm::rmvnorm(sizes[i], centres[i,], covs[[i]])
      X[idx:(idx+sizes[i]-1), ] <- pts
      yTrue[idx:(idx+sizes[i]-1)] <- i
      idx <- idx+sizes[i]
    }
    U <- matrix(0, n, 3); piVec <- sizes/n
    for(i in 1:3) U[,i] <- piVec[i]*mvtnorm::dmvnorm(X, centres[i,], covs[[i]])
    U <- U/rowSums(U)
    list(X = X, U = U, y = yTrue, k = 3)
  }

  genHyperbolas <- function(n, noise) {
    nBall <- floor(n*0.2); nHyp <- floor(n*0.2)
    nRem <- n-nBall-4*nHyp; nBall <- nBall+nRem
    ball <- mvtnorm::rmvnorm(nBall, c(0, 0), 0.2*diag(2))
    a <- 0.5; b <- 0.5; offS <- 1.2
    t <- seq(-2, 2, length.out = nHyp)

    hyp1 <- cbind(offS+a*cosh(t), b*sinh(t))+
      matrix(rnorm(nHyp*2, 0, noise), ncol = 2)
    hyp2 <- cbind(-offS-a*cosh(t), b*sinh(t))+
      matrix(rnorm(nHyp*2, 0, noise), ncol = 2)
    hyp3 <- cbind(a*sinh(t), offS+b*cosh(t))+
      matrix(rnorm(nHyp*2, 0, noise), ncol = 2)
    hyp4 <- cbind(a*sinh(t), -offS-b*cosh(t))+
      matrix(rnorm(nHyp*2, 0, noise), ncol = 2)

    X <- rbind(ball, hyp1, hyp2, hyp3, hyp4)
    yTrue <- c(rep(1, nBall), rep(2, nHyp), rep(3, nHyp),
               rep(4, nHyp), rep(5, nHyp))
    U <- calcH.Mem(X, yTrue, a, b, offS)
    list(X = X, U = U, y = yTrue, k = 5)
  }

  calcH.Mem <- function(X, yTrue, a, b, offS) {
    n <- nrow(X)
    U <- matrix(0, n, 5)
    bw.bw <- 0.5
    for(i in 1:n) {
      pt <- X[i,]
      wBall <- mvtnorm::dmvnorm(matrix(pt, nrow = 1), c(0, 0), 0.2*diag(2))*50
      tRange <- seq(-2, 2, length.out = 50)

      d1 <- calc.mDist(pt, cbind(offS+a*cosh(tRange), b*sinh(tRange)))
      d2 <- calc.mDist(pt, cbind(-offS-a*cosh(tRange), b*sinh(tRange)))
      d3 <- calc.mDist(pt, cbind(a*sinh(tRange), offS+b*cosh(tRange)))
      d4 <- calc.mDist(pt, cbind(a*sinh(tRange), -offS-b*cosh(tRange)))
      wgts <- c(wBall,
                   exp(-d1^2/(2*bw.bw^2)),
                   exp(-d2^2/(2*bw.bw^2)),
                   exp(-d3^2/(2*bw.bw^2)),
                   exp(-d4^2/(2*bw.bw^2)))
      wgts[yTrue[i]] <- wgts[yTrue[i]]*3
      if(sum(wgts) > 0) U[i,] <- wgts/sum(wgts) else U[i, yTrue[i]] <- 1
    }
    U
  }

  calc.mDist <- function(pt, curve) {
    dists <- sqrt(rowSums((curve-matrix(pt, nrow(curve), 2, byrow = TRUE))^2))
    min(dists)
  }

  genSpirals <- function(n, noise) {
    nSpiral <- floor(n/3); nRem <- n-3*nSpiral
    n1 <- nSpiral+(nRem > 0); n2 <- nSpiral+(nRem > 1); n3 <- nSpiral

    t1 <- seq(0, pi, length.out = n1)
    r1 <- 0.5+0.8*t1
    x1 <- r1*cos(t1)+rnorm(n1, 0, noise)
    y1 <- r1*sin(t1)+rnorm(n1, 0, noise)

    t2 <- seq(0, pi, length.out = n2)
    r2 <- 0.5+0.8*t2
    x2 <- r2*cos(t2+2*pi/3)+rnorm(n2, 0, noise)
    y2 <- r2*sin(t2+2*pi/3)+rnorm(n2, 0, noise)

    t3 <- seq(0, pi, length.out = n3)
    r3 <- 0.5+0.8*t3
    x3 <- r3*cos(t3+4*pi/3)+rnorm(n3, 0, noise)
    y3 <- r3*sin(t3+4*pi/3)+rnorm(n3, 0, noise)

    X <- rbind(cbind(x1, y1), cbind(x2, y2), cbind(x3, y3))
    yTrue <- c(rep(1, n1), rep(2, n2), rep(3, n3))
    U <- calcS.Mem(X, yTrue)
    list(X = X, U = U, y = yTrue, k = 3)
  }

  calcS.Mem <- function(X, yTrue) {
    n <- nrow(X)
    U <- matrix(0, n, 3)
    bw.bw <- 0.5
    for(i in 1:n) {
      pt <- X[i,]
      dists <- numeric(3)
      for(s in 1:3) {
        tDense <- seq(0, pi, length.out = 100)
        rDense <- 0.5+0.8*tDense
        angle <- if(s == 1) 0 else if(s == 2) 2*pi/3 else 4*pi/3
        xSpiral <- rDense*cos(tDense+angle)
        ySpiral <- rDense*sin(tDense+angle)
        dists[s] <- calc.mDist(pt, cbind(xSpiral, ySpiral))
      }
      for(s in 1:3) U[i,s] <- exp(-dists[s]^2/(2*bw.bw^2))
      U[i, yTrue[i]] <- U[i, yTrue[i]]*2
      distCenter <- sqrt(sum(pt^2))
      if(distCenter < 1) {
        centerWeight <- exp(-distCenter)
        U[i,] <- U[i,]*(1-0.5*centerWeight)+0.5*centerWeight/3
      }
      if(sum(U[i,]) > 0) U[i,] <- U[i,]/sum(U[i,]) else U[i, yTrue[i]] <- 1
    }
    U
  }

  genWedges <- function(n) {
    nWedges <- 8; nPerWedge <- floor(n/nWedges)
    nRem <- n-nWedges*nPerWedge; gapAng <- 0.15
    totalGap <- gapAng*nWedges
    wAngle <- (2*pi-totalGap)/nWedges
    inRad <- 1; outRad <- 4
    X <- matrix(0, 0, 2)
    yTrue <- c()

    for(i in 1:nWedges) {
      nCurr <- nPerWedge+(i <= nRem)
      sAngle <- (i-1)*(wAngle+gapAng)
      eAngle <- sAngle+wAngle
      angles <- runif(nCurr, sAngle, eAngle)
      radii <- runif(nCurr, inRad, outRad)
      x <- radii*cos(angles)+rnorm(nCurr, 0, 0.02)
      y <- radii*sin(angles)+rnorm(nCurr, 0, 0.02)
      X <- rbind(X, cbind(x, y))
      yTrue <- c(yTrue, rep(i, nCurr))
    }

    U <- calcW.Memb(X, yTrue, nWedges, wAngle, gapAng, inRad, outRad)
    list(X = X, U = U, y = yTrue, k = nWedges)
  }

  calcW.Memb <- function(X, yTrue, nWedges, wAngle, gapAng, inRad, outRad) {
    n <- nrow(X); U <- matrix(0, n, nWedges)
    for(i in 1:n) {
      pt <- X[i,]
      r <- sqrt(sum(pt^2))
      theta <- atan2(pt[2], pt[1])
      if(theta < 0) theta <- theta+2*pi
      for(j in 1:nWedges) {
        centerAngle <- (j-1)*(wAngle+gapAng)+wAngle/2
        distAngle <- abs(theta-centerAngle)
        distAngle <- min(distAngle, 2*pi-distAngle)
        angMemb <- exp(-distAngle^2/0.3^2)
        radPenal <- 1
        if(r < inRad+0.5) {
          radPenal <- (r-inRad)/0.5
        } else if(r > outRad-0.5) {
          radPenal <- (outRad-r)/0.5
        }
        radPenal <- max(0, min(1, radPenal))
        U[i,j] <- angMemb*radPenal
      }
      U[i, yTrue[i]] <- U[i, yTrue[i]]*3
      if(sum(U[i,]) > 0) U[i,] <- U[i,]/sum(U[i,]) else U[i, yTrue[i]] <- 1
    }
    U
  }

  genRings <- function(n, noise) {
    nRings <- 3; radii <- c(1, 2.5, 4); widths <- c(0.3, 0.4, 0.5)
    numRing <- floor(n/nRings); nRem <- n-nRings*numRing
    X <- matrix(0, 0, 2); yTrue <- c()
    for(i in 1:nRings) {
      nCurr <- numRing+(i <= nRem)
      angles <- runif(nCurr, 0, 2*pi)
      rVar <- runif(nCurr, -widths[i]/2, widths[i]/2)
      r <- radii[i]+rVar
      x <- r*cos(angles)+rnorm(nCurr, 0, noise)
      y <- r*sin(angles)+rnorm(nCurr, 0, noise)
      X <- rbind(X, cbind(x, y))
      yTrue <- c(yTrue, rep(i, nCurr))
    }
    U <- matrix(0, nrow(X), nRings)
    for(i in 1:nrow(X)) {
      pt <- X[i,]
      r <- sqrt(sum(pt^2))
      for(j in 1:nRings) {
        distFromRing <- abs(r-radii[j])
        U[i,j] <- exp(-distFromRing^2/(2*widths[j]^2))
      }
      if(sum(U[i,]) > 0) U[i,] <- U[i,]/sum(U[i,]) else U[i, yTrue[i]] <- 1
    }
    list(X = X, U = U, y = yTrue, k = nRings)
  }

  genWorms <- function(n, noise) {
    nWorms <- 4; nPerWorm <- floor(n/nWorms); nRem <- n-nWorms*nPerWorm
    X <- matrix(0, 0, 2); yTrue <- c()
    wormPars <- list(
      list(amp = 1.5, freq = 1, phase = 0, yoffS = 4),
      list(amp = 1.2, freq = 0.8, phase = pi/2, yoffS = 2),
      list(amp = 1.8, freq = 1.2, phase = pi, yoffS = -1),
      list(amp = 1.3, freq = 0.9, phase = 3*pi/2, yoffS = -4)
    )

    for(i in 1:nWorms) {
      nCurr <- nPerWorm+(i <= nRem)
      t <- seq(0, 2*pi, length.out = nCurr)
      params <- wormPars[[i]]
      x <- 2*(t-pi)
      y <- params$amp*sin(params$freq*t+params$phase)+params$yoffS
      width <- 0.3
      perpoffS <- runif(nCurr, -width/2, width/2)
      dx <- c(diff(x), x[nCurr]-x[nCurr-1])
      dy <- c(diff(y), y[nCurr]-y[nCurr-1])
      norm <- sqrt(dx^2+dy^2)
      perpX <- -dy/norm*perpoffS
      perpY <- dx/norm*perpoffS
      x <- x+perpX+rnorm(nCurr, 0, noise)
      y <- y+perpY+rnorm(nCurr, 0, noise)
      X <- rbind(X, cbind(x, y))
      yTrue <- c(yTrue, rep(i, nCurr))
    }

    U <- calcWorm.Mem(X, yTrue, wormPars)
    list(X = X, U = U, y = yTrue, k = nWorms)
  }

  calcWorm.Mem <- function(X, yTrue, wormPars) {
    n <- nrow(X); nWorms <- length(wormPars)
    U <- matrix(0, n, nWorms)
    for(i in 1:n) {
      pt <- X[i,]
      for(j in 1:nWorms) {
        params <- wormPars[[j]]
        tTest <- seq(0, 2*pi, length.out = 100)
        xTest <- 2*(tTest-pi)
        yTest <- params$amp*sin(params$freq*tTest+params$phase)+params$yoffS
        minDist <- calc.mDist(pt, cbind(xTest, yTest))
        U[i,j] <- exp(-minDist^2/(2*0.5^2))
      }
      if(sum(U[i,]) > 0) U[i,] <- U[i,]/sum(U[i,]) else U[i, yTrue[i]] <- 1
    }
    U
  }

  genRandom <- function(n, k, covType) {
    if(is.null(k)) k <- 20
    covMatrix <- switch(covType,
                        "spherical" = diag(2),
                        "diagonal" = diag(c(2, 0.5)),
                        "rotated" = {
                          rot <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)), 2, 2)
                          rot %*% diag(c(2, 0.5)) %*% t(rot)
                        },
                        "correlated" = matrix(c(1, 0.7, 0.7, 1), 2, 2)
    )

    centres <- matrix(runif(k*2, 0, 30), ncol = 2)
    clusSz <- sample(20:80, k, replace = TRUE)
    scFac <- n/sum(clusSz)
    clusSz <- round(clusSz*scFac)
    diff <- n-sum(clusSz)
    if(diff > 0) {
      addTo <- sample(1:k, diff, replace = TRUE)
      for(i in addTo) clusSz[i] <- clusSz[i]+1
    } else if(diff < 0) {
      while(diff < 0) {
        largest <- which.max(clusSz)
        clusSz[largest] <- clusSz[largest]-1
        diff <- diff+1
      }
    }

    datList <- list()
    for(i in 1:k) datList[[i]] <- mvtnorm::rmvnorm(clusSz[i], centres[i,], covMatrix)

    X <- do.call(rbind, datList); yTrue <- unlist(mapply(rep, 1:k, clusSz))
    U <- matrix(0, n, k); piVec <- clusSz/n

    for(j in 1:k) U[,j] <- piVec[j]*mvtnorm::dmvnorm(X, centres[j,], covMatrix)

    tDens <- rowSums(U)
    zDens <- tDens == 0

    if(any(zDens)) {
      for(i in which(zDens)) {
        dists <- sqrt(rowSums((centres-matrix(X[i,], k, 2, byrow = TRUE))^2))
        NN <- which.min(dists)
        U[i, NN] <- 1
      }
    }

    U[!zDens,] <- U[!zDens,]/tDens[!zDens]
    list(X = X, U = U, y = yTrue, k = k, centres = centres,
         clusSz = clusSz, covMatrix = covMatrix)
  }
  result <- switch(dataset,
                   "gaussian" = genGaussian(n, noise),
                   "hyperbolas" = genHyperbolas(n, noise),
                   "spirals" = genSpirals(n, noise),
                   "wedges" = genWedges(n),
                   "rings" = genRings(n, noise),
                   "worms" = genWorms(n, noise),
                   "random" = genRandom(n, k, covType)
  )
  return(result)
}



