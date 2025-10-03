## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 6,
  fig.height = 5,
  message = FALSE,
  warning = FALSE
)
library(mclust)
library(fclust)
library(ggplot2)
library(patchwork)
library(mvtnorm)
library(stats)
library(knitr)
library(np)
library(MASS)
library(rmarkdown)

## ----eval = FALSE-------------------------------------------------------------
# library(devtools)
# install_github("ghashti-j/FuzzySpec")
# library(FuzzySpec)

## ----fig.align='center'-------------------------------------------------------
set.seed(1)
data <- FuzzySpec::gen.fuzzy(n = 300, dataset = "spirals", noise = 0.15) # data generation
FuzzySpec::plot.fuzzy(data, plotFuzzy = TRUE, colorCluster = TRUE) # plot data generating process

## ----message = FALSE----------------------------------------------------------
W <- FuzzySpec::make.adjacency(
  data = data$X,
  method = "vw",           # variable-weighted distances
  isLocWeighted = TRUE,    # Locally-adaptive scaling
  scale = FALSE            # scaling not required for kernel methods
)

## -----------------------------------------------------------------------------
res <- FuzzySpec::fuzzy.spectral.clustering(
  W = W, k = 3, m = 1.5, method = "CM"           
)
res$u[1:5,]

## -----------------------------------------------------------------------------
acc <- FuzzySpec::clustering.accuracy(data$y, res$cluster)
cat("Clustering accuracy:", round(acc, 3), "\n")

## -----------------------------------------------------------------------------
far <- FuzzySpec::fari(data$U, res$u)
cat("FARI:", round(far, 3), "\n")

## ----fig.align='center'-------------------------------------------------------
resDF <- list(
  X = data$X, U = res$u, y = factor(res$cluster), k = 3
)
FuzzySpec::plot.fuzzy(resDF, plotFuzzy = TRUE, colorCluster = TRUE)

