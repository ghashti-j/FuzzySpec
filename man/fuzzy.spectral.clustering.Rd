\name{fuzzy.spectral.clustering}
\alias{fuzzy.spectral.clustering}

\title{Fuzzy Spectral Clustering with Normalized Eigenvectors}

\description{
Implementation of the FVIBES algorithm by Ghashti, Hare, and Thompson (2025). Performs spectral clustering on a similarity (adjacency) matrix and returns either fuzzy c-means memberships or Gaussian mixture posterior probabilities computed on the leading normalized eigenvectors.
}

\usage{
fuzzy.spectral.clustering(W = NULL,
                          k = NULL,
                          m = NULL,
                          method = "CM",
                          nstart = 10,
                          max.iter = 1000)
}

\arguments{
  \item{W}{A nonnegative \eqn{n \times n} similarity (adjacency) matrix. Diagonal entries are set to 0 internally. Required.}
  \item{k}{Integer number of clusters. Required.}
  \item{m}{Fuzzifier for c-means, only used when \code{method = "CM"}. When not provided, algorithm will set \code{m = 2}.}
  \item{method}{Clustering method applied to the spectral embedding with \code{"CM"} for fuzzy c-means with \pkg{fclust}, or \code{"GMM"} for Gaussian mixtures with \pkg{mclust}. Default is \code{"CM"}.}
  \item{nstart}{Number of random starts for \code{fclust::FKM} when \code{method = "CM"}.}
  \item{max.iter}{Maximum number of iterations for \code{fclust::FKM} when \code{method = "CM"}.}
}

\details{
Let \eqn{D} be the diagonal degree matrix with \eqn{D_{ii} = \sum_j W_{ij}}. The routine forms the symmetrically normalized similarity
\eqn{L = D^{-1/2} W D^{-1/2},} (Ng, Jordan, and Weiss, 2001) computes its top \eqn{k} eigenvectors, stacks them in \eqn{X \in \mathbb{R}^{n \times k}}, and row-normalizes to \eqn{Y} with
\eqn{Y_{i\cdot} = X_{i\cdot} / \|X_{i\cdot}\|_2}. Clustering is then performed in the rows of \eqn{Y}.

When \code{method = "CM"}, clustering uses c-means (Bezdek, 1981) with \code{fclust::FKM} on \eqn{Y} with fuzzifier \code{m}, number of starts \code{nstart}, and maximum iterations \code{max.iter}.
When \code{method = "GMM"}, clustering uses Gaussian mixture models (see McLachlan and Krishnan, 2008) with \code{mclust::Mclust} with \code{G = k} on \eqn{Y}.
}

\value{
A list with components:
  \item{cluster}{An integer vector of length \eqn{n} hard cluster labels.}
  \item{u}{An \eqn{n \times k} matrix of fuzzy cluster memberships: for \code{"CM"}, fuzzy c-means memberships \eqn{U}; for \code{"GMM"}, posterior probabilities \eqn{Z}.}
  \item{evecs}{The \eqn{n \times k} matrix \eqn{Y} of row-normalized leading eigenvectors, i.e., the spectral embedding.}
  \item{centers}{Cluster centers for the embedding matrix \eqn{Y}.}
}

\references{
  J.C. Bezdek (1981). \emph{Pattern Recognition with Fuzzy Objective Function Algorithms}. Plenum Press, New York.

  Ferraro, M.B., Giordani, P., and A. Serafini (2019). fclust: An R Package for Fuzzy Clustering. \emph{The R Journal, 11}.

  Ghashti, J. S., Hare, W., and J. R. J. Thompson (2025). Variable-weighted adjacency constructions for fuzzy spectral clustering. Submitted.

  McLachlan, G. and T. Krishnan (2008). \emph{The EM algorithm and extensions}, Second Edition. John Wiley & Sons.

  Ng, A., Jordan, M., and Y. Weiss (2001). On spectral clustering: Analysis and an algorithm. \emph{Advances in Neural Information Processing Systems, 14.}

  Scrucca, L., Fraley, C., Murphy, T.B., and A. E. Raftery (2023). \emph{Model-Based Clustering, Classification, and Density Estimation Using mclust in R}. Chapman and Hall/CRC.
}

\seealso{
\code{\link{make.adjacency}}, \code{\link{gen.fuzzy}}, \code{\link{plot.fuzzy}}, \code{\link{rNN.dist}}, \code{\link{find.radius}},
\code{\link{compute.sigma}}, \code{\link{compute.SNN}}, \code{\link[np]{npudensbw}}, \code{\link[fclust]{FKM}}, \code{\link[mclust]{Mclust}}
}

\examples{
set.seed(1)
d <- gen.fuzzy(n = 300,
               dataset = "spirals",
               noise = 0.18)

plot.fuzzy(d) # visualize data generating process

adj <- make.adjacency(data = d$X,
                      method = "vw",
                      isLocWeighted = TRUE,
                      isModWeighted = FALSE,
                      isSparse = FALSE,
                      ModMethod = NULL,
                      scale = FALSE,
                      sig = 1,
                      radius = NULL,
                      cv.method = "cv.ls") # vwla-id from paper

spectRes <- fuzzy.spectral.clustering(W = adj,
                                      k = 3,
                                      m = 1.5,
                                      method = "CM",
                                      nstart = 50,
                                      max.iter = 1000)


head(spectRes$u) # first 6 rows of U

plotDf <- list(
  X = d$X,
  y = factor(spectRes$cluster),
  U = spectRes$u,
  k = 3
)

plot.fuzzy(plotDf) # visualize results

clustering.accuracy(d$y, spectRes$cluster) # compare results
}

\keyword{clustering}
\keyword{spectral}
\keyword{fuzzy}
\keyword{graph}
\keyword{cmeans}
\keyword{fuzzy clustering}
