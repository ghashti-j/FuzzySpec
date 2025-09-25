\name{compute.sigma}
\alias{compute.sigma}

\title{Compute Locally-Adaptive Scaling Parameters from a Distance Matrix}

\description{
Derives pointwise scale parameters \eqn{\sigma_i} from a distance matrix,
based on the \eqn{r}-th nearest neighbor distances. This is useful for
constructing the adaptive similarity graphs proposed by Zelnik-Manor and Perona (2004).
}

\usage{
compute.sigma(distance, r = NULL)
}

\arguments{
  \item{distance}{An \eqn{n \times n} numeric (symmetric) distance matrix. Required.}
  \item{r}{Integer valued neighbourhood radius. If \code{NULL}, it is estimated adaptively with \code{\link{find.radius}}.}
}

\details{
For each observation \eqn{i}, the function sorts the distances \code{distance[i, ]}, excludes the zero self-distance, and takes the \eqn{(r+1)} smallest value, where \eqn{\sigma_i} reflects the distance to the \eqn{r}-th nearest neighbour for observation \eqn{i}.

}

\value{
A list with components:
  \item{sigma}{A numeric vector of length \eqn{n}, containing local scale parameters.}
  \item{radius}{The neighborhood radius \eqn{r} used.}
}

\references{
  Ghashti, J. S., Hare, W., and J. R. J. Thompson (2025). Variable-weighted adjacency constructions for fuzzy spectral clustering. Submitted.

  Zelnik-Manor, L., and P. Perona (2004). Self-tuning spectral clustering. \emph{Advances in Neural Information Processing Systems, 17}.
}

\seealso{
\code{\link{make.adjacency}}, \code{\link{gen.fuzzy}}, \code{\link{plot.fuzzy}}, \code{\link{rNN.dist}}, \code{\link{find.radius}},
 \code{\link{compute.SNN}}, \code{\link{fuzzy.spectral.clustering}}
}

\examples{
set.seed(1)
X <- matrix(rnorm(50), nrow = 10)
D <- as.matrix(dist(X))

res <- compute.sigma(D) # automatically determine r
res$sigma
res$radius

res2 <- compute.sigma(D, r = 3) # user-specified r
res2$sigma
}

\keyword{nearest neighbours}
\keyword{adaptive similarity}
\keyword{kernel}
