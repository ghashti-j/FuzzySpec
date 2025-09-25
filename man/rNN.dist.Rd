\name{rNN.dist}
\alias{rNN.dist}

\title{Compute r-Nearest Neighbours from a Distance Matrix}

\description{
Given a symmetric distance matrix, returns the indices of the \eqn{r} nearest neighbours for each observation.
}

\usage{
rNN.dist(D, r)
}

\arguments{
  \item{D}{An \eqn{n \times n} numeric (symmetric) distance matrix.}
  \item{r}{Integer indicating the number of nearest neighbours to extract for each observation.}
}

\details{
For each row \eqn{i} of \code{D}, the function orders the distances \eqn{D[i, \cdot]},
excludes the self-distance, and returns the indices of the first \code{r} smallest distances.
This provides the indices of the \code{r} nearest neighbours of observation \eqn{i}.
}

\value{
An \eqn{n \times r} integer matrix, where row \eqn{i} contains the indices of the \code{r} nearest neighbours of observation \eqn{i}.
}

\references{
  Ghashti, J. S., Hare, W., and J. R. J. Thompson (2025). Variable-weighted adjacency constructions for fuzzy spectral clustering. Submitted.
}

\seealso{
\code{\link{make.adjacency}}, \code{\link{gen.fuzzy}}, \code{\link{plot.fuzzy}}, \code{\link{find.radius}},
\code{\link{compute.sigma}}, \code{\link{compute.SNN}}, \code{\link{fuzzy.spectral.clustering}}
}

\examples{
set.seed(1)
X <- matrix(rnorm(20), nrow = 5)
D <- as.matrix(dist(X))
rNN.dist(D, r = 2)  # find 2 nearest neighbours for each row
}

\keyword{nearest neighbors}
\keyword{distance}
\keyword{graph}

