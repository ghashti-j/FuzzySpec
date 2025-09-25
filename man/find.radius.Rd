\name{find.radius}
\alias{find.radius}

\title{Adaptive Radius Selection by Natural Neighbours}

\description{
Implements a Natural Neighbour search to adaptively determine a neighbourhood
radius \eqn{r} from a general distance matrix. The algorithm increases \eqn{r}
until the number of points with zero in-degree in the \eqn{r}-nearest-neighbour
graph no longer descreases, and returns this radius \eqn{r}.
}

\usage{
find.radius(D)
}

\arguments{
  \item{D}{An \eqn{n \times n} numeric (symmetric) distance matrix.}
}

\details{
This procedure is adapated from the \emph{Natural Neighbor (NaN)}
algorithm by Zhu, Feng and Huang (2016). The algorithm works as follows:

\enumerate{
  \item For each integer \eqn{r}, build the directed \eqn{r}-nearest-neighbour
  (rNN) graph from \code{D}.
  \item Compute the in-degree counts, i.e., how many times each observation
  appears among others' first \eqn{r} neighbours.
  \item Track the number of zero in-degree points; if this number stops
  decreasing when \eqn{r} increases, the algorithm stops and returns the
  current \eqn{r}, interpreted as the natural neighbour radius.
}

This provides a parameter-free way to adaptively set the neighbourhood size.
}

\value{
A single integer \eqn{r}, the selected neighbourhood radius.
}

\seealso{
\code{\link{make.adjacency}}, \code{\link{gen.fuzzy}}, \code{\link{plot.fuzzy}}, \code{\link{rNN.dist}},
\code{\link{compute.sigma}}, \code{\link{compute.SNN}}, \code{\link{fuzzy.spectral.clustering}}
}

\references{
Zhu, Q., Feng, J., and J. Huang (2016). Natural neighbor: A self-adaptive neighborhood method without parameter K. \emph{Pattern recognition letters, 80}, 30-36.
}

\examples{
set.seed(1)
X <- matrix(rnorm(100), nrow = 20)
D <- as.matrix(dist(X))
r <- find.radius(D) # Estimate the natural neighbour radius
r
rNN.dist(D, r) # use selected r for rNN.dist
}

\keyword{nearest neighbors}
