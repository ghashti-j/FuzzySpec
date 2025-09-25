\name{compute.SNN}
\alias{compute.SNN}

\title{Shared-Nearest-Neighbours (SNN) Similarity from a Similarity Matrix}

\description{
Builds a Shared-Nearest-Neighbours (SNN) similarity matrix from an input similarity matrix \code{similarity}. For each pair of observations, the SNN score is the fraction of shared indices among their top-\eqn{r} neighbour lists.
}

\usage{
compute.SNN(similarity, r)
}

\arguments{
  \item{similarity}{An \eqn{n \times n} numeric \emph{similarity} matrix. The diagonal is assumed to correspond to self-similarity and is ignored when forming neighboir lists.}
  \item{r}{Integer value number of nearest neighbours per observation used to compute SNN overlap.}
}

\details{
For each observation \eqn{i}, the function forms its neighbour set by ordering
\code{similarity[i, ]} in decreasing order, dropping \code{i} itself, and
retaining the first \eqn{r} indices. For a pair \eqn{(i, j)}, the SNN similarity
is \deqn{\mathrm{SNN}(i,j) = \frac{|N_r(i) \cap N_r(j)|}{r},} i.e., the size of
the intersection of their neighbour sets divided by \eqn{r}. The result is
symmetric with ones on the diagonal.
}

\value{
An \eqn{n \times n} symmetric numeric matrix \code{SNN.S} with entries in \eqn{[0,1]}.
}

\references{
  Ghashti, J. S., Hare, W., and J. R. J. Thompson (2025). Variable-weighted adjacency constructions for fuzzy spectral clustering. Submitted.

  Jarvis, R. A., and A. E. Patrick (1973). Clustering using a similarity measure based on shared near neighbors. \emph{IEEE Transactions on Computers, 22}(11), 1025-1034.
}

\seealso{
\code{\link{make.adjacency}}, \code{\link{gen.fuzzy}}, \code{\link{plot.fuzzy}}, \code{\link{rNN.dist}}, \code{\link{find.radius}},
\code{\link{compute.sigma}}, \code{\link{compute.SNN}}, \code{\link{fuzzy.spectral.clustering}}
}

\examples{
set.seed(1)
X <- matrix(rnorm(50), nrow = 10)
D <- as.matrix(dist(X))
r <- 3
S <- exp(-D^2)
SNN <- compute.SNN(S, r)
head(SNN, 5)

# inspect average SNN similarity to nearest neighbour
rowMeans(SNN - diag(diag(SNN)))
}

\keyword{nearest neighbours}
\keyword{shared neighbours}

