\name{clustering.accuracy}
\alias{clustering.accuracy}

\title{Clustering Accuracy with Optimal Label Matching}

\description{
Computes the fraction of correctly classified observations between two label vectors after optimally matching cluster labels using \code{Thresher::matchLabels}.
}

\usage{
clustering.accuracy(A, B)
}

\arguments{
  \item{A}{An integer or character vector of cluster labels of length \eqn{n}.}
  \item{B}{An integer or character vector of cluster labels of length \eqn{n}.}
}

\details{
The function creates the contingency table \code{table(A, B)}, permutes columns to best align labels using \code{\link[Thresher]{matchLabels}},
and returns the sum of the diagonal divided by \eqn{n}.

Inputs must have equal length and the same number of unique labels; otherwise an error is given.
}

\value{
A single numeric value in \eqn{[0, 1]}: the accuracy after optimal label matching.
}

\references{
  K. R. Coombes (2025). \emph{Thresher: Threshing and Reaping for Principal Components}. R package version 1.1.5.
}

\seealso{
\code{\link{fuzzy.spectral.clustering}}, \code{\link{gen.fuzzy}}, \code{\link{plot.fuzzy}}, \code{\link[Thresher]{matchLabels}}
}

\examples{
set.seed(1)
n <- 200
k <- 3
A <- sample.int(k, n, replace = TRUE) # assumed true clustering labels
perm <- sample.int(k) # assumed predicted labels (sampled by permutating)
B <- perm[A]
flips <- sample.int(n, 20) # add some error a few errors
B[flips] <- sample.int(k, length(flips), replace = TRUE)

clustering.accuracy(A, B)
}

\keyword{clustering}
\keyword{classification}
\keyword{accuracy}
\keyword{confusion matrix}
