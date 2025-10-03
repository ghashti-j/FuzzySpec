\name{fari}
\alias{fari}

\title{Frobenius Adjusted Rand Indice for Comparing Two Partition Matrices}

\description{
Computes fuzzy generalizations of the  Adjusted Rand Index based on Frobenius inner products of membership matrices. These measures extends the Adjusted Rand Index to compare fuzzy partitions.
}

\usage{
fari(a, b)
}

\arguments{
  \item{a}{An \eqn{n \times G_1} matrix of hard or fuzzy cluster memberships, where each row sums to 1.}
  \item{b}{An \eqn{n \times G_2} matrix of hard or fuzzy cluster memberships, where each row sums to 1.}
}

\value{
A single numeric value
  \item{fari}{The Frobenius Adjusted Rand index between \code{a} and \code{b}.}
}

\references{
Andrews, J.L., Browne, R. and C.D. Hvingelby (2022). On Assessments of Agreement Between Fuzzy Partitions. \emph{Journal of Classification, 39}, 326–342.

J.L. Andrews, FARI (2013). GitHub repository, \href{https://github.com/its-likeli-jeff/FARI}{https://github.com/its-likeli-jeff/FARI}
}


\examples{
set.seed(1)
a <- matrix(runif(600), nrow = 200, ncol = 3)
a <- a / rowSums(a)
b <- matrix(runif(600), nrow = 200, ncol = 3)
b <- b / rowSums(b)

fari(a, b)
}

\keyword{cluster}
\keyword{fuzzy}
\keyword{ari}
\keyword{rand index}
