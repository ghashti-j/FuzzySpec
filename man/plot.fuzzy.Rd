\name{plot.fuzzy}
\alias{plot.fuzzy}

\title{Plot 2D Fuzzy Data with Optional Uncertainty Sizing and Cluster Colouring}

\description{
Creates a ggplot of 2D points with optional colouring by hard labels and optional observation-size mapping to fuzzy uncertainty.
}

\usage{
\method{plot}{fuzzy}(x, plotFuzzy = TRUE, colorCluster = TRUE, ...)
}

\arguments{
  \item{x}{A list as returned by \code{\link{gen.fuzzy}}, containing \code{X}, \code{U}, \code{y}, and \code{k}.}
  \item{plotFuzzy}{Logical; if \code{TRUE}, map observation size to uncertainty \eqn{1-\max_j U_{ij}}.}
  \item{colorCluster}{Logical; if \code{TRUE}, colour points by the hard cluster label \code{y}.}
  \item{...}{Additional arguments (currently unused).}
}

\details{
The plotting aesthetics can be modified as follows:
\itemize{
  \item If \code{plotFuzzy} and \code{colorCluster} are both \code{TRUE} (default), the plot contains cluster coloured observations that are size scaled by uncertainty.
  \item If only \code{plotFuzzy} is \code{TRUE}, the plot contains monochrome coloured observations that are size scaled by uncertainty.
  \item If only \code{colorCluster} is \code{TRUE}, the plot contains cluster coloured observations with fixed size.
  \item If \code{plotFuzzy} and \code{colorCluster} are both \code{FALSE}, the plot is monochrome coloured observations with fixed size.
}
}

\value{
A \code{ggplot} object.
}

\references{
  H. Wickham (2016). \emph{ggplot2: Elegant Graphics for Data Analysis}. Springer-Verlag New York.
}

\seealso{
\code{\link{gen.fuzzy}}, \code{\link[ggplot2]{ggplot}}
}

\examples{
set.seed(1)
d1 <- gen.fuzzy(n = 600, dataset = "gaussian", seed = 1)
p1 <- plot.fuzzy(d1)
p1 # default

p2 <- plot.fuzzy(d1, plotFuzzy = TRUE, colorCluster = FALSE)
p2 # only uncertainty sizing, monochrome

p3 <- plot.fuzzy(d1, plotFuzzy = FALSE, colorCluster = TRUE)
p3 # only coloured by cluster, no uncertainty sizing
}

\keyword{visualization}
\keyword{clustering}
\keyword{fuzzy}
