\name{make.adjacency}
\alias{make.adjacency}

\title{A General Framework for Adjacency Matrix Construction}

\description{
Builds an \eqn{n \times n} adjacency matrix from data using Euclidean or variable-weighted distances, with optional locally-adaptive scalings and optional variable weighting by shared-nearest-neighbours (SNN), similarity ranks (SIM), or both. The options reproduce a family of adjacencies including the variants described by Ghashti, Hare and Thompson (2025).
}

\usage{
make.adjacency(data,
                method = "vw",
                isLocWeighted = FALSE,
                isModWeighted = FALSE,
                isSparse = FALSE,
                ModMethod = NULL,
                scale = FALSE,
                sig = 1,
                radius = NULL,
                cv.method = "cv.ls")
}

\arguments{
  \item{data}{Numeric matrix or data frame of size \eqn{n \times p}. Required.}
  \item{method}{Distance construction  \code{"eu"} for Euclidean; \code{"vw"} for variable-weighted scaling (see Details). Defaults to \code{"v"}.}
  \item{isLocWeighted}{Logical. If \code{TRUE}, use locally-adaptive scalings \eqn{\sigma_i} (Zelnik–Perona) to self-tune the kernel; otherwise use a global scale \code{sig}. Defaults to \code{FALSE}.}
  \item{isModWeighted}{Logical. If \code{TRUE}, apply a weighting matrix \eqn{M} based on SNN and/or SIM (see Details). Defaults to \code{FALSE}}
  \item{isSparse}{Logical. If \code{TRUE} and \code{isModWeighted = TRUE}, then SNN (\code{ModMethod = "snn"}) and/or SIM (\code{ModMethod = "sim"}) arguments sparisfy the adjacency matrix (see Details). Defaults to \code{FALSE}.}
  \item{ModMethod}{One of \code{"snn"}, \code{"sim"}, or \code{"both"} when \code{isModWeighted=TRUE}.}
  \item{scale}{Logical; standardize columns of \code{data} before distance construction.}
  \item{sig}{Positive numeric value for global kernel width, used only when \code{isLocWeighted=FALSE}.}
  \item{radius}{Integer \eqn{r}. If \code{NULL}, \eqn{r} is estimated with \code{\link{find.radius}} on the constructed distance matrix.}
  \item{cv.method}{Bandwidth selector for \code{method="vw"} passed to \code{np::npudensbw}; one of \code{"cv.ml"} or \code{"cv.ls"}.}
}

\details{
\strong{Step 1: Distance.}
\itemize{
  \item \code{method="eu"}: \eqn{D_{ij}=\|x_i-x_j\|_2}.
  \item \code{method="vw"}: compute product-kernel bandwidths \code{h} via \code{np::npudensbw}, set feature-weights \eqn{w_j=1/h_j^2}, rescale data as \eqn{\tilde{x}_{ij}=\sqrt{w_j}\,x_{ij}}, then \eqn{D_{ij}=\|\tilde{x}_i-\tilde{x}_j\|_2}. (Variable-weighted metric.)
}

\strong{Step 2:Similarity kernel.}
\itemize{
  \item \emph{Locally-adapative scaling} from Zelnik-Manor: if \code{isLocWeighted=TRUE}, compute \eqn{\sigma_i} as the distance to the \eqn{r}-th neighbour with \code{\link{compute.sigma}} and set
  \deqn{S_{ij}=\exp\!\Big(-D_{ij}^2 / (\sigma_i \sigma_j)\Big),\quad S_{ii}=1.}
  \item \emph{Global scale}: if \code{isLocWeighted=FALSE},
  \deqn{S_{ij}=\exp\!\Big(-D_{ij}^2 / \sigma^2\Big),\quad S_{ii}=1.}
}

\strong{Step 3: Weighting Matrix (optional).}

Let \eqn{\text{SNN}_{ij}} be the shared-\eqn{r}-NN overlap fraction (see \code{\link{compute.SNN}}) for observations \eqn{i} and \eqn{j}, and
\eqn{\rho_i} be the \eqn{(r+1)}-largest entry of observation \eqn{i} in matrix \eqn{S}. Define \eqn{\text{SIM}_{ij}=\sqrt{\rho_i \rho_j}}.
For \code{isModWeighted=TRUE} we have the following options

\itemize{
  \item \code{ModMethod="snn"}: \eqn{M_{ij} = \begin{cases}
        0.5\,(1+\text{SNN}_{ij}), & \text{if } \texttt{isSparse=FALSE},\\
        \text{SNN}_{ij}, & \text{if } \texttt{isSparse=TRUE}.
      \end{cases}}
  \item \code{ModMethod="sim"}: \eqn{M_{ij} = \begin{cases}
        0.5\,(1+\text{SIM}_{ij}), & \text{if } \texttt{isSparse=FALSE},\\
        \text{SIM}_{ij}, & \text{if } \texttt{isSparse=TRUE}.
      \end{cases}}
  \item \code{ModMethod="both"}: \eqn{M_{ij} = \begin{cases}
        0.25\,(1+\text{SIM}_{ij})\,(1+\text{SNN}_{ij}), & \text{if } \texttt{isSparse=FALSE},\\
        \text{SIM}_{ij}\cdot \text{SNN}_{ij}, & \text{if } \texttt{isSparse=TRUE}.
      \end{cases}}
}
The returned adjacency is \eqn{W = S \circ M} when \code{isModWeighted = TRUE}, otherwise \eqn{W=S}.
These choices align with the table of named adjacencies (see Mapping below).
}

\section{Mapping to named adjacencies}{

Relating to the paper by Ghashti, Hare, and Thompson (2025), let "vw" denote the variable-weighted distance \code{method="vw"} and \code{"eu"} for the traditional squared Euclidean distance; "la" denotes locally-adaptive scaling when \code{isLocWeighted=TRUE} (Zelnik-Manor and Perona, 2004);
"id" denotes identity for \code{isModWeighted = FALSE}, and "sim", "snn" and "simsnn" denote weightings for \eqn{M} described above.

To reproduce results from the 2025 paper, below are a few examples of adjacency construction:
\itemize{
  \item \code{vw-id}: \eqn{\exp(-D^2/\sigma^2)} with \code{method="vw"}, \code{isLocWeighted=FALSE}, \code{isModWeighted=FALSE}.
  \item \code{vwla-id}: \eqn{\exp(-D^2/(\sigma_i\sigma_j))} with \code{method="vw"}, \code{isLocWeighted=TRUE}, \code{isModWeighted=FALSE}.
  \item \code{vw-sim}: \eqn{0.5\,\exp(-D^2/\sigma^2)\,(1+\text{SIM})} with \code{method="vw"}, \code{isLocWeighted=FALSE}, \code{isModWeighted=TRUE}, \code{ModMethod="sim"}, \code{isSparse=FALSE}.
  \item \code{vw-snns}: \eqn{0.5\,\exp(-D^2/\sigma^2)\,(1+\text{SNN})} with \code{method="vw"}, \code{isLocWeighted=FALSE}, \code{isModWeighted=TRUE}, \code{ModMethod="snn"}, \code{isSparse=TRUE}.
  \item \code{vw-simsnns}: \eqn{0.25\,\exp(-D^2/\sigma^2)\,(1+\text{SIM})(1+\text{SNN})} with \code{method="vw"}, \code{isLocWeighted=FALSE}, \code{isModWeighted=TRUE}, \code{ModMethod="both"}, \code{isSparse=TRUE}.
  }

Also note that:
\itemize{
  \item If \code{radius} is \code{NULL}, \eqn{r} is chosen adaptively via \code{\link{find.radius}} on the constructed distance matrix.
  \item \code{method="vw"} requires \code{\link[np]{npudensbw}} for variable weighted bandwidths, with default \code{np::npudensbw(data, bwmethod = cv.method, nmulti = 3)} (Hayfield and Racine, 2008).
  }
}

\section{Notes}{
\itemize{
  \item When \code{r} is determined by \code{\link{find.radius}}, we implement a modified version of the Natural Neighbors algorithm from Zhu et al. (2016).
  \item SNN is a modified version of the Shared Nearest Neighbors algorithm from Jarvis and Patrick (1973).
  \item More informating on locally-adaptive scalings are seen in Zelnik-Manor and Perona (2004).
}
}


\value{
An \eqn{n \times n} numeric adjacency matrix \eqn{W} with ones on the diagonal.
}


\references{
  Ghashti, J. S., Hare, W., and J. R. J. Thompson (2025). Variable-weighted adjacency constructions for fuzzy spectral clustering. Submitted.

  Hayfield, T., and J. S. Racine (2008). Nonparametric Econometrics: The np Package. \emph{Journal of Statistical Software 27}(5).

  Jarvis, R. A., and A. E. Patrick (1973). Clustering using a similarity measure based on shared near neighbors. \emph{IEEE Transactions on Computers, 22}(11), 1025-1034.

  Zelnik-Manor, L., and P. Perona (2004). Self-tuning spectral clustering. \emph{Advances in Neural Information Processing Systems, 17}.

  Zhu, Q., Feng, J., and J. Huang (2016). Natural neighbor: A self-adaptive neighborhood method without parameter K. \emph{Pattern Recognition Letters, 80}, 30-36.
}

\seealso{
\code{\link{gen.fuzzy}}, \code{\link{plot.fuzzy}}, \code{\link{rNN.dist}}, \code{\link{find.radius}},
\code{\link{compute.sigma}}, \code{\link{compute.SNN}}, \code{\link{fuzzy.spectral.clustering}}, \code{\link[np]{npudensbw}}
}

\examples{
set.seed(1)
X <- scale(matrix(rnorm(200), 100, 2))


W1 <- make.adjacency(X,
                    method = "eu",
                    isLocWeighted = TRUE) # "eula-id" named adjacency

W2 <- make.adjacency(X,
                    method = "vw",
                    isLocWeighted = TRUE) # "vwla-id" named adjacency

# compare W(xi,xj) i,j = 1,...,5 for eu/vw pair W1 and W2
W1[1:5,1:5]
W2[1:5,1:5]

W3 <- make.adjacency(X,
                    method = "eu",
                    isLocWeighted = TRUE,
                    isModWeighted = TRUE,
                    ModMethod = "snn",
                    isSparse = FALSE) # "eula-snn" named adjacency

W4 <- make.adjacency(X,
                    method = "vw",
                    isLocWeighted = TRUE,
                    isModWeighted = TRUE,
                    ModMethod = "snn",
                    isSparse = FALSE) # "vwla-snn" named adjacency

# compare W(xi,xj) i,j = 1,...,5 for eu/vw pair W3 and W4
W3[1:5,1:5]
W4[1:5,1:5]
}

\keyword{similarity}
\keyword{adjacency}
\keyword{variable weighting}
\keyword{spectral clustering}
\keyword{nearest neighbors}
