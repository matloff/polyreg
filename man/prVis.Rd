\name{prVis}
\alias{prVis}
\alias{addRowNums}

\title{Polynomial-Based Manifold Exploration }

\description{
Polynomial-based alternative to t-SNE, UMAP etc.}

\usage{
prVis(xy, labels = FALSE, deg = 2, scale = FALSE, nSubSam = 0, 
    nIntervals = NULL, saveOutputs = FALSE, cex = 0.5)
addRowNums(np, savedPrVisOut) 
}

\arguments{
  \item{xy}{Data frame with labels, if any, in the last column.}
  \item{labels}{If TRUE, have class labels.}
  \item{deg}{Degree of polynomial.}
  \item{scale}{If TRUE, call \code{scale} on nonlabels data before 
     generating polynomial terms.}
  \item{nSubSam}{Number of random rows of \code{xy} to sample; 0 means
     use the full dataset.} 
  \item{nIntervals}{If labels column is continuous, discretize into this
     many levels.}
  \item{saveOutputs}{Save outputs for use in \code{addRowNums}.}
  \item{cex}{Point size for plot.}
  \item{np}{Number of points to label in plot.}
  \item{savedPrVisOut}{Output save from a previous call, so can avoid
     duplicate computation.}
}

\details{

A number of "nonlinear" analogs of Principle Components Analysis (PCA)
have emerged, such as ICA, t-SNE, UMAP and so on. Intuitively, an
approach based on polynomials may be effective too.  Specifically,
\code{prVis} first expands \code{xy} to polynomial terms, then applies
PCA to the result.

}

\value{
}

\examples{


}
