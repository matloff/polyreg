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

Once a plot is displayed, \code{addRowNums} can be used to add
row-number IDs of random points, to gain further insight into the data.

}

\value{
If \code{saveOutputs} is set, an R list is returned, with comppnents
\code{gpOut}, the generated polynomial matrix, and \code{prout}, the
return value from the call to \code{prcomp}.
}

\examples{

getPE()  # prgeng data, included in pkg
# may want to predict wage; look at some predictors
pe1 <- pe[,c(1,2,6:16)]
z <- prVis(pe1,nSubSam=2000,saveOutputs=T,labels=FALSE)  
# get a bunch of streaks; why?
# call addRowNums() (not shown); discover that points on the same streak
# tend to have same combination of sex, education and occupation; moving
# along a streak mainly consists of variying age

print('see data/SwissRoll for another example')

}

