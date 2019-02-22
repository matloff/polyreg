\name{xvalPoly,xvalKf,xvalDnet}
\alias{xvalPoly}
\alias{xvalKf}
\alias{xvalDnet}

\title{
Cross validation on Polynomial regression
}
\description{
Separate the dataset to training and test set, and obtain the mean absolute error (for linear regression) or accuracy (for logistic regression).
}
\usage{
xvalPoly(xy, maxDeg, maxInteractDeg = maxDeg, use = "lm",
               pcaMethod = NULL, pcaPortion = 0.9,
               glmMethod = "one",
               nHoldout = min(10000, round(0.2*nrow(xy))),
               pcaLocation='front',
               yCol = NULL, cls=NULL, startDeg=1)
}

\arguments{
  \item{xy}{Data matrix or dataframe with response variable in the last column.}
  \item{maxDeg}{Max degree for power terms.}
  \item{maxInteractDeg}{Max degree of interaction terms.}
  \item{startDeg}{Starting degree. Enables user to skip lower degree models. (But defaults to 1.)}
  \item{use}{Can be "lm" for using linear regreesion, and "glm" for
     using logistic regression.}
  \item{pcaMethod}{Should be either NULL (no PCA), "prcomp" (base R implementation of PCA), or "RSpectra" (recent implementation of PCA).}
  \item{pcaPortion}{If pcaMethod is provided, use components up to this
     proportion of total variance.}
  \item{pcaLocation}{Is PCA performed before or after the (polynomial) model matrix is constructed?}
  \item{glmMethod}{For classification problems.  If there are more than
     two classes, this can be "all" for All-vs-All method,
     or "one" for One-vs-All method.}
  \item{nHoldout}{Size of test test.}
  \item{yCol}{If provided, the column where y is found (i.e., other than the last colum, where it is expected).}
  \item{cls}{cluster (for clusterEvalQ()).}

}

\details{

   The \code{xvalPoly} function divides the data to training and test sets,
   and use \code{polyFit} to generate models using training data and use
   \code{polyFit.predict} to generate results on test set, and compare the
   results. The \code{xvalKf} does the same for \code{kerasformula}
   neural networks package.
}
\value{
The return value of \code{xvalPoly} is an R vector of mean absolute
error (for \code{lm}) or probability of correct classification
(for \code{glm}).  The i-th element of the vector is for degree i.
}

\examples{
y <- mtcars[,1] # make mpg the outcome
Xy <- as.data.frame(cbind(mtcars[,-1], y)) # make y column the last column
xvalPoly(Xy,2,use="lm") # not using pca
xvalPoly(Xy,2, use="lm", pcaMethod="prcomp", 0.8) # using pca on X first
}

