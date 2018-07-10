\name{polyFit}
\alias{polyFit}
\alias{predict.polyFit}

\title{Polynomial Fit}

\description{
Fit a polynomial regression using linear regreesion or logistic regreesion for classification, and predict on a new dataset.
}

\usage{
polyFit(xy, deg, maxInteractDeg, use = "lm", pcaMethod=FALSE, pcaPortion = 0.9)
predict.polyFit(object, newdata)
}

\arguments{
  \item{xy}{Data matrix or dataframe with response variable at the last column.}
  \item{maxDeg}{The max degree to generate polynomial terms.}
  \item{maxInteractDeg}{The max degree of dummny and nondummy interaction predictor variables.}
  \item{use}{Can be "lm" for using linear regreesion, and "glm" for using logistic regression.}
  \item{pcaMethod}{default is NULL, can be either "prcomp" (use the prcomp() function to compute PCA) or "RSpectra" (use sparse Matrix and compute eigenvalues/vectors to compute PCA).}
  \item{pcaPortion}{The portion of principal components to be used.}
}

\details{

   The \code{polyFit} function takes in a dataframe and use \code{getPoly}
   function to generate polynomial terms of predictor variables, then it fits
   the generated data to linear or logistic regreesion models. If
   \code{pcaMethod} is not \code{NULL}, the principal component analysis is performed
   using the indicated portion of principal components before fitting the
   models. When logistic regression for classification is indicated,
   All-vs-All or One-vs-All methods can be applied to deal with
   multiclass problem.
}
\value{
The return value of \code{polyFit()} is an \code{polyFit} object.
}
\examples{
y <- mtcars[,1]
data <- cbind(mtcars[,-1], y) # make y column the last column
f <- polyFit(data,2)
pred <- predict(f,data[,-ncol(data)])
# e.g. Honda Civic prediction is 30.26006

}
