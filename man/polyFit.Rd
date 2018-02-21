\name{polyFit}
\alias{polyFit}

\title{Polynomial Fit}

\description{
Fit a polynomial regression using linear regreesion or logistic regreesion for classification.
}

\usage{
polyFit(xy, maxDeg, maxInteractDeg = maxDeg, use = "lm", trnProp = 0.8,
        pcaMethod = FALSE, pcaPortion = 0.9, glmMethod = "all")
}

\arguments{
  \item{xy}{Data matrix or dataframe with response variable at the last column.}
  \item{maxDeg}{The max degree to generate polynomial terms.}
  \item{maxInteractDeg}{The max degree of dummny and nondummy interaction predictor variables.}
  \item{use}{Can be "lm" for using linear regreesion, and "glm" for using logistic regression.}
  \item{trnProp}{The portion of data to be used as training set.}
  \item{pcaMethod}{Whether to use PCA (can be TRUE or FALSE).}
  \item{pcaPortion}{The portion of principal components to be used.}
  \item{glmMethod}{When the response variable has more than 2 classes, it can be "all" for All-vs-All method, and "one" for One-vs-All method.}
}

\details{

   The \code{polyFit} function takes in a dataframe and use \code{getPoly}
   function to generate polynomial terms of predictor variables, then it fits
   the generated data to linear or logistic regreesion models. If
   \code{pcaMethod == TRUE}, the principal component analysis is performed
   using the indicated portion of principal components before fitting the
   models. When logistic regression for classification is indicated,
   All-vs-All or One-vs-All methods can be applied to deal with
   multiclass problem.
}
\value{
The return value of \code{polyFit} is an R list of mean absolute error (for \code{lm}) or accuracy (for \code{glm}). The i-th element of the list is for degree = i.
}
\examples{
y <- mtcars[,1]
data <- cbind(mtcars[,-1], y) # make y column the last column
pf1 <- polyFit(data,2,1,"lm",0.8,FALSE)
pf2 <- polyFit(data,5,1,"lm",0.8,TRUE,0.9)

mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
pf3 <- polyFit(mydata,2,1,"glm",0.8,FALSE)
pf4 <- polyFit(mydata,2,1,"glm",0.8,FALSE,0.9,"one")
}
