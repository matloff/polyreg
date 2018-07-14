\name{polyFit}
\alias{polyFit}
\alias{predict.polyFit}

\title{Polynomial Fit}

\description{
Fit polynomial regression using a linear or logistic model; predict on 
new data.}

\usage{
polyFit(xy,deg,maxInteractDeg,use='lm',pcaMethod=NULL,pcaPortion=0.9,
   glmMethod='one',printTimes=TRUE,polyMat=NULL,cls=NULL,dropout=0) 
predict.polyFit(object, newdata)
}

\arguments{
  \item{xy}{Data frame with response variable in the last 
     column.  In the classification case, response is class ID.}
  \item{maxDeg}{The max degree for polynomial terms.}
  \item{maxInteractDeg}{The max degree of interaction  terms}
  \item{use}{Set to lm for linear regression or 'glm' for 
     using logistic regression.}
  \item{pcaMethod}{NULL for no PCA. For PCA, can be either 'prcomp' 
     (use the \code{prcomp} function) or 'RSpectra' (use
     the \code{eigs} function in the \pkg{RSpectra} package, on a sparse
     matrix from the \pkg{Matrix} package).}
  \item{pcaPortion}{The portion of total variance of principal components 
     to be used.}
  \item{polyMat}{If TRUE, \code{getPoly} has already been called and
     thus need not be called now.}
  \item{cls}{Virtual cluster, for parallel computation (One Vs. All
     case).}
  \item{dropout}{If nonzero, delete randomly chosen columns from the
     polynomial terms matrix.  Number of deleted columns determined by
     the value of this argument.}
  \item{newdata}{Data frame, one row for each "X" to be predicted.  Must
     have the same column names as in \code{xy} (without "Y").}
}

\details{

   The \code{polyFit} function calls \code{getPoly}
   to generate polynomial terms from predictor variables, then fits
   the generated data to a linear or logistic regression model. 
   (Powers of dummy variables will not be generated, other than degree
   1, but interaction terms will calculated.
   
   If \code{pcaMethod} is not \code{NULL}, a principal component analysis 
   is performed before fitting the models. The argument \code{pcaPortion}
   determines how many principal components to be used, according to
   proportion of total variance.  When logistic regression for 
   classification is indicated, All-vs-All or One-vs-All methods can 
   be applied to deal with multiclass problem.
}
\value{
The return value of \code{polyFit()} is an \code{polyFit} object.  The
orginal arguments are retained, 
}
\examples{
y <- mtcars[,1]
data <- cbind(mtcars[,-1], y) # make y column the last column
f <- polyFit(data,2)
pred <- predict(f,data[,-ncol(data)])
# e.g. Honda Civic prediction is 30.26006

}
