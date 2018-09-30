\name{polyFit}
\alias{polyFit}
\alias{predict.polyFit}

\title{Polynomial Fit}

\description{
Fit polynomial regression using a linear or logistic model; predict 
new data.}

\usage{
polyFit(xy,deg,maxInteractDeg,use='lm',pcaMethod=NULL,pcaLocation='front',
   pcaPortion=0.9,glmMethod='one',cls=NULL,del0cols=TRUE) 
predict.polyFit(object,newdata)
}

\arguments{
  \item{xy}{Data frame with response variable in the last 
     column.  In the classification case, response is class ID, stored
     in a vector, not as a factor.}
  \item{maxDeg}{The max degree for polynomial terms.}
  \item{maxInteractDeg}{The max degree of interaction  terms}
  \item{use}{Set to 'lm' for linear regression, 'glm' for 
     logistic regression, or 'mvrlm' for multivariate-response \code{lm}.}
  \item{pcaMethod}{NULL for no PCA. For PCA, can be either 'prcomp' 
     (use the \code{prcomp} function) or 'RSpectra' (use
     the \code{eigs} function in the \pkg{RSpectra} package.}
  \item{pcaLocation}{In case PCA is applied, specify 'front' to have PCA
     calculated before forming polynomials, otherwise 'back.}
  \item{pcaPortion}{If less than 1.0, use as many principal components
     so as to achieve this portion of total variance.  Otherwise, use
     this many components.  In the 'RSpectra' case, this value must be
     an integer of 1 or more.}
  \item{cls}{Virtual cluster, for parallel computation (One Vs. All
     case).}
  \item{del0cols}{Delete all-0 columns in polynomial data frame.
     See \code{\link[=exportAttribute]{getPoly}} 
     documentation.}
  \item{newdata}{Data frame, one row for each "X" to be predicted.  Must
     have the same column names as in \code{xy} (without "Y").}
}

\details{

   The \code{polyFit} function calls \code{getPoly} to generate
   polynomial terms from predictor variables, then fits the generated
   data to a linear or logistic regression model.  (Powers of dummy
   variables will not be generated, other than degree 1, but interaction
   terms will calculated.
   
   If \code{pcaMethod} is not \code{NULL}, a principal component
   analysis is performed before or after generating the polynomials. 
   
   When logistic regression for classification is indicated, with more
   than two classes, All-vs-All or One-vs-All methods, coded
   \code{'all'} and \code{'one'}, can be applied to deal with multiclass
   problem.  Multinomial logit (\code{'multilog'}) is also available.

   Under the 'mvrlm' option in a classification problem, \code{lm} is
   called with multivariate response, using \code{cbind} and dummy
   variables for class membership as the response.  Since predictors are
   used to form polynomials, this should be a reasonable model, and is
   much faster than 'glm'.

}
\value{
The return value of \code{polyFit()} is an \code{polyFit} object.  The
orginal arguments are retained, along with the fitted models and so on.

The prediction function \code{predict.polyFit} returns the predicted
value(s) for \code{newdata}. In the classification case, these will be
the predicted class labels, 1,2,3,...
}

\examples{

getPE(Dummies=T)  # prgeng data
pe1 <- pe[,c(1,2,4,6,7,12:16,3)]   # select some predictors
pfout <- polyFit(pe1,2) 
# predict worker like pe1[1,] but age 42 with MS degree
newdata <- pe1[1,-11] 
newdata[1,1] <- 42 
newdata[1,4] <- 1 
predict(pfout,newdata[,-11])  # 81022.27

}
