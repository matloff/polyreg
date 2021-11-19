\name{polyFit}
\alias{polyFit}
\alias{predict.polyFit}

\title{Polynomial Fit}

\description{
Fit polynomial regression using a linear or logistic model; predict
new data.}

\usage{
polyFit(xy, deg, maxInteractDeg=deg, use = "lm", pcaMethod=NULL,
     pcaLocation='front', pcaPortion=0.9, glmMethod="one", 
     return_xy=FALSE, returnPoly=FALSE, noisy=TRUE)
\method{predict}{polyFit}(object, newdata, ...)
}

\arguments{
  \item{xy}{Data frame with response variable in the last
     column.  Latter is numeric, except in the classification case: 
     Categorical variables (> 2 levels) must be passed as factors 
     or character variables if \code{use} is 'glm'; an integer vector
     must be used for the for the 'mvrlm' case, or for the 2-class case
     (0s and 1s).}
  \item{deg}{The max degree for polynomial terms.  A term such as uv,
  for instance, is considered degree 2.}
  \item{maxInteractDeg}{The max degree of interaction terms.}
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
  \item{glmMethod}{Defaults to 'one,' meaning the One Versus All Method.
  Use 'all' for All Versus All.}
  \item{newdata}{Data frame, one row for each "X" to be predicted.  Must
     have the same column names as in \code{xy} (without "Y").}
  \item{object}{An item of class 'polyFit' containing output. Can be 
     used with predict().}
  \item{return_xy}{Return data? Default: FALSE}
  \item{returnPoly}{return polyMatrix object? Defaults to FALSE since may be quite large.}
  \item{noisy}{Logical: display messages?}
  \item{...}{Additional arguments for getPoly().}
}

\details{

   The \code{polyFit} function calls \code{getPoly} to generate
   polynomial terms from predictor variables, then fits the generated
   data to a linear or logistic regression model.  (Powers of dummy
   variables will not be generated, other than degree 1, but interaction
   terms will calculated.)

   If \code{pcaMethod} is not \code{NULL}, a principal component
   analysis is performed before or after generating the polynomials.

   When logistic regression for classification is indicated, with more
   than two classes, All-vs-All or One-vs-All methods, coded
   \code{'all'} and \code{'one'}, can be applied to deal with multiclass
   problem.  

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
value(s) for \code{newdata}. It also contains probability for each class as
an attribute named \code{prob}. In the classification case, these will be
the predicted class labels, 1,2,3,...
}

\examples{

N <- 125
xyTrain <- data.frame(x1 = rnorm(N), 
                      x2 = rnorm(N),
                      group = sample(letters[1:5], N, replace=TRUE),
                      score = sample(100, N, replace = TRUE) # final column is y
                      )

pfOut <- polyFit(xyTrain, 2)

# 4 new test points
xTest <- data.frame(x1 = rnorm(4), 
                    x2 = rnorm(4),
                    group = sample(letters[1:5], 4, replace=TRUE))
  
predict(pfOut, xTest) # returns vector of 4 predictions


# spot checks
stopifnot(length(predict(pfOut, xTest)) == nrow(xTest))

}
