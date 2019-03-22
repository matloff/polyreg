\name{polyFit}
\alias{polyFit}
\alias{predict.polyFit}

\title{Polynomial Fit}

\description{
Fit polynomial regression using a linear or logistic model; predict
new data.}

\usage{
polyFit(xy, deg, maxInteractDeg=deg, use = "lm", pcaMethod=NULL,
     pcaLocation='front', pcaPortion=0.9, glmMethod="one", cls=NULL)
\method{predict}{polyFit}(object, newdata, ...)
}

\arguments{
  \item{xy}{Data frame with response variable in the last
     column.  In the classification case, response is class ID, stored
     in a vector, not as a factor. Categorical variables (> 2 levels) should be passed as factors, not dummy variables or integers, to ensure the polynomial matrix is constructed properly.}
  \item{deg}{The max degree for polynomial terms.}
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
  \item{glmMethod}{Defaults to "one."}
  \item{cls}{Virtual cluster, for parallel computation (One Vs. All
     case).}
  \item{newdata}{Data frame, one row for each "X" to be predicted.  Must
     have the same column names as in \code{xy} (without "Y").}
  \item{object}{An item of class 'polyFit' containing output. Can be used with predict().}
  \item{...}{Additional arguments for getPoly().}
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

data(prgeng) # US Census engineer wage data 
indx <- seq(1, nrow(prgeng), 10)
xy <- prgeng[indx, c(1:4,6,5)] # Y value (wageinc here) must be last
xy <- toFactors(xy, 2:4)  # convert education, occupation, and sex to factors
                          # (don't pass categorical variables as integer codes when fitting)

# regression example
pfout <- polyFit(xy, 2)  # try degree 2 model

# say we want to predict wage for a woman, age 35, some college,
# occupation 101, worked 52 weeks
newx <- data.frame(age=35,educ=11,occ=101,sex=2,wkswrkd=52) 
predict(pfout,newx)  # 38830.54 

# classification example, logistic
xy <- xy[,c(1,2,4:6,3)]  # Y value (occupaton here) must be last
pfout <- polyFit(xy,2,use='glm')  # try degree 2 model
# say we want to predict occupation for a woman, age 35, college grad,
# occupation 101, worked 52 weeks; and for a man with the same
# characteristics
newx <- data.frame(age=c(35,35),
                   educ=c(13,13),
                   sex=c(2,1),
                   wageinc=c(62000,62000),
                   wkswrkd=c(52,52)) 
preds <- predict(pfout,newx)  
levels(xy$occ)[preds]  # translate from class number; 102, 102



}
