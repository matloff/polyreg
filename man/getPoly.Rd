\name{getPoly}
\alias{getPoly}
\alias{polyMatrix}

\title{Get polynomial terms}

\description{
Generate polynomial terms of predictor variables for a
data frame or data matrix.
}

\usage{
getPoly(xdata = NULL, deg = 1, maxInteractDeg = deg,
        Xy = NULL, standardize = FALSE,
        noisy = TRUE, intercept = FALSE, returnDF = TRUE, 
        modelFormula = NULL, retainedNames = NULL, ...)
}



\arguments{
  \item{xdata}{Data matrix or data frame without response variable. Categorical variables (> 2 levels) should be passed as factors, not dummy variables or integers, to ensure the polynomial matrix is constructed properly.}
  \item{deg}{The max degree of power terms. Default 1 so just returns model matrix by default.}
  \item{maxInteractDeg}{The max degree of nondummy interaction terms. x1 * x2 is degree 2. x1^3 * x2^2 is degree 5. Implicitly constrained by deg. For example, if deg = 3 and maxInteractDegree = 2, x1^1 * x2^2 (i.e., degree 3) will be included but x1^2 * x2^2 (i.e., degree 4) will not.}
  \item{Xy}{ the dataframe with the response in the final column (provide xdata or Xy but not both).Categorical variables (> 2 levels) should be passed as factors, not dummy variables or integers, to ensure the polynomial matrix is constructed properly.}
  \item{standardize}{ standardize all continuous variables? (Default: FALSE.)}
  \item{noisy}{Output progress updates? (Default: TRUE.)}
  \item{intercept}{Include intercept? (Default: FALSE.)}
  \item{returnDF}{Return a data.frame (as opposed to model.matrix)? (Default: TRUE.)}
  \item{modelFormula}{ Internal use. Formula used to generate the training model matrix. Note: anticipates that polynomial terms are generated using internal functions of library(polyreg). Also, providing modelFormula bypasses deg and maxInteractDeg.}
  \item{retainedNames}{ Internal use. colnames of polyMatrix object$xdata. Requires modelFormula be inputted as well.}
  \item{...}{ additional arguments to be passed to model.matrix() via polyreg:::model_matrix(). Note na.action = "na.omit".}
}

\details{

   The \code{getPoly} function takes in a data frame or data matrix and
   generates polynomial terms of predictor variables.

   Note the subtleties involving dummy variables.  The square, cubic and
   so on terms are the same as the original variable, and the various
   duplicates must be eliminated.

   Similarly, after dummy variable are created from a categorical
   variable having more than two levels, the resulting columns will be
   orthogonal to each other.  In almost
   all cases, this argument should be set to TRUE at the training stage, and
   then in predictions one should use the vector of names in the
   component in the return value;
   \code{predict.polyFit} does the latter automatically.

   % Also, this function treats predictor variables that have only two kinds
   % of values categorical variables, so the data needs preprocessing to
   % multiple dummy variables if needed.
}

\value{
The return value of \code{getPoly} is a \code{polyMatrix}
object.  This is an S3 class containing a model.matrix \code{
xdata} of the generated polynomial terms. The predictor
variables have column names V1, V2, etc. The object also contains
\code{modelFormula}, the formula used to construct the model matrix, and 
\code{XtestFormula}, the formula which should be used out-of-sample 
(when y_test is not available).
}

\examples{
N <- 125
rawdata <- data.frame(x1 = rnorm(N), 
                      x2 = rnorm(N),
                      group = sample(letters[1:5], N, replace=TRUE),
                      z = sample(c("treatment", "control"), N, replace=TRUE),
                      result = sample(c("win", "lose", "tie"), N, replace=TRUE))
head(rawdata)

P <- length(levels(rawdata$group)) - 1 + 
     length(levels(rawdata$z)) - 1 + 
     length(levels(rawdata$result)) - 1 + 
     sum(unlist(lapply(rawdata, is.numeric)))

# quadratic polynomial, includes interactions 
# since maxInteractDeg defaults to deg
X <- getPoly(rawdata, 2)$xdata 
ncol(X) # 40

# cubic polynomial, no interactions
X <- getPoly(rawdata, 3, 1)$xdata
ncol(X) # 13

# cubic polynomial, interactions
X <- getPoly(rawdata, 3, 2)$xdata
ncol(X) # 58

# cubic polynomial, interactions
X <- getPoly(rawdata, 3)$xdata
ncol(X) # 101

# making final column the response variable, y
# results in TRUE (fewer columns)
ncol(getPoly(Xy=rawdata, deg=2)$xdata) < ncol(getPoly(rawdata, 2)$xdata)

# preparing polynomial matrices for crossvalidation
# getPoly() returns a polyMatrix() object containing XtestFormula
# which should be used to ensure factors are handled correctly out-of-sample
Xtrain <- getPoly(rawdata[1:100,],2)
Xtest <- getPoly(rawdata[101:125,], 2, modelFormula = Xtrain$XtestFormula)
    
