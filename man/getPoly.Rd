\name{getPoly}
\alias{getPoly}
\alias{polyMatrix}

\title{Get polynomial terms}

\description{
Generate polynomial terms of predictor variables for a
data frame or data matrix.
}

\usage{
getPoly(xydata, deg, maxInteractDeg = deg)
}

\arguments{
  \item{xydata}{Data matrix or data frame without response variable.}
  \item{deg}{The max degree of power terms.}
  \item{maxInteractDeg}{The max degree of nondummy interaction terms. x1 * x2 is degree 2. x1^3 * x2^2 is degree 5. Implicitly constrained by deg. For example, if deg = 3 and maxInteractDegree = 2, x1^1 * x2^2 (i.e., degree 3) will be included but x1^2 * x2^2 (i.e., degree 4) will not.}
  \item{Xy}{ the dataframe with the response in the final column (provide xdata or Xy but not both). Factors with more than two levels should not be inputted as integers.}
  \item{modelFormula}{ Internal use. Formula used to generate the training model matrix. Note: anticipates that polynomial terms are generated using internal functions of library(polyreg). Also, providing modelFormula bypasses maxDeg and maxInteractDeg.}
  \item{standardize}{ standardize all continuous variables? (Default: FALSE.)}
  \item{intercept}{Include intercept? Default: FALSE.}
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
object.  This is an S3 class containing a data frame \code{
xdata} of the generated polynomial terms. The predictor
variables have column names V1, V2, etc. The object also contains
a vector \code{retainedNames}, the names of the non-all-0 columns.
}

\examples{

x1 <- 1:4
z1 <- c(0,1,0,0)
z2 <- c(1,0,0,0)
z3 <- c(0,0,1,0)
xz <- cbind(x1,z1,z2,z3)
getPoly(xz,2)
# xdata component of output:
#   V1 V2 V3 V4 V5 V9 V10 V11
# 1  1  0  1  0  1  0   1   0
# 2  2  1  0  0  4  2   0   0
# 3  3  0  0  1  9  0   0   3
# 4  4  0  0  0 16  0   0   0
# V1-V4 are copies of x1,z1,z2,z3, 1st-degree terms; V5 is x1^2; V9 is
# product x1 * z2; etc.; note that, V6-V8 were  not retained

}

