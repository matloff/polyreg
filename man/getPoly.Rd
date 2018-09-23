\name{getPoly}
\alias{getPoly}
\alias{polyMatrix}

\title{Get polynomial terms}

\description{
Generate polynomial terms of predictor variables for a dataframe or data matrix.
}

\usage{
getPoly(xydata, deg, maxInteractDeg = deg, del0cols=TRUE)
}

\arguments{
  \item{xydata}{Data matrix or dataframe without response variable.}
  \item{deg}{The max degree to generate polynomial terms.}
  \item{maxInteractDeg}{The max degree of dummy and nondummy interaction 
        predictor variables.}
   \item{del0cols}{Delete all-0 columns from output.  If TRUE, determine
        such columns during computation; if a character vector, delete
        the colums with those names.}
}

\details{

   The \code{getPoly} function takes in a dataframe or data matrix and
   generates polynomial terms of predictor variables. 
   
   Note the subtleties involving dummy variables.  The square, cubic and
   so on terms are the same as the original variable, and the various
   duplicates must be eliminated.  
   
   Similarly, after dummy variable are created from a categorical
   variable having more than two levels, the resulting columns will be
   orthogonal to each other.  In most cases, the result all-0 columns
   should be eliminated, hence the argument \code{del0cols}.  In typical
   usage, one would set this argument to TRUE at the training stage, and
   then for all predictions use the vector of names in the
   \code{retainedNames} component in the return value of this function;
   \code{predict.polyFit} does the latter automatically, so
   \code{getPoly} must set \code{del0cols} to TRUE in the training
   phase.

   Also, this function treats predictor variables that have only two kinds
   of values categorical variables, so the data needs preprocessing to
   multiple dummy variables if needed.
}
\value{
The return value of \code{getPoly} is a \code{polyMatrix} 
object.  This is an S3 class containing a dataframe \code{
xdata} of the generated polynomial terms. The predictor 
variables have column names V1, V2, etc. The object also contains 
a vector \code{endCols}, whose i-th element is the end column of polynomial 
terms for i-th degree.
}
\examples{
x1 <- 1:3
x2 <- 2:4
x3 <- 3:5
y <- c(1,2,4)
xy1 <- cbind(x1,x2,x3,y)
getPoly(xy1,2)
# output is the poly term matrix xy, plus endCols; for this case, the
# latter is (4,14), meang hat the first-degree terms end with V4 and the
# second-degree ones go from V5 to V14

d1 <- 2:4
d2 <- 7:9
d3 <- c(TRUE, FALSE, FALSE)
d4 <- c(TRUE, TRUE, FALSE)
dy <- c(1,2,4)
xy2 <- cbind(d1,d2,d3,d4,dy)
getPoly(xy2,3)
}
