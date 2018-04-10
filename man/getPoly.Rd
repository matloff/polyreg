\name{getPoly}
\alias{getPoly}
\alias{polyMatrix}

\title{Get polynomial terms}

\description{
Generate polynomial terms of predictor variables for a dataframe or data matrix.
}

\usage{
getPoly(xydata, deg, maxInteractDeg = deg)
}

\arguments{
  \item{xydata}{Data matrix or dataframe without response variable.}
  \item{deg}{The max degree to generate polynomial terms.}
  \item{maxInteractDeg}{The max degree of dummny and nondummy interaction predictor variables.}
}

\details{

   The \code{getPoly} function takes in a dataframe or data matrix and
   generates polynomial terms of predictor variables, and it gets rid of
   redundant terms for categorical variables. For example, for a predictor
   variable with only 0 and 1, its square and cubic terms are the same.
   Also, this function treats predictor variables that have only two kinds
   of values categorical variables, so the data needs preprocessing to
   multiple dummy variables if needed.
}
\value{
The return value of \code{getPoly} is a \code{polyMatrix} 
object.  This is an S3 class containing a dataframe of the generated 
polynomial terms. The predictor variables have column 
names V1, V2, etc. The object also contains a list, whose i-th 
element is the end column of polynomial terms for i-th degree.
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
