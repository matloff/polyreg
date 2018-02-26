\name{polyFit}
\alias{polyFit}

\title{Polynomial Fit}

\description{
Fit a polynomial regression/classification model using linear 
or logistic regreesion. 
}

\usage{
polyFit(xy, maxDeg, maxInteractDeg = maxDeg, use = "lm", trnProp = 0.8,
        pcaMethod = FALSE, pcaPortion = 0.9, glmMethod = "all")
}

\arguments{
  \item{xy}{Data matrix or dataframe with response variable in the last column.}
  \item{maxDeg}{Max degree for power terms.}
  \item{maxInteractDeg}{Max degree of interaction terms.}
  \item{use}{Can be "lm" for using linear regreesion, and "glm" for 
              using logistic regression.}
  \item{trnProp}{Proportion of data to be used as training set.}
  \item{pcaMethod}{If TRUE, use PCA.}
  \item{pcaPortion}{If pcaMethod is TRUE, use components up to this
                    proportion of total variance.}
  \item{glmMethod}{For classification problems.  If there are more than
                   two classes, this can be "all" for All-vs-All method,
                   or "one" for One-vs-All method.}
}

\details{

   The \code{polyFit} function takes in a dataframe and uses the \code{getPoly}
   function to generate polynomial terms of the predictor variables. 
   Then it fits the generated data to linear or logistic regression models. If
   \code{pcaMethod == TRUE}, the principal component analysis is performed
   using the indicated portion of principal components before fitting the
   models. When logistic regression for classification is indicated,
   All-vs-All or One-vs-All methods can be applied to deal with
   multiclass problem.
}
\value{
The return value of \code{polyFit} is an R vector of mean absolute 
error (for \code{lm}) or probability of correct classification
(for \code{glm}).  The i-th element of the vector is for degree i.
}
\examples{
y <- mtcars[,1]
data <- cbind(mtcars[,-1], y) # make y column the last column
pf1 <- polyFit(data,2,2,"lm",0.8,FALSE)
pf2 <- polyFit(data,5,3,"lm",0.8,TRUE,0.9)

}
