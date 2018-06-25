\name{xvalPoly,xvalKf,xvalDnet}
\alias{xvalPoly}
\alias{xvalKf}
\alias{xvalDnet}

\title{
Cross validation on Polynomial regression
}
\description{
Separate the dataset to training and test set, and obtain the mean absolute error (for linear regression) or accuracy (for logistic regression).
}
\usage{
xvalPoly(xy, maxDeg, maxInteractDeg = maxDeg, use = "lm", pcaMethod =
   FALSE, pcaPortion = 0.9, glmMethod = "all",nHoldout=10000,yCol=NULL)
xvalKf(xy, nHoldout = min(10000, round(0.2 * nrow(xy))), yCol = NULL, 
    units, activation, dropout) 
xvalDnet(x,y,hidden,output="\\"sigm\\"",numepochs=3,pcaMethod = FALSE, 
    pcaPortion=0.9,scaleXMat=TRUE,nHoldout=min(10000, 
    round(0.2*nrow(x)))) 
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
  \item{nHoldout}{Size of test test.}
  \item{units}{Vector specifying number of units per layer. There is one
     element for each hidden layer, then an element for the dummy layer
     that produces the predictions; specify this as NA.}
  \item{activation}{Vector specifying the activiation functions. One
     element for each hidden layer, then one function producing the
     final prediction, e.g. \code{'linear'} for regression and
     \code{'softmax'} for classification.}
  \item{units}{Vector specifying number of units per layer.}
  \item{x}{Numerical matrix of predictor values.}
  \item{y}{Numerical vector (regression case) or matrix (classification
     case) of response values.  In the classification, if there are
     \code{q} classes then the matrix will have \code{q} columns,
     exactly one 1 per row.}
}

\details{

   The \code{xvalPoly} function divides the data to training and test sets,
   and use \code{polyFit} to generate models using training data and use
   \code{polyFit.predict} to generate results on test set, and compare the
   results. The \code{xvalKf} does the same for \code{kerasformula}
   neural networks package.
}
\value{
The return value of \code{xvalPoly} is an R vector of mean absolute
error (for \code{lm}) or probability of correct classification
(for \code{glm}).  The i-th element of the vector is for degree i.
}

\examples{
y <- mtcars[,1]
data <- cbind(mtcars[,-1], y) # make y column the last column
pf1 <- xvalPoly(data,2,2,"lm",0.8,FALSE)
pf2 <- xvalPoly(data,5,3,"lm",0.8,TRUE,0.9)
}

