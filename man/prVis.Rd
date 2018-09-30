\name{prVis}
\alias{prVis}
\alias{addRowNums}

\title{Polynomial-Based Manifold Exploration }

\description{
Polynomial-based alternative to t-SNE, UMAP etc.}

\usage{
prVis(xy, labels = FALSE, deg = 2, scale = FALSE, nSubSam = 0, 
    nIntervals = NULL, saveOutputs = FALSE, cex = 0.5)
addRowNums(np, savedPrVisOut) 
}

\arguments{
  \item{xy}{Data frame with labels, if any, in the last column.}
  \item{labels}{If TRUE, have class labels.}
  \item{deg}{Degree of polynomial.}
  \item{scale}{If TRUE, call \code{scale} on nonlabels data before 
     generating polynomial terms.}
  \item{nSubSam}{Number of random rows of \code{xy} to sample; 0 means
     use the full dataset.} 
  \item{nIntervals}{If labels column is continuous, discretize into this
     many levels.}
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
