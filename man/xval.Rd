\name{polyFit}
\alias{polyFit}

\title{Crossvalidate a Polynomial Regression}

\description{
Crossvalidate a polynomial regression using a linear or logistic model; convenience function for polyFit().
}

\usage{
xvalPoly(reps=5, mean_loss=function(y, yhat){mean(abs(y - yhat))}, stat=median, pTesting=0.1, seed=NULL, ...)
}

\arguments{
\item{reps}{Number of cross-validations to perform (separate cross-validations, not k-folds). Default: 5.}
\item{mean_loss}{Function for average out-of-sample loss. Default: mean absolute error (MAE).}
\item{stat}{Function defining the statistic that will summarize mean_loss across reps. Default: median.}
\item{pTesting}{Proportion of holdout data. Default: 0.1.}
\item{seed}{Optional seed. Default: seed not set by this function.}
\item{...}{Additional arguments for polyFit().}
}

\examples{
N <- 125
xy <- data.frame(x1 = rnorm(N), 
                 x2 = rnorm(N),
                 group = sample(letters[1:5], N, replace=TRUE),
                 score = sample(100, N, replace = TRUE) # final column is y
                 )
xvalPoly(xy=xy, deg=2, noisy=FALSE)                      
                      
}



