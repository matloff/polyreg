\name{misc}
\alias{toFactors}
\alias{toDummies}

\title{Miscellaneous}

\description{Utilities}

\usage{
   toFactors(df,cols)
   toDummies(df,col)
}

\arguments{
  \item{df}{A data frame.}
  \item{cols}{A vector of column numbers.}
  \item{col}{A column number.}
}

\details{

   The \code{toFactors} function converts each \code{df} column in 
   \code{cols} to a factor, returns new version of \code{df}.

   The \code{toDummies} function converts the \code{df} column 
   \code{col} to a matrix of dummy variables, returns new version of 
   \code{df}, with the original column removed.

}

