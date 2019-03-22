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
   Should be used on categorical variables stored as integer codes 
   before calling the library's main functions, 
   including \code{getPoly}, \code{FSR}, or \code{polyFit}.

   The \code{toDummies} function converts the \code{df} column 
   \code{col} to a matrix of dummy variables, returns new version of 
   \code{df}, with the original column removed. Convenience function.
   Should NOT be used on categorical variables with more than two levels 
   in conjunction with the library's main functions, 
   including \code{getPoly}, \code{FSR}, or \code{polyFit}.

}

