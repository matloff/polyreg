\name{misc}
\alias{toFactors}

\title{Miscellaneous}

\description{Utilities}

\usage{
   toFactors(df,cols)
}

\arguments{
  \item{df}{A data frame.}
  \item{cols}{A vector of column numbers.}
}

\details{

   The \code{toFactors} function converts each \code{df} column in 
   \code{cols} to a factor, returns new version of \code{df}. 
   Should be used on categorical variables stored as integer codes 
   before calling the library's main functions, 
   including \code{getPoly}, \code{FSR}, or \code{polyFit}.

}

