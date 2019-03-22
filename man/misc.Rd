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

\examples{
data(prgeng) # US Census engineer wage data 
xy <- prgeng[, c(1:4,6,5)] # Y value (wageinc here) must be last for polyFit() and FSR()
xy <- toFactors(xy, 2:4)  # convert education, occupation, and sex to factors
                          # (don't pass categorical variables as integer codes when fitting)
}
