\name{prgeng}
\docType{data}
\alias{prgeng}
\title{Silicon Valley programmers and engineers}
\description{
This data set is adapted from the 2000 Census (5\% sample, person records). It is restricted to programmers and engineers in the Silicon Valley area.

The variable codes, e.g. occupational codes, are available from the Census Bureau, at http://www.census.gov/prod/cen2000/doc/pums.pdf. (Short code lists are given in the record layout, but longer ones are in the appendix Code Lists.)
The variables are:

age, with a U(0,1) variate added for jitter

yrentry, year of entry to the U.S. (0 for natives)

wageinc, wage income

wkswrkd, number of weeks worked

categorical variables:

educ: 01-09 code no college; 10-12 means some college; 13 is a bachelor's degree, 14 a master's, 15 a professional degree, and 16 is a doctorate

occ, occupation

birth, place of birth

sex, 1 for male, 2 for female

}
\usage{data(prgeng)}
\source{Public Use Microdata Sample: 2000 Census of Population and Housing.}
\references{
  Census 2000, Public Use Microdata Sample, (PUMS), United States, prepared by the U.S. Census Bureau, 2003. 
}
\keyword{datasets}

\examples{
data(prgeng) # US Census engineer wage data 
xy <- prgeng[, c(1:4,6,5)] # Y value (wageinc here) must be last for polyFit() and FSR()
xy <- toFactors(xy, 2:4)  # convert education, occupation, and sex to factors
                          # (don't pass categorical variables as integer codes when fitting)
}

