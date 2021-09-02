\name{prgeng}
\alias{prgeng}
\alias{peDumms}
\alias{peFactors}
\alias{pef}

\title{
Silicon Valley programmers and engineers data
}
\usage{
data(prgeng)
data(peDumms)
data(peFactors)
}

\description{
This data set is adapted from the 2000 Census (5\% sample, person
records).  It is mainly restricted to programmers and engineers in the 
Silicon Valley area.  (Apparently due to errors, there are some from
other ZIP codes.)

There are three versions:

   \itemize{

    \item{\code{prgeng}, the original data, with categorical variables,
    e.g. Occupation, in their original codes}

    \item{\code{peDumms}, same but with categorical variables
    converted to dummies; due to the large number of levels the birth
    and PUMA data is not included}

    \item{\code{peFactors}, same but with categorical variables
    converted to factors}

    \item{\code{pef}, same as \code{peFactors}, but having only columns
    for age, education, occupation, gender, wage income and weeks
    worked.  The education column has been collapsed to Master's degree,
    PhD and other.}

    }

The variable codes, e.g. occupational codes, are available from the Census
Bureau, at
\url{http://www.census.gov/prod/cen2000/doc/pums.pdf}.
(Short code lists are given in the record layout, but longer ones are in
the appendix Code Lists.)

The variables are:

   \itemize{

   \item{\code{age}, with a U(0,1) variate added for jitter}

   \item{\code{cit}, citizenship; 1-4 code various categories of
   citizens; 5 means noncitizen (including permanent residents)}

   \item{\code{educ}: 01-09 code no college; 10-12 means some college;
   13 is a bachelor's degree, 14 a master's, 15 a professional degree and
   16 is a doctorate}

   \item{\code{occ}, occupation}

   \item{\code{birth}, place of birth}

   \item{\code{wageinc}, wage income}

   \item{\code{wkswrkd}, number of weeks worked}

   \item{\code{yrentry}, year of entry to the U.S. (0 for natives)}

   \item{\code{powpuma}, location of work} 

   \item{\code{gender}, 1 for male, 2 for female}

   }

}


