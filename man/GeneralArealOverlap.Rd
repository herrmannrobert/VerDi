% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{GeneralArealOverlap}
\alias{GeneralArealOverlap}
\title{General areal overlap}
\usage{
GeneralArealOverlap(data)
}
\arguments{
\item{data}{matrix or data frame including two colums - one for each population's distribution.}
}
\description{
Function calculates interspecific overlap among two populations.
}
\details{
The general areal overlap is defined as the area occupied by both populations devided by 
the area occupied by either one or the other population.
}
\examples{
#Example shows how two populations overlap in one of five water layers

df <- cbind(c(1,1,1,0,0),
            c(0,0,1,1,1))

GeneralArealOverlap(df)
}
