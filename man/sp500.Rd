% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sp500}
\alias{sp500}
\title{S&P 500 daily log-returns}
\format{
An object of class \code{"xts"}.
}
\source{
Realized library of the Oxford-Man Institute
}
\usage{
data(sp500)
}
\description{
Daily data on S&P 500 collected from the realized library of
the Oxford-Man Institute \insertCite{heber_2009}{rumidas}.
}
\examples{
head(sp500)
summary(sp500)
plot(sp500)
}
\references{
\insertAllCited{}
}
\keyword{datasets}
