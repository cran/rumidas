% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/midas_functions.R
\name{MEM_X_loglik}
\alias{MEM_X_loglik}
\title{MEM-X log-likelihood (with skewness parameter)}
\usage{
MEM_X_loglik(param, x, daily_ret, z)
}
\arguments{
\item{param}{Vector of starting values.}

\item{x}{Dependent variable, usually the realized volatility. It must be positive and "xts" object.}

\item{daily_ret}{Daily returns, which must be an "xts" object, and with the same length of x.}

\item{z}{Additional daily variable which must be an "xts" object, and with the same length of x.}
}
\value{
The resulting vector is the log-likelihood value for each \eqn{i,t}.
}
\description{
Obtains the log-likelihood of the base MEM, with an asymmetric term linked to past negative returns and an additional
X part (for instance, the VIX).
}
\examples{
\donttest{
start_val<-c(alpha=0.10,beta=0.8,gamma=0.05,delta=0.01)
real<-(rv5['2009/2010'])^0.5		# realized volatility
r_t<-sp500['2009/2010']
z<-vix['2009/2010']
sum(MEM_X_loglik(start_val,real,r_t,z))
}
}
\references{
\insertAllCited{}
}
\keyword{internal}
