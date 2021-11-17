\name{LRfit.object}
\alias{LRfit.object}

\title{
  Estimated while-alive loss rates
}
\description{
This class of objects is returned by the \code{LRfit} functions.
Objects of this class have methods for the functions \code{print},
\code{summary}, and \code{plot}.
}
\value{
All numerical results of an object \code{obj}
are contained in \code{obj$content}, which
is organized in a two-dimensional array with each column containing the results
for one group. Below are the row variables.
\item{t}{A vector of follow-up times \eqn{\tau}.}
\item{llr}{A vector of log-loss rate estimates at \code{t}.}
\item{se_llr}{A vector of standard error estimates for the log-loss rates
in \code{llr}.}
\item{lmuR}{A vector of log-raw cumulative loss estimates at \code{t}.}
\item{se_lmuR}{A vector of standard error estimates for the log-raw
cumulative losses in \code{lmuR}.}
\item{lmuD}{A vector of log-RMST estimates at \code{t}.}
\item{se_lmuD}{A vector of standard error estimates for the log-RMST in \code{lmuD}.}
\item{St}{A vector of overall survival probabilities at \code{t}.}
\item{...}{}
}

\seealso{
\code{\link{LRfit}}, \code{\link{summary.LRfit}}, \code{\link{plot.LRfit}}.
}
\keyword{LRfit}
