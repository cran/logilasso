\name{traceplot.logilasso}
\alias{traceplot}
\alias{traceplot.logilasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the trajectories.}
\description{
The trace of all components of beta is plotted against all lambdas which
were considered in the solution path. The function can be applied to an object of
class \code{logilasso}. This object can be of various subclasses such as
\code{cvlogilasso}, \code{predlogilasso} or
\code{predlogilassospez}. Possibly the solution path only consists of
one lambda, then the components of the corresponding beta are plotted.
}
\usage{
traceplot(beta, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{An object of class \code{logilasso}. }
  \item{\dots}{Additional arguments to the plot function.}
}
\examples{
library(gRbase)
data(reinis)

fit <- logilasso(reinis,lambdainit=,lambdamin=0.1)
traceplot(fit)

fit2 <- logilasso(reinis,lambdainit=1,lambdamin=0.1,cvfold=3)
traceplot(fit2)

fit3 <- levelcv(reinis,lambdainit=1,lambdamin=0.1,to.which.int=3,cvfold=3)
traceplot(fit3)

traceplot(predict(fit3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
