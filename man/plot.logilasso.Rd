\name{plot.logilasso}
\alias{plot.logilasso}
\alias{plot.predlogilasso}
\alias{plot.predlogispez}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot function}
\description{
Plot function for objects of class \code{logilasso},
\code{predlogilasso} and \code{predlogispez}.
}
\usage{
\method{plot}{logilasso}(x, grenze = 0,...)
\method{plot}{predlogilasso}(x, grenze = 0, ...)
\method{plot}{predlogispez}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{logilasso}, \code{predlogilasso} or
    \code{predlogispez}.}
  \item{grenze}{Up to which fraction of the maximal absolute component
    of beta, the individual components are treated as zero. If 1, then
    all beta components are treated as zero, if 0, then all components
    which do not equal zero are treated as nonzeros.}
  \item{\dots}{Additional arguments to the plot function.}
}
\details{
For objects of class \code{predlogispez} the trajectories of the
components of the beta vectors are plotted against the lambdas. For
objects of class \code{predlogilasso} the corresponding undirected
graphical model (CIG) is plotted.

For objects of the class \code{logilasso}, first the \code{predict()}
function is applied and the resulting object of either class
\code{predlogispez} or \code{predlogilasso} is plotted as described above.
}
\examples{
library(gRbase)
data(reinis)

fit <- logilasso(reinis,lambdainit=1,lambdamin=0.1)
plot(fit)

fitcv <- logilasso(reinis,lambdainit=1,lambdamin=0.1,cvfold=3)
plot(fitcv)

levellogi <- levelcv(reinis,lambdainit=1,lambdamin=0.1,to.which.int=3,cvfold=3)
plot(levellogi)

pred <- predict(fit,lambda=0.3)
plot(pred)
}
\keyword{hplot}

