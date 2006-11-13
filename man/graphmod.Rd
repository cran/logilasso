\name{graphmod}
\alias{graphmod}
\alias{graphmod.predlogilasso}
\alias{graphmod.cvlogilasso}
\alias{graphmod.levellogilasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to plot a graphical model.}

\description{Function to plot a graphical model for objects of class
\code{predlogilasso}, \code{cvlogilasso} and \code{levellogilasso}. For
details about graphical models, see reference.  }

 \synopsis{graphmod(obj, grenze = 0, nnames = NULL, ...)}

 \usage{
graphmod(obj, grenze = 0, nnames = NULL, ...)

\method{graphmod}{predlogilasso}(obj, grenze = 0, nnames = NULL,... )

\method{graphmod}{cvlogilasso}(obj, grenze=0, nnames = NULL, lambda = NULL, ... )
\method{graphmod}{levellogilasso}(obj, grenze=0, nnames = NULL, lambda = NULL, to.which.int = NULL, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{Object of class \code{predlogilasso}, \code{cvlogilasso} or
    \code{levellogilasso}.}
  \item{grenze}{Components of beta which are smaller than
    \code{grenze} times the biggest component of beta are truncated to
    zero.}
  \item{nnames}{Character vector. Names of the nodes in the graphical
    model.}
  \item{lambda}{For object of class \code{cvlogilasso} or
    \code{levellogilasso}: Should not be specified except if the
    graphical model for a specific value of \code{lambda} should be
    drawn. Otherwise this value is assessed by cross-validation.}
  \item{to.which.int}{For objects of class \code{levellogilasso}: Should
  not be specified except if the graphical model for a specific
  interaction number is desired. Otherwise assessed by cross-validation.}
  \item{\dots}{Additional arguments to be passed to the plotting function.}
}
\references{Steffen L. Lauritzen, Graphical Models, Oxford University Press, 1996}
\examples{
library(gRbase)
data(reinis)

fit <- logilasso(reinis,lambdainit=1)
pred <- predict(fit,lambda=0.5)
# Graphical model of an object of class predlogilasso
graphmod(pred,nnames=as.character(c(1:6)))

fit <- logilasso(reinis,lambdainit=1,cvfold=3)
# Graphical model of an object of class cvlogilasso
graphmod(fit,nnames=c("one","two","three","four","five","six"))}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
