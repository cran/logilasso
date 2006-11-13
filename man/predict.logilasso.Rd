\name{predict.logilasso}
\alias{predict.logilasso}
\alias{predict.cvlogilasso}
\alias{predict.levellogilasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predicts the interaction vector beta of a loglinear interaction
  model.}
\description{
Predicts the interaction vector(s) beta of a loglinear model
\code{log(p)=X*beta}, fitted either
by \code{logilasso} or \code{levelcv}. If \code{lambda} is specified,
this lambda is taken to predict the beta for this value of lambda. If no
value for \code{lambda} is specified, then the optimal value calculated
by cross-validation is taken for objects of class
\code{cvlogilasso}. For objects of class \code{logilasso} where no
cross-validation was performed, the whole solution path for all lambdas
is returned.
}
\usage{
\method{predict}{logilasso}(object, lambda = NULL, ...)
\method{predict}{cvlogilasso}(object, lambda = NULL, ...)
\method{predict}{levellogilasso}(object, lambda = NULL, to.which.int =
NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{levellogilasso}, \code{cvlogilasso} or class \code{logilasso}. }
  \item{lambda}{Value for the penalization parameter \code{lambda},
    for which the corresponding beta should be calculated.}
  \item{to.which.int}{The number of factors the model should be
    predicted for.}
  \item{...}{Additional arguments to predict function.}
}

\value{ Is either an object of class \code{predlogilasso} if a value for
  \code{lambda} was specified or if the optimal lambda can was
  assessed by cross-validation. Otherwise, if no value for \code{lambda}
  was specified and at the same time \code{cvfold} was chosen to be 1
  (no cross-validation) it is of class \code{predlogispez}. The
  difference between these two classes is described below.
  
  \item{beta}{A predicted value for beta if the object is of class
    \code{predlogilasso}.  For the class \code{predlogispez} this is a
    matrix consisting of the columns \code{beta} for the whole solution
    path.}

\item{lambda}{The lambda(s) for which the beta(s) was/were calculated.}
  \item{probs}{Probabilities according to the model probs=exp(X*beta)}
  \item{nls}{Negative likelihood score. For details see \url{http://stat.ethz.ch/~dahinden/Paper/BMC.pdf}}
  \item{betapath}{The whole solution originally calculated path. For objects of class
    \code{predlogispez} this equals \code{beta}.}
  \item{lambdapath}{The \code{lambdas} corresponding to the value \code{betapath}.}
  \item{losspath}{\code{nls} for the whole solution path.}  }

\seealso{\code{\link{logilasso}}}
\examples{
library(gRbase)
data(reinis)

fit <- logilasso(reinis,lambdainit=1,lambdamin=0.1)
pred1 <- predict(fit,lambda=0.5)
pred2 <- predict(fit)

fitcv <- logilasso(reinis,lambdainit=1,lambdamin=0.1,cvfold=3)
predcv1 <- predict(fitcv)

levellogi <- levelcv(reinis,lambdainit=1,lambdamin=0.1,to.which.int=3,cvfold=3)
predlevel <- predict(levellogi)

## Methods plot and graphmod exist for all predicted models
## Except for pred2, there is no graphmod method, because no
## lambda was specified

plot(predcv1)
graphmod(predcv1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

