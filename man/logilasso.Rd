\name{logilasso}
\alias{logilasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fits a log-linear model and/or performs cross-validation.}

\description{Fits a log-linear interaction model (log(p)=X*beta)
  penalizing the log-likelihood function assuming a multinomial
  sampling scheme. Penalization can be chosen as either an l1, l2 or a
  group l1 penalty. In addition cross-validation is performed if
  cvfold is chosen larger than 1. The returned objects are of class
  \code{logilasso} if no cross-validation has been performed or of
  class \code{cvlogilasso}, which is a subclass of \code{logilasso},
  if cross-validation has been performed.}
  

\usage{logilasso(Y, combX = NULL, to.which.int=NULL, epsilon = 0.1,
lambdainit = 1, lambdamin = 0.1, trace = 1, method = "groupl1",
cvfold=1, Newton = TRUE, stopkrit = 1e-8, sse = NULL, X = NULL)}

%- maybe also 'usage' for other objects documented here.

\arguments{ There are 2 different ways of initializing the function:
1. \code{Y} as a contingency table.  2. \code{Y} together wit a
combination matrix \code{combX} (see example).

\item{Y}{The vector containing the counts for each cell in the
 contingency table. Either given as table or as vector together with a
 matrix \code{combX}, where each component of the vector corresponds
 to a combination of factors in the corresponding row of
 \code{combX}.}

\item{combX}{Matrix of dimension length(Y) x number of factors. For
 each component of \code{Y} the corresponding combinations of the
 level is given as a row entry. See example.}

\item{to.which.int}{Up to which interaction should the solution be
 calculated. to.which.int=n indicate that interaction involving n
 factors are considered.}

\item{epsilon}{The step length of lambda. In each step the
 penalization parameter lambda is decreased by epsilon.}

\item{lambdainit}{The upper bound for lambda, where the solution path
 for beta starts.}

\item{lambdamin}{The lower bound for lambda, where the solution path ends.}

\item{trace}{Defines what is printed out during the calculation.
 0 = nothing,
 1 = Current cvfold 
 2 = additionally points for each lambda
 3 = Every 10th step writes the active set
 4 = Additionally the new active set are written out, whenever a
     component enters the active set}

\item{method}{Is either "groupl1", "l1" or "l2", depending on the penalization of 
 the coefficients.}

\item{cvfold}{If cvfold is larger than 1, cross-validation is performed.}

\item{Newton}{Logical. If \code{Newton=TRUE}, Newton steps are
 performed, otherwise the function \code{optim} is used.}

\item{stopkrit}{Convergence tolerance; the smaller the more precise, see
 details below.}

\item{sse} {.Random.seed vector.} 

\item{X}{The design matrix \code{X} which is used for fitting a
log-linear model. Should not be specified by the user.}
}

\details{ For the convergence criteria see chapter 8.2.3.2 of Gill et
 al. (1981).  \emph{Practical Optimization}, Academic Press.
  
  Dimitri P. Bertsekas (2003) \emph{Nonlinear Programming}, Athena
  Scientific.

  For the resulting objects of class \code{logilasso} and
  \code{cvlogilasso} the methods \code{plot}, \code{predict} and
  \code{traceplot} are available. If Cross-Validation was performed
  (\code{cvfold>1}), then in addition the method \code{graphmod} is
  applicable, which plots a graphical model.}

\value{A \code{cvlogilasso}
  object is returned in case cvfold is larger
  than 1. A \code{logilasso} object is returned in case cvfold is equal
  to 1. The class \code{cvlogilasso} is a subclass of
  \code{logilasso}.

\item{loss}{A loss matrix of dimension cvfold x number of assessed
 lambdas in the solution path. For each part of the data left out, the
 loss for all lambdas of the solution path is calculated. The
 \code{loss} is \code{NULL} in case of \code{cvfold=1}}

\item{path}{A matrix. In the second row the newly active or newly
 inactive components of beta are listed to the corresponding lambda in
 the first row of the matrix. Is \code{NULL} for the class
 \code{cvlogilasso}}

\item{betapath}{A matrix of dimension length(beta) x number of
 assessed lambda in the solution path. The columns consist of the
 betas for the different lambdas in \code{lambdapath}. Is \code{NULL}
 for the class \code{cvlogilasso}.}

\item{lambdapath}{A vector of all lambdas for which the solution was
 calculated.}

\item{X}{The design matrix used to fit the log-linear model.}

\item{nrfac} {Number of factors.}  
}
\references{Corinne Dahinden, Giovanni Parmigiani, Mark Emerik and Peter
Buehlmann available at \url{http://stat.ethz.ch/~dahinden/Paper/BMC.pdf}}
\author{Corinne Dahinden, \email{dahinden@stat.math.ethz.ch} }
\examples{
## Use logilasso on the reinis dataset provided in the
## package gRbase
library(gRbase)
data(reinis)

## Fit a log-linear model for lambdas between 1 and 0.1
## No cross-validation is performed
fit <- logilasso(reinis,lambdainit=1,lambdamin=0.1)

### Different initialization: Y and combX
### 5 factors: All have 2 levels
Y     <- c(4,1,3,2,9)
combX <- rbind(c(1,0,1,1,0),c(1,0,0,1,1),c(0,1,0,0,1),c(0,0,1,0,0),c(1,1,0,0,1))
### 4 observations wit level 1 of factor 1, level 0 of factor two, level 1 of factor
### 3 and so on.
### The rows of combX correspond to a the levels of the five factors
### Must be numeric with 0/1/2/... and so on.
fit2 <- logilasso(Y,combX)

## Trajectories from lambdainit to lambda optimal
plot(fit2)
traceplot(fit2)

## Predict functions
pred <- predict(fit,lambda=0.3)

## Perform 3-fold cross-validation
fitcv <- logilasso(reinis,lambdainit=1,lambdamin=0.1,cvfold=3)

## Plots a graphical model with the lambda calculated by cross-validation.
plot(fitcv)
graphmod(fitcv)
  
predcv <- predict(fitcv)
}
\keyword{models}
\keyword{regression}
