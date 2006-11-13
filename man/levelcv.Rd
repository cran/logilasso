\name{levelcv}
\alias{levelcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-Validation for different interaction levels}
\description{
  Performs cross-validation for interactions involving \code{k,
    k+1,...,k+u} factors which can be specified by the arguments
  \code{from.which.int=k} and \code{to.which.int=k+u}. For each
  interaction order \code{k+j} a \code{cvfold} cross-validation is
  performed. Returns an object of class \code{levellogilasso} which is a
  subclass of \code{cvlogilasso}.  }

\usage{levelcv(Y ,combX = NULL , from.which.int = 1, to.which.int=NULL,
epsilon = 0.1, lambdainit =1, lambdamin = 0.1, trace = 1,
method = "groupl1", cvfold=10, Newton = TRUE, stopkrit = 1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  The function can be initialized in 2 different ways:
  1. \code{Y} as a contingency table.
  2. \code{Y} together wit a combination matrix \code{combX} (see example).

  \item{Y}{The vector containing the counts for each cell in the
 contingency table. Either given as table or as vector together with a
 matrix \code{combX}, where each component of the vector corresponds
 to a combination of factors in the corresponding row of
 \code{combX}.}

\item{combX}{Matrix of dimension length(Y) x number of factors. For
 each component of \code{Y} the corresponding combinations of the
 level is given as a row entry. See example.}

  \item{from.which.int}{Numeric. The samllest number of factors whose interaction
                        is considered. For example if it is =1, the
			algorithm starts considering all main effects models. }
  \item{to.which.int}{Up to which interaction should the solution be
    calculated. to.which.int=n indicate that interaction involving n
    factors are considered.}
  \item{epsilon}{The step length for lambda. In each step, lambda
    diminished by epsilon.}
  \item{lambdainit}{The upper bound for lambda, where the solution path
    for beta starts.}
  \item{lambdamin}{The lower bound for lambda, where the solution path ends.}
  \item{trace}{Defines what is printed out during the calculation.
    0 = nothing,
    1 = Current cvfold 
    2 = additionally points for each lambda
    3 = Every 10th step writes the active set
    4 = Additionally the new active set are written out, whenever a
        component enters the active set
      }
  \item{method}{Is either "groupl1", "l1" or "l2", depending on the
                desired penalization of the coefficients.}
  \item{cvfold}{If cvfold is larger than 1, cross-validation is performed.}
  \item{Newton}{Logical. If \code{Newton=TRUE}, Newton steps are
    performed, otherwise the function \code{optim} is used.}
  \item{stopkrit}{Convergence tolerance; the smaller the more precise.}
}
\details{
  A \code{levellogilasso} object is returned. This is a list of objects
  of class \code{cvlogilasso} for each interaction from
  \code{from.which.int} to \code{to.which.int}. For objects of that
  class, the methods \code{plot}, \code{predict}, \code{traceplot} and
  \code{graphmod} are available.

}
\value{
  The elements of the list
  of class \code{cvlogilasso} each consists of the following entries:
  \item{loss}{A loss matrix of dimension cvfold x number of assessed
lambdas in the solution
  path. For each part of the data left out, the loss for all lambdas of
  the solution path is calculated. The \code{loss} is \code{NULL} in
  case of \code{cvfold=1}}
  \item{path}{A matrix. In the second row the newly active or newly
    inactive components of beta are listed to the corresponding lambda
    in the first row of the matrix. Is \code{NULL} for the class \code{cvlogilasso}}
  \item{betapath}{A matrix of dimension length(beta) x lambda in the
    solution path. The columns consist of the betas for the different
    lambdas in \code{lambdapath}. Is \code{NULL} for the class
    \code{cvlogilasso}.}
  \item{lambdapath}{A vector of all lambdas for which the solution was
    calculated.}
  \item{X}{The design matrix used to fit the log-linear model.}

  \item{nrfac} {Number of factors.}  
}
\references{Corinne Dahinden, Giovanni Parmigiani, Mark Emerik and Peter
Buehlmann available at \url{http://stat.ethz.ch/~dahinden/Paper/BMC.pdf}}

\examples{
library(gRbase)
data(reinis)
levellogi <- levelcv(reinis,lambdainit=1,lambdamin=0.1,from.which.int=1,to.which.int=3,cvfold=3)

## The methods predict, plot, traceplot and graphmod are defined for the resulting
## object of class cvlogilasso
predlevel <- predict(levellogi)
plot(predlevel)
traceplot(predlevel)
graphmod(predlevel)

### Different initialization: Y and combX
### 5 factors: All have 2 levels
### nlev would be c(2,2,2,2,2)
Y     <- c(4,1,3,2,9)
combX <- rbind(c(1,0,1,1,0),c(1,0,0,1,1),c(0,1,0,0,1),c(0,0,1,0,0),c(1,1,0,0,1))
### 4 observations wit level 1 of factor 1, level 0 of factor two, level 1 of factor
### 3 and so on.
levellogi2 <- levelcv(Y,combX=combX,from.which.int=1,to.which.int=3,cvfold=3)
}
\keyword{models}

