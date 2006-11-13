\name{logilasso-package}
\alias{logilasso-package}
\docType{package}
\title{
logilasso
}
\description{
Fits loglinear models to data arising from sparse contingency tables
by penalizing the multinomial likelihood.}
\details{
\tabular{ll}{
Package: \tab logilasso\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-11-13\cr
Depends: \tab methods, dynamicGraph\cr
License: \tab GPL\cr
Built: \tab R 2.4.0; ; 2006-11-13  09:38:41; unix\cr
}
Please note that this is an *early test release*.

The best entry point for the package are the examples in the help
file of the function \code{\link{logilasso}}.
Index:
\preformatted{
  logilasso                 Fits a loglinear model or/and performs
                            cross-validation
  levelcv                   Performs cross-validation for the specified
                            number of interactions
  traceplot                 Plots the solution path from lambdamax to lambdamin for
                            all components of the solution vector beta
  graphmod                  Plots a graphical model
  plot.logilasso            Plot method for a logilasso object
  predict.logilasso         Predict method for a logilasso object
}
}
\author{
  Corinne Dahinden
  
  Maintainer: Corinne Dahinden <dahinden@stat.math.ethz.ch>

}
\keyword{package}
