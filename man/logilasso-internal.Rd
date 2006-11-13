\name{logilasso-internal}
\alias{f.createX}
\alias{f.rbindspez}
\alias{f.xy}
\alias{f.nabla}
\alias{f.nabla2}
\alias{f.fn}
\alias{f.gradient2}
\alias{f.to.minimize2}
\alias{f.group}
\alias{f.newtonstep}
\alias{f.nabla2active}
\alias{f.adplot}
\alias{graphmod.default}
\alias{traceplot.default}
\title{Internal logilasso functions}
\description{Internal logilasso functions}
\usage{
f.createX(nlev, to.which.int = NULL )
f.rbindspez(x,y)
f.xy(y,X,to.which.int=NULL)
f.nabla(tXY,X,expXbeta,n)
f.nabla2(Xactive,expXbeta,n)
f.fn(betactive, n, tXactive, Xactive, tXY, lambda, anzahl , tester,
     methode , active)
f.gradient2(beta,n,tXactive,Xactive,tXY,lambda,anzahl,tester,
            methode,active,...)
f.to.minimize2(Xbeta, Y, n, expXbeta)
f.group(x,anzahl,tester,pen0=FALSE,summieren=TRUE)
f.newtonstep(betactive, nabla2active, nablactive, lambda,
             anzahlactive, testera, method)
f.nabla2active(nabla2, betactive, lambda, anzahlactive, tester, method)
f.adplot(adjac,nnames=NULL, directed = FALSE,...)

}
\details{These are not to be called by the user.}
\keyword{internal}
