#### 3 Mal predictfunktionen logilasso und cvlogilasso  ####
predict.logilasso <- function(object , lambda = NULL, ...)
{
  ## Purpose: Takes an output of class logilasso from logilasso and computes 
  ##          the corresponding beta for the desired lambda. If no lambda
  ##          is specified, beta and probs for all lambdas are computed.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##
  ## object: Output from logilasso
  ## lambda: Attention! Absolute lambda!   
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date:  6 Feb 2006, 14:15

  ##klasse is an additional class, in case lambda is not specified
  klasse <- NULL
  if(is.null(lambda))
    {
      lambda <- object$lambda
      klasse <- c("predlogispez","logilasso")
    }
  else
    klasse <- c("predlogilasso","logilasso")
     
  lambdaf         <- object$lambdapath
  beta            <- object$betapath
  X               <- object$X

  if((length(lambda)==1) && (max(lambdaf)<lambda | min(lambdaf)>lambda))
    warning("Lambda is outside range which was considered in the fit.")
  
  ii    <- which(abs(lambdaf-lambda)==min(abs(lambdaf-lambda)))
  betar <- beta[,ii,drop=FALSE]
  probs <- exp(X%*%betar)

  ret <- list(beta = betar , lambda = lambdaf[ii] , probs = probs ,
              X = X, combX = object$combX, Y = object$Y, 
              corresp = object$corresp ,
              nrfac = object$nrfac , to.which.int = object$to.which.int, 
              method = object$method , stopkrit = object$stopkrit,
              epsilon = object$epsilon, Newton = object$Newton, nls = NULL,
              betapath = object$beta, lambdapath = lambdaf,
              losspath = NULL)
  class(ret) <- klasse
  ret
}

predict.levellogilasso <- function(object, lambda = NULL, to.which.int = NULL, ...)
  {
    ## Purpose: Takes the output from levellogilasso, and computes  
    ##          the "best" beta and the corresponding lambda as well 
    ##          as the best negative likelihood score nls.
    ## 
    ## ----------------------------------------------------------------------
    ## Arguments: An object.
    ## ----------------------------------------------------------------------
    ## Author: Corinne Dahinden, Date: 20 Jan 2006, 17:28

    ##if no cross-validation was performed,levelcv makes no sense
    ##if no cv and no lambda and to.which.int is specified then brake

    if(object[[1]]$cvfold==1)
       {
         if(is.null(to.which.int))
           stop("The function levelcv makes no sense if no cv is performed.")

         else
           {
             warning("The function levelcv makes no sense if no cv is performed.")
             i <- which(unlist(lapply(object,function(x) x$to.which.int))==to.which.int)
           }

       }
    else
      {
        ret <- lapply(object,function(x) min(apply(x$loss,2,mean)))
        i   <- which.min(ret)
      }
    
    object <- object[[i]]

    ret <- predict(object,lambda, ...)
    ret
  }

predict.cvlogilasso <- function(object,lambda = NULL, ...)
  {
    lambdam<- object$lambda
    loss   <- object$loss
    X      <- object$X

    if(is.null(lambda))
      {
        ##remove all columns from loss, which equal the element 1,1. These lambda
        ##are chosen too big!
        indize <- !apply(loss,2,function(x,y) any(x==y),loss[1,1])
                                        #if length of lambdam is 1, then cv made no sense, as it was only evaluated
                                        #at one point... However, we calculate loss all the same without dropping.
        if(length(lambdam)!=1)
          {
            loss   <- loss[,indize,drop=FALSE]
            lambdam <- lambdam[indize]
          }
        
        if(length(loss)==0)
          stop("Lambda was initially chosen too big. Decrease lambda.")
        meanloss <- apply(loss,2,mean)
        ##Negative Likelihood score as described in the paper
        nls      <- min(meanloss)
        ind.max  <- which.min(meanloss)

        lambdafinal <- lambdam[ind.max]
      }
    else
      lambdafinal <- lambda
    #Otherwise error, because lambda suits to specify lambdamin and lambdainit
    object$lambda <- NULL
    #Otherwise error, because X and combX are specified
    combX        <- object$combX
    object$combX  <- NULL
    object$cvfold <- NULL

    topred <- object[names(object)%in%names(formals(logilasso))]
    
    ##To make sure the solution for lambdafinal is computed, we compute the
    ##solution to lambdafinal-0.1
    defaultwarn <- options()$warn
    options(warn=-1)
    rest  <- do.call(logilasso,c(topred,list(lambdamin=lambdafinal, 
                                   lambdainit=(lambdafinal+topred$epsilon))))
    options(warn=defaultwarn)
    predi <- predict(rest,lambdafinal, ...)
    beta  <- predi$beta
    probs <- exp(X%*%beta)

    if(is.null(lambda))
      lambdapath <- c(lambdam[(ind.max-1):ind.max])
    else
      {
        lambdapath <- c(lambdafinal+topred$epsilon,lambdafinal)
        meanloss <- NULL
      }

    ret <- list(beta=beta,lambda=lambdafinal,probs=probs,
                X=object$X, combX=combX, Y=object$Y,
                corresp = object$corresp,
                nrfac = object$nrfac, to.which.int=object$to.which.int,
                method = object$method , stopkrit=object$stopkrit,
                epsilon=object$epsilon, Newton = object$Newton,
                nls = nls, betapath = rest$betapath, lambdapath = lambdapath,
                losspath = meanloss)
    class(ret) <- c("predlogilasso","logilasso")
    ret    
  }

  
traceplot <- function(beta,...)
  UseMethod("traceplot")

traceplot.default <- function(beta,lambda, ...)
  
  ## Purpose: Plots the trace of all components of beta for the lambdas.
  ## ----------------------------------------------------------------------
  ## Arguments: Object of class logilasso, predlogilasso, 
  ##            cvlogilasso or just a matrix of the path.
  ## lambda: If obj is a matrix, then lambda can be specified as x axis
  ## ...: For plot function
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date:  4 Apr 2006, 14:38
  
  {
    if(!exists("xlab"))
      xlab <- "Lambda"
    if(!exists("ylab"))
      ylab <- "Different Traces"
      
    beta <- as.matrix(beta)
    beta <- beta[-1,,drop=FALSE]
    maxi <- max(apply(beta,2,max))
    mini <- min(apply(beta,2,min))
    
    lambda <- lambda[1:ncol(beta)]
    plot(0,xlim=c(min(lambda),max(lambda)),ylim=c(mini,maxi),type="n",
         xlab=xlab,ylab=ylab,...)
    if(ncol(beta)==1)
      {
        for(i in 1:nrow(beta))
          points(lambda,beta[i],col=i)
      }
    for(i in 1:nrow(beta))
      lines(lambda,beta[i,],col=i)

    text(min(lambda),beta[,ncol(beta)],as.character(2:(nrow(beta)+1)))
  }

traceplot.logilasso <- function(beta,...)
  {
    if(!any(class(beta)%in%c("predlogilasso","predlogispez")))
      obj    <- predict(beta, ...)
    else
      obj <- beta
    lambda <- obj$lambdapath
    beta   <- obj$betapath
    traceplot.default(beta,lambda,...)
  }

plot.logilasso <- function( x, grenze = 0,...)
{
  ## Purpose: Zuerst predict, dann entweder traceplot oder graphmod, falls von
  ##          cvlogilasso kommt.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: Of class logilasso
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date: 29 May 2006, 20:29

  ## Entweder predlogilasso oder predlogispez
  betavec <- predict(x)

  if(any(class(betavec)=="predlogispez"))
    traceplot.default(betavec$beta,betavec$lambda,...)
  else
    plot(betavec,grenze,...)

}

plot.predlogispez <- function(x, ...)
  traceplot.default(x$beta,x$lambda,...)


plot.predlogilasso <- function(x,grenze=0,...)
  {
    ## grenze ist bis da wo man das beta noch als ungleich 0 anschaut im Verhältnis zum
    ## maximalen Beta Koeffizient.
    if(!("predlogilasso"%in%class(x)))
     stop("No applicable method for object of this class.")
    
    beta         <- x$beta
    X            <- x$X
    nrfac        <- x$nrfac
    to.which.int <- x$to.which.int

    admat <- matrix(0,nrfac,nrfac)

    tester    <- attr(X,"assign")
    nlev      <- c(table(tester[which(tester%in%c(2:(nrfac+1)))]))+1
    model     <- f.createX(nlev,to.which.int)$model

    tot <- max(abs(beta[-1]))

    for(i in which(abs(beta)>(grenze*tot)))
      {
        if(is.null(i))
          break
        einer <- which(model[,i]==1)
        if(length(einer)>1)
          {
            for(j in 1:(length(einer)-1))
              {
                rest <- einer[(j+1):length(einer)]
                for(k in rest)
                  admat[einer[j],k] <- 1
              }
          }
      }
    admattot <- admat+t(admat)
    f.adplot(admattot,...)
  }

graphmod <- function(obj, grenze=0, nnames = NULL, ... )
  UseMethod("graphmod")

graphmod.predlogilasso <- function(obj, grenze=0, nnames = NULL,... )
  {
    plot.predlogilasso(obj,grenze,nnames,...)
  }

graphmod.cvlogilasso <- function(obj, grenze=0, nnames = NULL, lambda = NULL, ... )
  {
    obj   <- predict(obj,lambda)
    graphmod.predlogilasso(obj,grenze,nnames,...)
  }
graphmod.levellogilasso <- function(obj, grenze=0, nnames = NULL, lambda = NULL, to.which.int = NULL, ... )
  {
    obj   <- predict(obj,lambda,to.which.int)
    graphmod.predlogilasso(obj,grenze,nnames,...)
  }










