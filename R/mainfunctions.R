logilasso <- function(Y,
                      combX = NULL ,
                      to.which.int=NULL, 
                      epsilon = 0.1 ,
                      lambdainit = 1, lambdamin = 0.1 , 
                      trace = 1 ,
                      method = "groupl1"  ,
                      cvfold=1,
                      Newton = TRUE,
                      stopkrit = 1e-8, sse = NULL,X=NULL)
  {
    ## Purpose:  Returns a list with components beta, lambda and path.
    ##           Computes the whole solution path from lambda to lambdamin
    ##           If one wants the result for a single lambda, predict.log
    ##           ilasso has to be computed with the lambda.
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## Y: A contingency table. Or a vector corresponding to combX.
    ## combX: Is the combination matrix, where combX and y correspond.
    ##        This is in the rows the different combinations of exon/intron.
    ## There are 2 ways to specify the input.
    ## 1. Y in a table. 
    ## 2. Y together with combination matrix. 
    ## 3. Y together with the design matrix X. Only for experienced users
    ## to.which.int: The maximal number of interactions actually considered
    ##               in the model.
    ## epsilon: Step length.
    ## lambdainit: A starting value for lambda. Has to be specified for l2!
    ## lambdamin: Minimal lambda for which the solution path is computed
    ## trace: 0/1/2/3/4/5 with 0 = nothing,
    ##        1 = Additionally cvfold is written out
    ##        2 = additionally points
    ##        3 = Every 10th step writes the active set
    ##        4 = Additionally the new active set are written out
    ##        5 = Additionally the inner loop is documented.
    ## method: Is either groupl1, l1 or l2, depending on the penalization of 
    ##         the coefficients.
    ## cvfold: If cvfold is bigger than 1, cv is performed.
    ## Newton: If Newton Steps should be performed at all.
    ## stopkrit: The stoping criteria.
    ## sse: Random Seed, mainly used for levelcv. Do not change!! Don't give
    ##      a single number!!
    ## ----------------------------------------------------------------------
    ## Author: Corinne Dahinden, Date: 19 Jan 2006, 15:07


    if(!is.null(sse))
      .Random.seed <- sse
 
    if(!method%in%c("l1","groupl1","l2"))
      stop("Wrong method argument, currently only l1, groupl1 and l2 are implemented.")
    
    ## This is standard and only for me as user should X be given!!
    if(is.null(X))
      {
        ##If y is given as table.
        if(is.table(Y))
          {
            if(!is.null(combX))
              stop("Either specify Y as a table or together with a combination matrix.")
            nlev   <- dim(Y)
            Y      <- c(Y)
            X      <- round(f.createX(nlev,to.which.int)$mm,15)
          }
        
        ##combX must be a matrix or NULL
        if(!is.null(combX))
          {
            if(!is.matrix(combX)) #combX is not a matrix
              stop("combX argument wrong!")
            
            if(length(Y)!=nrow(combX))
              stop("combX and Y do not correspond.")
            
            comb  <- f.xy(Y,combX,to.which.int)
            Y     <- comb$y
            nlev  <- comb$nlev
            X     <- comb$modelmatrix
          }
  
        nrfac <- length(nlev)
      }

    else
      {
        
        if(!is.matrix(X))
          stop("X must be a matrix.")
        if(!is.null(combX))
          stop("Either specify X or combX, but not both!")
    
        if(is.table(Y))
          {
            warning("Y is coerced to a vector. Do you really want to specify Y together with the design matrix X?")
            Y <- c(Y)
          }
        if(nrow(X)!=length(Y))
          stop("X and Y do not correspond.")
        nrfac <- NULL
      }

    tX <- t(X)
    ##Is an index which columns belong together
    tester    <- attr(X,"assign")

    if(method == "l1" || is.null(tester))
      tester  <- 1:ncol(X)
    
    Yfinal    <- Y
    n         <- sum(Y)
    
    if(cvfold>1)
      {
        Ycv  <- rep(1:length(Y),Y)
        splitted <- split(sample(1:n),rep(1:cvfold,length=n))
        
        ##For CV we have to initialize these values:
        loss        <- NULL
        lambdafinal <- NULL
      }
    
    ## For l2, this remains for all times: these 3 conditions hold!!
    active      <- 1:ncol(X)
    anzahltmp   <- table(tester[active])
    ## diffenez is length of diffmenge
    diff0       <- TRUE
    kritfall    <- FALSE
    
### START ###    
    for(cf in 1:cvfold)
      {
        if(cvfold>1)
          {
            if(trace>0)
              cat("Run ", cf, " of ",cvfold,"\n")
            good    <- Ycv[-splitted[[cf]]]
            out     <- Ycv[splitted[[cf]]]
            
            Y    <- rep(0,length(Y))
            Yout <- rep(0,length(Y))
            Y[as.integer(names(table(good)))]    <- table(good)
            Yout[as.integer(names(table(out)))]  <- table(out)
          }
        Y         <- Y/sum(Y)
        n <- sum(Y)
        ##Computation of some needed values: Xbeta, tXY
        tXY       <- crossprod(X,Y)
        ##Would not need to compute, but is computed for the maximal lambda.
        beta      <- c(-log(nrow(X)),rep(0,ncol(X)-1))
        Xbeta     <- rep(beta[1],nrow(X))
        expXbeta  <- exp(Xbeta)
        
        nablaneu  <-  f.nabla(tXY=tXY,X=X,expXbeta=expXbeta,n=n)

        if(method!="l2"&&is.null(lambdainit))
          lambda <- max(abs(nablaneu))

        else
          lambda <- lambdainit
        
        ##If a start value for lambda is given, the corresponding beta has to be
        ##computed. This is just a first approximation!!!
        if(method == "l2"| (lambdainit!=1))
          {
            if(method!="l2")
              {
                active <- unique(c(1,which(abs(nablaneu)>lambda)))
                active <- which(tester%in%tester[active])
                anzahltmp <- c(table(tester[active]))
              }
            opt   <- optim(beta[active], f.fn, f.gradient2,
                           method = "CG",
                           n=n, tXactive = tX[active,,drop=FALSE],
                           Xactive = X[,active,drop=FALSE], tXY = tXY,
                           lambda = lambda , anzahl=anzahltmp, tester=tester[active],
                           methode=method,active=active)#,control=list(reltot=stopkrit))$par
            
            beta[active]   <- opt$par
            fval           <- opt$value
            Xbeta          <- crossprod(tX[active,,drop=FALSE],beta[active])
            expXbeta       <- exp(Xbeta)
          }
        ##If no start value is given, we just start at lambda_init
        else
          {
            active    <- 1
            lambda    <- max(abs(nablaneu))
            anzahltmp <- 1
            fvalopen     <- f.to.minimize2(Xbeta=Xbeta,Y=Y,n=n,expXbeta)
            if(method!="l2")
              fval <- fvalopen+lambda*f.group(beta[active],anzahltmp, tester[active])
            else
              fval <- fvalopen+lambda*sum(beta[-1]^2)
            
          }
        
        
###lambdatot, betatot and path are given out in the end
###lambdatot are all lambdas
###betatot is a matrix with columns beta
###path is a 2 x length(lambda) matrix with first row = active set
###and second row are the corresponding lambdas


### if lambda=lambdamin then nothing is done and only the optim function is used.
        nrlambda  <- floor((lambda-lambdamin)/(epsilon-stopkrit))+1

        lambdatot <- vector("numeric",length=nrlambda)
        betatot   <- matrix(ncol=nrlambda,nrow=ncol(X))
        lambdatot[1]<- lambda
        betatot[,1] <- beta
        path        <- rbind(active,rep(lambda,length(active)))
        
###########################################################################
#                                                                         #
#             Actual Optimization Steps                                   #
#                                                                         #
###########################################################################

        if(nrlambda>=2) for(counter in 1:(nrlambda-1))
          {
            fvalold             <- fval
            
            ##Just for the trace!
            ##For 1: . and for all 10s print lambda
            ##For 2: . active set, lambda for all 10s
            if(trace>1)
              {
                cat(".")
                if(trace>2)
                  {
                    if(method != "l2" && ((counter%%10)==1))
                      cat("\n","The active set is ",active,
                          " lambda is ",lambda,"\n")
                  }
              }
### Start ####
            lambda       <- lambda-epsilon
            
            ## New active set, because lambda was changed ##
            ## Bestimmt ob es ein kritischer Fall ist.
            nablaneu  <- f.nabla(tXY=tXY,X=X,expXbeta=expXbeta,n=n)               
            if(method!="l2")
              {
                ##New = Old active: This is good to calculate the path afterwards,
                ##see if it has changed:
                activeold <- active
                ##Maybe new active set!
                activeneu <- c(which(abs(nablaneu)> lambda))
                diffmenge <- setdiff(activeneu,active)
                diff0     <- length(diffmenge)==0
                ##Final active set
                if(!diff0)
                  {
                    active     <- sort(c(active,diffmenge))
                    active     <- which(tester%in%tester[active])
                    anzahltmp  <- c(table(tester[active]))
                    kritfall   <- any(table(tester[tester%in%tester[diffmenge]])>1)
                  }
                else
                  kritfall <- FALSE
                
                delta      <- sqrt(stopkrit)*(1+sqrt(sum(beta[active][-1]^2)))/length(active)
              }
            
            ##diff0 is automatically TRUE for the "l2" loss.
            ##Group lasso ist ausgenommen, falls diff0 nicht TRUE.
            if( !kritfall ) #Kein neues dazu, oder nur solche, die nicht stören
              {
                if(!diff0) # gehe ein wenig in richtige Richtung
                  {
                    if(trace>3)
                      cat("New active Set, optimization for lambda=",lambda,"New active set: ",
                          diffmenge,"\n")
                    
                    beta[diffmenge]<- sign(-nablaneu[diffmenge])*delta
                    Xbeta          <- Xbeta+crossprod(tX[diffmenge,,drop=FALSE],beta[diffmenge])
                    expXbeta       <- exp(Xbeta)
                    if(Newton)
                      nablaneu[active]<- f.nabla(tXY=tXY[active],X=X[,active,drop=FALSE],expXbeta=expXbeta,n=n)
                   
                 }
               if(Newton)
                 {
                   nabla2         <- f.nabla2(Xactive= X[,active,drop=FALSE],
                                              expXbeta=expXbeta , n=n)
                   nabla2active   <- f.nabla2active(nabla2,beta[active], lambda,
                                                    anzahlactive = anzahltmp,
                                                    tester = tester[active],method)
                   beta[active] <- f.newtonstep(betactive = beta[active],
                                                nabla2active = nabla2active,
                                                nablactive = nablaneu[active] ,
                                                lambda = lambda,
                                                anzahlactive=anzahltmp,
                                                tester = tester[active],
                                                method = method)
                   Xbeta        <- crossprod(tX[active,,drop=FALSE],beta[active])
                   expXbeta     <- exp(Xbeta)
                   fvalopen     <- f.to.minimize2(Xbeta=Xbeta,Y=Y,n=n,expXbeta)
                   if(method!="l2")
                     fval <- fvalopen+lambda*f.group(beta[active],anzahltmp, tester[active])
                   else
                     fval <- fvalopen+lambda*sum(beta[-1]^2)
                 }
             }
           
           if(kritfall || !Newton)
             {
               if(trace>2)
                 {
                   if(!Newton)
                     cat("No Newton, optim instead")
                   else
                     cat("New active Set, optimization with optim: ",lambda," New active set: ",                       
                     diffmenge,"\n")
                 }
               ##l2 only in this loop if !Newton, then above if is conducted.
               opt   <- optim(beta[active], f.fn, f.gradient2,
                              method = "CG",
                              n=n, tXactive = tX[active,,drop=FALSE],
                              Xactive = X[,active,drop=FALSE], tXY = tXY,
                              lambda = lambda , anzahl=anzahltmp, tester=tester[active],
                              methode=method,active=active)
               
               beta[active] <- opt$par
               fval <- opt$value
               Xbeta          <- crossprod(tX[active,,drop=FALSE],beta[active])
               expXbeta       <- exp(Xbeta)
             }
           
           if(method != "l2")
             {
               ##Set small betas equal to 0 and update active
               hvindex              <- which(abs(beta[active]) < delta)
               if(length(hvindex)>0)
                 {
                   beta[active][hvindex]<- 0
                   active <- active[-hvindex]
                   active <- unique(c(1,which(tester%in%tester[active])))
                   Xbeta        <- crossprod(tX[active,,drop=FALSE],beta[active])
                   expXbeta     <- exp(Xbeta)
                   anzahltmp    <- c(table(tester[active]))
                   
                   fvalopen         <- f.to.minimize2(Xbeta=Xbeta,Y=Y,n=n,expXbeta)
                   fval <- fvalopen+lambda*f.group(beta[active],anzahltmp, tester[active])
                   
                 }
             }
           
            
###############################################################################
#                                                                             #
#    Up to here it would be the normal optimization algorithm in the paper!!  #
#                                                                             #
###############################################################################
        
#################################################################################
#                                                                               #
#      Newtonsteps or optim if no stopcriteria has been reached                 #
#                                                                               #
#################################################################################

####Indicator if the optim function or a Newton step is performed
###First step is a Newton step. If Newton = FALSE, no NS!!
           useoptim <- FALSE
####Counter how many times no improvement is achieved.
           noimprove <- 0
           do.again <- TRUE

           starten <- sqrt(sum((betatot[,counter][active]-beta[active])^2))/(1+sqrt(sum(beta[active])^2))>sqrt(stopkrit)||(fvalold-fval)/(1+fval)>stopkrit

           while(starten && (do.again
                 ||
                 (sqrt(sum((betaold[active]-beta[active])^2))/(1+sqrt(sum(beta[active])^2))>sqrt(stopkrit)
                  ||(fvalold-fval)/(1+fval)>stopkrit)))
             {
               
               betaold       <- beta
               fvalold       <- fval
               Xbetaold      <- Xbeta
               expXbetaold   <- expXbeta
               activeoldtmp  <- active
               anzahltmpold  <- anzahltmp
               if( useoptim || !Newton) # use no Newton or useoptim
                 {
                   beta[active]   <- optim(beta[active], f.fn, f.gradient2,
                                           method = "CG",
                                           n=n, tXactive = tX[active,,drop=FALSE],
                                           Xactive = X[, active, drop=FALSE], tXY = tXY,
                                           lambda = lambda , anzahl=anzahltmp, tester = tester[active],
                                           methode=method,active=active,control=list(reltot=stopkrit))$par   
                   Xbeta        <- crossprod(tX[active,,drop=FALSE],beta[active])
                   expXbeta     <- exp(Xbeta)
                   
                 }
               else
                 {
                   nablaneu[active] <- f.nabla(tXY=tXY[active],X=X[,active,drop=FALSE],expXbeta=expXbeta,n=n)
                   nabla2       <- f.nabla2(Xactive= X[,active,drop=FALSE],
                                            expXbeta=expXbeta , n=n)
                   nabla2active <- f.nabla2active(nabla2,beta[active] ,
                                                  lambda ,
                                                  anzahlactive = anzahltmp,
                                                  tester = tester[active],
                                                  method)
                   beta[active] <- f.newtonstep(betactive =beta[active],
                                                nabla2active =nabla2active,
                                                nablactive=nablaneu[active],
                                                lambda = lambda,
                                                anzahlactive = anzahltmp,
                                                tester = tester[active],
                                                method = method)
                  
                   Xbeta    <- crossprod(tX[active,,drop=FALSE],beta[active])
                   expXbeta <- exp(Xbeta)
                   if(method=="l2")
                     fval <-  f.to.minimize2(Xbeta=Xbeta,Y=Y,n=n,expXbeta)+lambda*sum(beta[-1]^2)
                 }
               
               
               if(method!="l2")
                 {
                   hvindex              <- which(abs(beta[active]) < delta)
                   if(length(hvindex)>0)
                     {
                       beta[active][hvindex]<- 0
                       Xbeta        <- crossprod(tX[active,,drop=FALSE],beta[active])
                       expXbeta     <- exp(Xbeta)
                       active<- active[-hvindex]
                       active<- unique(c(1,which(tester%in%tester[active])))
                       anzahltmp <- c(table(tester[active]))
                     }
                   fval <-  f.to.minimize2(Xbeta=Xbeta,Y=Y,n=n,expXbeta)+lambda*f.group(beta[active],anzahltmp,tester[active])
                 }

               if(trace>4)
                 {
                   indnewton <- (!useoptim && Newton)
                   cat("Old: ",fvalold," New:",fval," Last optimization step was ", c("Newton","optim")[c(indnewton,!indnewton)]," in inner loop.", "\n")
                 }
               
               ##Don't do a step into the wrong direction
               if(fvalold<=fval)
                 {
                   useoptim<- !useoptim
                   noimprove <- noimprove+1
                   if(fvalold<fval )
                     {
                       do.again <- TRUE
                       beta      <- betaold
                       fval      <- fvalold
                       Xbeta     <- Xbetaold
                       expXbeta  <- expXbetaold
                       active    <- activeoldtmp
                       anzahltmp <- anzahltmpold
                     }
                 }
               else  {noimprove <- 0
                      do.again <- FALSE}
               
               if(noimprove>=2)
                 {
                   do.again <- FALSE
                   if(trace>1)
                     cat("|")
                 }
             }
           
           
           if(cvfold==1)
             {
               if( method %in% c("l1","groupl1"))
                 {
                   if(!setequal(active,activeold))
                     {              
                       new.minus.old <- setdiff(active,activeold)
                       old.minus.new <- setdiff(activeold,active)
                       
                                        #If new.minus.old not equal to old.minus.new
                       path <- cbind( path, rbind(new.minus.old,
                                                  rep(lambda,length(new.minus.old))),
                                     -rbind(old.minus.new,
                                            rep(lambda,length(old.minus.new))))
                     }
                 }
             }
            betatot[,counter+1]   <- beta
            lambdatot[counter+1]  <- lambda
         }
       
       
       if(cvfold>1)
         {
           loss.temp <- -crossprod(Yout,X)%*%betatot/sum(Yout)
           loss      <- f.rbindspez(loss,loss.temp)
           if(ncol(loss)==length(loss.temp))
             lambdafinal <- lambdatot
         }
       
     }
   if(is.null(to.which.int)) to.which.int <- nrfac
   
   if(cvfold==1)
     {
       if(is.vector(betatot))
         betatot <- matrix(betatot,ncol=1)
       
       to.return <- list(loss = NULL, path = path , betapath = betatot ,
                         lambdapath = lambdatot , X = X , combX = combX,
                         Y = Yfinal , cvfold= cvfold, trace = trace,
                         corresp = tester , nrfac = nrfac , 
                         to.which.int = to.which.int ,  method = method,
                         stopkrit=stopkrit,
                         epsilon = epsilon,
                         Newton = Newton)
       class(to.return) <- "logilasso"
     }
   else
     {
       to.return <- list(loss=loss, path = NULL, betapath = NULL,
                         lambdapath=lambdafinal , X=X , combX = combX,
                         Y=Yfinal , cvfold= cvfold, trace = trace,
                         corresp = tester , nrfac = nrfac , 
                         to.which.int = to.which.int ,  method = method,
                         stopkrit= stopkrit,
                         epsilon = epsilon,
                         Newton = Newton)
       class(to.return) <- c("cvlogilasso","logilasso")
     }
   to.return
 }

levelcv     <-   function(Y ,combX = NULL ,
                           from.which.int = 1, to.which.int=NULL, 
                           epsilon = 0.1 ,
                           lambdainit =1, lambdamin = 0.1 , 
                           trace = 1 ,
                           method = "groupl1"  ,
                           cvfold=10,
                           Newton = TRUE,
                           stopkrit = 1e-8)

  ## Purpose: Computes for all levels of interaction the logilasso.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date: 20 Sep 2006, 19:14

  {
    rnorm(1)
    saveseed <- .Random.seed
    ret  <- list()
    ind  <- to.which.int-from.which.int+1
    
    ret[[ind]] <- logilasso(Y,combX=combX,
                            to.which.int=to.which.int,
                            epsilon=epsilon,
                            lambdainit=lambdainit,lambdamin=lambdamin,
                            trace=trace, method=method,
                            cvfold=cvfold, Newton = Newton,
                            stopkrit = stopkrit,sse=saveseed)
    
    nrfac  <- ret[[ind]]$nrfac
    Y      <- ret[[ind]]$Y
    X      <- ret[[ind]]$X

    tester <- attr(X,"assign")
   

    if(is.null(to.which.int))
      to.which.int <- nrfac

    if(to.which.int>from.which.int)
      {
        for(i in (from.which.int):(to.which.int-1))
          {
            Xt <- X[,c(1:sum(choose(nrfac,0:i)))]
            attributes(Xt)$assign <- attr(X,"assign")[c(1:sum(choose(nrfac,0:i)))]
            ret[[i-from.which.int+1]] <- logilasso(Y,to.which.int=i,
                                                   epsilon=epsilon,lambdamin=lambdamin,trace=trace,
                                                   lambdainit=lambdainit,
                                                   method=method,
                                                   cvfold=cvfold,Newton = Newton,
                                                   stopkrit = stopkrit,sse=saveseed,X=Xt)
            ret[[i-from.which.int+1]]$nrfac <- nrfac
          }
      }

    class(ret) <- c("levellogilasso","cvlogilasso","logilasso")
    ret
  }





