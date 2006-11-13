#### Nun: 3 Hilfsfunktionen in logilasso ########

f.createX <- function( nlev , to.which.int = NULL )
{
  ## Purpose: Computes the design matrix 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date: 14 Jun 2006, 17:04
  
  le <- length(nlev)
  
  if(is.null(to.which.int))
    to.which.int <- le

  defaultcontrast <- options()$contrasts
  options(contrasts=c("contr.poly","contr.poly"))
  on.exit(options(contrasts=defaultcontrast))
  
  cumpr <- cumprod(nlev)
  tot   <- prod(nlev)
  

  fakt <- list()
  fakt[[1]] <- gl(nlev[1],1,tot)

  if(length(nlev)>1)
    {
      for(i in 2:length(nlev))
        {
          fakt[[i]] <- gl(nlev[i],cumpr[i-1],tot)
        }
    }
  
  names(fakt) <- as.character(1:length(nlev))
  dd          <- data.frame(fakt)

  combination.matrix <- apply(dd,2,as.integer)-1

  if(to.which.int>1)
    {
      mm <- model.matrix(formula(paste("~.^",to.which.int,sep="")),data=dd)
      model <- cbind(rep(0,length(nlev)),attr(terms(formula(paste("~.^",to.which.int,sep="")),data=dd),"factors"))
    }
  if(to.which.int==1)
    {
      mm <- model.matrix(~.,data=dd)
      model <- cbind(rep(0,length(nlev)),attr(terms(~.,data=dd),"factors"))
    }
    

  tester <- attr(mm,"assign")+1
  ncc <- nrow(mm)
  
  mm <- apply(mm,2,scale,center=FALSE)*sqrt(ncc/(ncc-1))

  attributes(mm)$assign <- tester
  list(mm=mm,combination.matrix=combination.matrix,model=model)
 
}

f.rbindspez <- function(x,y)
{
  ## Purpose: Takes 2 vectors or a matrix and a vector as input of not
  ##          necessarily the same length and makes rbind resulting in
  ##          a matrix with ncol=dimension of the longer one.
  ##          The shorter is expanded by reproducing the first element
  ##          as many times as needed
  ## ----------------------------------------------------------------------
  ## Arguments: x and y
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date: 19 Jan 2006, 10:18


  #Wenn eines Null ist, gebe einfach das andere zurück
  if(is.null(x))
    return(y)
  if(is.null(y))
    return(x)
    
  if(!is.matrix(x))
    x <- matrix(x,nrow=1)
      
  if(!is.matrix(y))
    y <- matrix(y,nrow=1)
      
  n1 <- dim(x)[2]
  n2 <- dim(y)[2]
  
  if(n1>n2)
    y <- cbind(matrix(rep(y[,1],(n1-n2)),ncol=(n1-n2),byrow=F),y)
  
  if(n1<n2)
    x <- cbind(matrix(rep(x[,1],(n2-n1)),ncol=(n2-n1),byrow=F),x)
  
  toret <- rbind(x,y)
  toret
}


f.xy <- function(y,X,to.which.int=NULL)
  #Takes y and the combination matrix X and gives out the final Y and the model matrix X.
  {
    if(is.null(to.which.int))
      to.which.int <- ncol(X)

    defaultcontrast <- options()$contrasts
    options(contrasts=c("contr.poly","contr.poly"))
    on.exit(options(contrasts=defaultcontrast))
    #X must be given in levels of the form 0/1/2 and so on.
    
    nlev <- round(apply(X,2,function(x) length(unique(x))),3)
    
    cumpr <- cumprod(nlev)
    tot   <- prod(nlev)
    
    fakt <- list()
    fakt[[1]] <- gl(nlev[1],1,tot)
    
    for(i in 2:length(nlev))
      {
        fakt[[i]] <- gl(nlev[i],cumpr[i-1],tot)
      }
    names(fakt) <- as.character(1:length(nlev))
    dd          <- data.frame(fakt)

    if(to.which.int==1)
      {
        model <- cbind(rep(0,length(nlev)),attributes(terms(formula(paste("~.")),data=dd))$factors)
        mm    <- model.matrix(formula(paste("~.")),data=dd)
      }
    else
      {
        model <- cbind(rep(0,length(nlev)),attributes(terms(formula(paste("~.^",to.which.int,sep="")),data=dd))$factors)
        mm    <- model.matrix(formula(paste("~.^",to.which.int,sep="")),data=dd)
      }

    tester <- attributes(mm)$assign+1
    ncc    <- nrow(mm)
    mm     <- apply(mm,2,scale,center=FALSE)*sqrt(ncc/(ncc-1))

    attributes(mm)$assign <- tester
    
    mati <- apply(dd,2,as.integer)-1
    
    u <- apply(mati,1,paste,collapse="")
    x <- apply(X,1,paste,collapse="")

    index <- match(x,u)

    if(is.null(index))
      stop("Levels of second argument must be 0,1 and so on.")

    yneu        <- rep(0,nrow(mati))
    yneu[index] <- y
    
    list(y=yneu,nlev=nlev,modelmatrix=mm,combination.matrix=mati,model=model)
  }

f.nabla <- function(tXY,X,expXbeta,n)
  (-tXY + crossprod(X,expXbeta)*n)/n
#calculates the gradient without penalty

f.nabla2 <- function(Xactive,expXbeta,n)
  {
    #Calculates ony the second derivative of the function to minimize
    #without penalty
    crossprod((Xactive*c(expXbeta)),Xactive)*n
  }


### Function to optimize, mainly for the general case. In optim.
f.fn <- function(betactive, n, tXactive, Xactive, tXY, lambda, anzahl , tester, methode , active)
  {
    ## Purpose: Objective function. This function has to be minimized
    ##          in the step algorithm. This function calculates the
    ##          function value. Total!! Is used in optim .
    ##          Very general.
    ## ----------------------------------------------------------------------
    ## Arguments: Clear!
    ## anzahl: Welche levels zusammengehoeren werden in einer Zahl zusammen-
    ##         gefasst. Output von f.createX
    ## ----------------------------------------------------------------------
    ## Author: Corinne Dahinden, Date: 19 Jan 2006, 15:16

    if(methode == "l2")
      return((f.to.minimize(betactive,n,tXactive,tXY,active=active))/n+lambda*sum(betactive[-1]^2))
 
    else
      (f.to.minimize(betactive,n,tXactive,tXY,active=active))/n+lambda*f.group(betactive,anzahl,tester)
  }

f.gradient2 <- function(beta,n,tXactive,Xactive,tXY,lambda,anzahl,tester,
                        methode,active,...)
  {
    ##Is used in optim
    gradwopen <- -tXY[active]+crossprod(Xactive,exp(crossprod(tXactive,beta)))*n
    if(methode=="l2")
      to.add <- 2*beta[active]*lambda
    else
      {
        L      <- rep(f.group(beta,anzahl,tester,summieren=FALSE,pen0=TRUE),
                      times=anzahl)
        to.add <- lambda*beta*c(0,1/L[-1])
        ##if L is 0, then we just take the negative gradient sign. Should not happen!!
        to.add[to.add%in%c(-Inf,Inf,NaN)] <- sign(-gradwopen[to.add%in%c(-Inf,Inf,NaN)])*lambda
      }
    gradwopen+to.add
  }

########### 2 Functions to calculate the raw value without penalty #########
# First is just for the standard case. The second is general. 
#Must depend on betactive and active at least.

f.to.minimize2 <- function(Xbeta, Y, n, expXbeta)
  {#schnell!
   #If the fval has to be computed quick, this function is taken.
   #In the preopt function and in logilasso it is applied.
    -sum(Y*Xbeta) +n*sum(expXbeta)
  }

f.to.minimize <- function(betactive, n, tXactive, tXY, active)
  {#langsam
   #Is the normal f.to.minimize function. This function depends on betactive
   #active and a number of other parameters.
    -crossprod(tXY[active],betactive)+n*sum(exp(crossprod(tXactive,betactive)))
  }    

########### Is used everywhere to compute fval ############

f.group <- function(x,anzahl,tester,pen0=FALSE,summieren=TRUE)
  {
  ## Purpose: Calculates the group lasso norm.
  ##          Is used everywhere to compute the fval: f.fn (minimize),
  ##          f.nabla2active for second derivative. 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## anzahl: Die Anzahl, welche zu einem bestimmten betasubtyp gehoeren.
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date:  2 Jun 2006, 14:13

    if(all(anzahl==1))
      norm <- abs(x)
    else
    norm <- tapply(x,tester,function(x) sqrt(sum(x^2)))

    if(summieren)
      {
        if(pen0) return(sum(norm))
        else  #pen0 == FALSE
          return(sum(norm[-1]))
        
      }

    else  #summieren == FALSE
      {
        if(pen0)
          return(norm)
        else
          return(norm[-1])
      }


  }

f.newtonstep <- function(betactive, nabla2active, nablactive, lambda,
                         anzahlactive, testera, method)
  {
    ## Purpose: Performs a newtonstep and gives out the new beta.
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Corinne Dahinden, Date: 20 Jan 2006, 12:15

    if(method=="l2")
      {
        gradient     <- nablactive+2*betactive*lambda
        gradient[1]  <- 0
        schrittweite <- 1
      }

    else
      gradient <- f.gradientfull(nablactive, betactive, lambda,
                                 anzahlactive, testera, method)                                
      
    #Deltabeta, but one might not be able to do the full step,
    #but stop when beta gets zero.
    deltabeta    <- try(solve(nabla2active,-gradient))

    if(class(deltabeta)=="try-error")
      return(betactive)
    
    if(method !="l2")
      {
        #How far one can go (in a negative fraction of deltabeta)
        #in each direction. If 0, then one can go as far
        #as necessary, else the negative fraction is in x.
        x            <- betactive[(betactive/deltabeta)<0]/deltabeta[(betactive/deltabeta)<0]
        #Gives Inf, if only 0, but with the min it gives the maximal
        #step size.
        schrittweite <- min(abs(x),1)
      }
    
    e    <- betactive + schrittweite*deltabeta
    e   
  }
######## Next two functions calculate the FULL gradient and second derivative
## correspond to the functions above!! nabladefault plus penalty=f.gradientfull
## nabla2default plus second derivative of penalty. In f.newtonstep.

f.gradientfull <- function(nablactive , betactive , lambda , anzahlactive ,
                           tester, method)
{
  ## Purpose: Computes the whole gradient vector! For the function of
  ##          interest plus the penalty of interest.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date:  6 Feb 2006, 18:18

  if(method=="l2")
    to.add <- 2*lambda*c(0,betactive[-1])
  else
    {
      L      <- rep(f.group(betactive,anzahlactive,tester,summieren=FALSE,pen0=TRUE),
                    times=anzahlactive)
      to.add <- lambda*betactive*c(0,1/L[-1])
      #if L is 0, then we just take the negative gradient sign. Should not happen!!
      to.add[to.add%in%c(-Inf,Inf,NaN)] <- sign(-nablactive[to.add%in%c(-Inf,Inf,NaN)])*lambda
    }
  
  nablactive+to.add
 
}


f.nabla2active <- function( nabla2, betactive, lambda,  anzahlactive, tester, method)
{
  ## Purpose: Is full nabla2 with group penalty second derivative.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Corinne Dahinden, Date:  2 Jun 2006, 16:37


  if(method=="l2")
    return(nabla2+2*lambda*diag(c(0,rep(1,length(betactive)-1))))
  
  if(all(anzahlactive==1))  
    return(nabla2)


  L <- f.group(betactive,anzahlactive,tester,summieren=FALSE,pen0=TRUE)
  

  betactive <- matrix(betactive/rep(L,times=anzahlactive),nrow=1)
  betactive[is.na(betactive)] <- 0
  a <- crossprod(betactive,betactive)

  b      <- matrix(0,ncol=length(betactive),nrow=length(betactive))
  
  for(i in unique(tester)[-1])
    b[tester==i,tester==i] <- a[tester==i,tester==i]/L[i]
  
  b[is.na(b)] <- 0
  to.add <- diag(c(0,1/rep(L,times=anzahlactive)[-1]))
  to.add[to.add%in%c(-Inf,Inf,NaN)] <- 0
  nabla2active <- nabla2+to.add*lambda- b*lambda
  nabla2active
}



f.adplot <- function(adjac,nnames=NULL,directed=FALSE,...)
  {
    if(!is.null(nnames))
      rownames(adjac) <- colnames(adjac) <- nnames

    if(!directed)
      {
        if(any(adjac!=t(adjac)))
          {
            adjact <- adjac+t(adjac)
            index <- which(adjact>0)
            adjac <- matrix(0,ncol=ncol(adjact),nrow=nrow(adjact))
            adjac[index] <- 1
          }
        trueadjacency <- adjac
        edgemode <- "undirected"
      }
    else
      {
        trueadjacency <- adjac
        edgemode <- "directed"
      }
    truegraphnel  <- new("graphAM",trueadjacency,edgemode=edgemode)
    Nel           <- as(truegraphnel,"graphNEL")
    
    plot(Nel,...)
  }
