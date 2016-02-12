#' Fit zero-inflated beta regression with random effects
#'
#' @param logistic.cov the covariates in logistic component
#' @param beta.cov the covariates in beta component
#' @param Y the response variable in the regression model
#' @param subject.ind the variable for subjects
#' @param time.ind the variable for time points
#' @param quad.n Gaussian quadrature points
#' @param verbose print the fitting process
#' @return abu,theta,phi,converge
#' @export
#' @examples
#' \dontrun{
#' mssq(logistic.cov=logistic.cov,beta.cov=beta.cov,Y=Y,subject.ind=subject.ind,time.ind=time.ind)
#' }


mssq <- function(X,tc,l,N,S,
                                           max.diff=0.001,
                                           max.EM.ita = 2000,
                                           sum.by.markers=NA,
                                           species.ind=NA,
                                           species.names=NA,
                                           estimate.phi=TRUE,
                                           filter.empty.marker = FALSE,
                                           verbose=FALSE){
  ## matrix for summation by markers
  ## 2 markers per species, 5 species
  ## ncol: nspecies
  ## nrow: total nmarkers
  #1 0 0 0 0
  #1 0 0 0 0
  #0 1 0 0 0
  #0 1 0 0 0
  #0 0 1 0 0
  #0 0 1 0 0
  #...
  
  ## Empty marker: the counts across all samples are zeros,colSums(X)=0
  ## MetaPhlAn: still use that marker in calculation (use l_jk)
  ## Our model: remove that marker since phi_jk=0 and l_jk is not used in calculation
  if (filter.empty.marker){
    ### remove empty markers
    mind <- which(colSums(X,na.rm=TRUE)!=0)
    X <- X[,mind,drop=FALSE]
    l <- l[colnames(X),,drop=FALSE]
    species.names <- species.names[mind]
  }
  
  if(all(!is.na(species.names))){
    species.ind = rep(NA,length(species.names))
    k = 1
    for (spe in unique(species.names)){
      species.ind[which(species.names == spe)] <- k
      k = k + 1
    }
    sum.by.markers <- matrix(0,nrow=ncol(X),ncol=S)
    for (i in 1:S){
      sum.by.markers[species.ind==i,i] <- 1
    }
  }
  
  ## replace NA in X with 0
  na.ind <- is.na(X)
  X[na.ind] <- 0
  
  ## initial values
  l  <- as.matrix(l)
  tc <- as.matrix(tc)
  phi <- as.matrix(rep(1,ncol(X)))
  theta <- matrix(1,ncol=S,nrow=N)
  theta <- sweep(theta, 1, rowSums(theta), FUN="/")
  ita.EM <- 0
  diff <- 100
  converge <- TRUE
  
  while(ita.EM < max.EM.ita & diff > max.diff){
    ita.EM <- ita.EM + 1
    ### Save the value before each iteration
    theta.temp  <- theta
    phi.temp <- phi
    
    ### !!!!!!!!!!!!!!!!!!!!!!!!
    ### MUST ESTIMATE THETA FIRST
    ### @@@@@@@@@@@@@@@@@@@@@@@@
    
    ### estimate theta
    mat.temp <- tc %*% t(phi*l)
    mat.temp[na.ind] <- 0
    theta <- (X %*% sum.by.markers) / ( mat.temp %*% sum.by.markers)
    ### ( mat.temp %*% sum.by.markers) can be zero
    theta[is.na(theta)] <- 0
    
    ### normalize theta
    theta <- sweep(theta, 1, rowSums(theta,na.rm=TRUE), FUN="/")
    
    
    ### estimate phi
    if (estimate.phi){
      tc.aug <- tc[,rep(1,ncol(X))]
      tc.aug[na.ind] <- NA
      theta.aug <- theta[,species.ind,drop=FALSE]
      phi <- colSums(X,na.rm=TRUE) / (l*colSums(tc.aug*theta.aug,na.rm=TRUE))
      ### (l*colSums(tc.aug*theta.aug,na.rm=TRUE)) can be zero
      ### the column of tc.aug can be all NAs, marked as outliers
      ### then colSums(tc.aug*theta.aug) can be zero
      phi[is.na(phi)] <- 0
    }
    
     
    
    ### use relative difference
    ### when phi.temp=0, set it to NA 
    ### so that max will not include it
    phi.temp[phi.temp==0] <- NA
    theta.temp[theta.temp==0] <- NA
    diff <- max(c(abs((phi-phi.temp)/phi.temp),
                  abs((theta-theta.temp)/theta.temp)),na.rm=TRUE)
    if (verbose){
      cat(ita.EM,'-> max relative difference: ', diff, '\n')
      if(FALSE){
        sum(is.na((X %*% sum.by.markers)))
        sum(is.na(theta))
        sum(is.na(tc))
        colSums(theta)
        colSums(phi)
        colSums(tc.aug)
        (l*colSums(tc.aug*theta.aug,na.rm=TRUE))[which((l*colSums(tc.aug*theta.aug,na.rm=TRUE))==0),]
        X[,which((l*colSums(tc.aug*theta.aug,na.rm=TRUE))==0)]
        l[which((l*colSums(tc.aug*theta.aug,na.rm=TRUE))==0),]
      }
    }
  }
  ## convergence:  An integer code. 0 indicates successful convergence
  if (ita.EM > max.EM.ita-1 ){converge <- FALSE;cat('Not converge ',diff,'\n')}
  
  ## 
  rownames(phi) <- colnames(X)
  colnames(theta) <- unique(species.names)
  
  return(list(abu=theta*100,theta=theta,phi=phi,converge=converge))
} 