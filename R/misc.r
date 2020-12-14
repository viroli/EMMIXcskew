#
#  EM algorithm for Mixture of Multivariate Canonical Fundamental Skew t-distributioins
#  Package: EMMIX-cskew
#  Version: 0.9-1
#
#  Code by S. X. Lee
#  Updated on 07 Sep 2015
#
# Lee S.X. and Mclachlan, G.J. (2015) Finite mixtures of canonical fundamental 
# skew t-distributions: the unification of the restricted and unrestricted 
# skew t-mixture models. Statistics and Computing. doi:10.1007/s11222-015-9545-x
#
# The following code are adapted from other packages
#
################################################################################
#  misc.r
#                           Miscellaneous Tools
#
################################################################################

#skewness
skewness <- function (x, na.rm = FALSE){
    if (is.matrix(x))
        apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
    if (na.rm) x <- x[!is.na(x)]
	       n <- length(x)
        (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x))
        sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
}

#mahalanobis
mahalanobis. <- function(x, center, cov, inverted=FALSE, ...) {
    x <- if(is.vector(x)) matrix(x, ncol=length(x)) else as.matrix(x)
    x <- t(sweep(x, 2, center))# = (x - center)
    retval <- colSums(x * if(inverted) cov%*%x else solve(cov, x, ...))
    names(retval) <- rownames(x)
    retval
}

mahalanobis <- function(x, center, cov, inverted=FALSE, ...){
    x <- if(is.vector(x)) matrix(x, ncol=length(x)) else as.matrix(x)
    x <- sweep(x, 2, center)# = (x - center)
    if(!inverted)
	cov <- solve(cov, ...)
    retval <- rowSums((x%*%cov) * x)
    names(retval) <- rownames(x)
    retval
}

#
isPosDef <- function(M) { 
    if ( all(M == t(M) ) ) { 
        if (  all(eigen(M)$values > 0) ) {TRUE}
        else {FALSE} 
    }else {FALSE}  
} 

#
is.whole <- function(a) { 
   (is.numeric(a) && floor(a)==a) ||
   (is.complex(a) && floor(Re(a)) == Re(a) && floor(Im(a)) == Im(a))
}


#multivariate normal density
dmn <- function(x, mu=rep(0,d), Sigma, log=FALSE) {
      d  <- if(is.matrix(Sigma)) ncol(Sigma) else 1
      if(d>1 & is.vector(x)) x <- matrix(x, 1, d)
      n  <- if(d==1)  length(x) else nrow(x) 
      X  <- t(matrix(x, nrow=n, ncol=d)) - mu
      Q  <- apply((solve(Sigma)%*% X)* X, 2, sum) 
      logDet <- sum(logb(abs(diag(qr(Sigma)[[1]]))))
      logPDF <- as.vector(Q + d*logb(2*pi)+logDet)/(-2)
      if(log) logPDF else exp(logPDF)
}

#multivariate t density
dmt <- function(dat, mu, sigma, dof = Inf, log = FALSE) {
    if (dof == Inf)  return(dmn(dat, mu, sigma, log = log))
    dof <- round(dof)
    d  <- if(is.matrix(sigma)) ncol(sigma) else 1
    n  <- if(d==1) length(dat) else nrow(dat)
    X <- t(matrix(dat, nrow = n, ncol = d)) - mu
    Q <- apply((solve(sigma) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(sigma)$qr))))
    logPDF <- (lgamma((dof + d)/2) - 0.5 * (d * logb(pi * dof) + logDet)
               - lgamma(dof/2) - 0.5 * (dof + d) * logb(1 + Q/dof))
    if(log) logPDF else exp(logPDF)
}


#error rate
error.rate <- function(clust1,clust2){
    if((n=length(clust1))!=length(clust2))  stop("error: length not equal")
    if( (g=length(table(clust1)))!=length(table(clust2))) stop("the number of clusters are not equal")
    permute<-function(a){
        n<-length(a)
        if(n==1) f<-a
        else{
            nm<-gamma(n)
            f<-array(0,c(n,n*nm))
            j<-1
        
            for(i in a){
                f[1, (j-1)*nm+1:nm]<-i
                f[-1,(j-1)*nm+1:nm]<-permute(setdiff(a,i))
                j<-j+1
            }
        }    
        f
    }
       
    id<-1:n
    cmb<-permute(1:g)
    nperm<-ncol(cmb)
    rate<-rep(0,nperm)
    
    for(i in 1:nperm){
        tmp<-rep(0,g)
        tc<-rep(0,n)
        for(j in 1:g)
            tc[clust2==j]=cmb[j,i]
    
        for(j in 1:g){  
            tmp1<-0 
            for(k in (1:g)[-j])
                tmp1<-tmp1+length(intersect(id[clust1==j],id[tc==k]))
            tmp[j]<-tmp1
        }
        rate[i]<-sum(tmp)/n
    }
    min(rate)
}


