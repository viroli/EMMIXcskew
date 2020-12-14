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

################################################################################
#  density.r
#             density functions
#
################################################################################

dcfust <- function(dat, mu=NULL, sigma=NULL, delta=NULL, dof=1, known=NULL) { 
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
        dof <- known$dof
    } 
    if(missing(dat)) stop("dat must be provided")     
    if(is.null(delta)) stop("delta is missing") else q <- ncol(delta)
    n <- nrow(dat);  p <- ncol(dat);           
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)
    if(is.null(mu)) mu <- rep(0,p)
    if(is.null(sigma)) sigma<-diag(p)
    mu <- as.numeric(mu)
    Ydiff <- t(dat) - mu%*% matrix(1,1,n)               
    Omega <- sigma + delta %*% t(delta)                 
    InvOmega <- solve(Omega)                           
    Lambda <- diag(q) - t(delta) %*% InvOmega %*% delta 
    Q <- t(delta) %*% InvOmega %*% Ydiff               
    eta <- mahalanobis(dat, as.numeric(mu), InvOmega, T)
    N <- sqrt((dof+p)/abs(dof+eta))                     
    T1 <- dmt(dat, as.numeric(mu), Omega, dof);  T2 <- T1 
    for (i in 1:n) T2[i] <- pmt(N[i]*Q[,i], rep(0,q), Lambda, dof+p)
    ST <- 2^q * T1 * T2
    return(ST)
}


dfmcfust <- function(dat, mu=NULL, sigma=NULL, delta=NULL, dof=NULL, pro=NULL, known=NULL) {
    n <- nrow(dat);  p <- ncol(dat); g <- length(pro)
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)  
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
        dof <- known$dof
        pro <- known$pro
    }
    if(missing(dof)) stop("dof is missing.")
    if(g==1) return(dcfust(dat, known=as.fmcfust(1, mu, sigma, delta, dof, pro))) 
    if(is.null(dof)) dof <- rep(1,g)
    if(is.null(pro)) pro <- rep(1/g,g)
    ST <- rep(0, length=n)
    for (i in 1:g) ST <- ST + pro[i]*dcfust(dat, as.numeric(mu[[i]]), sigma[[i]], delta[[i]], dof[i])
    return(ST)
}

dmsn <- function(dat, mu=NULL, sigma=NULL, delta=NULL, known=NULL) {
    if(missing(dat)) stop("dat must be provided")
    n <- nrow(dat);  p <- ncol(dat);                
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
    } 
    if(is.null(mu)) mu <- rep(0,p)
    if(is.null(delta)) delta <- rep(0,p)
    if(is.null(sigma)) sigma<-diag(p)
    mu <- as.numeric(mu)
    delta <- as.numeric(delta)
    Delta <- diag(p); diag(Delta) <- delta            
    Ydiff <- t(dat) - mu%*% matrix(1,1,n)             
    Omega <- sigma + Delta %*% t(Delta)               
    InvOmega <- solve(Omega)                          
    Lambda <- diag(p) - Delta %*% InvOmega %*% Delta  
    Q <- Delta %*% InvOmega %*% Ydiff                 
    N1 <- dmn(dat, mu=as.numeric(mu), Sigma=Omega, log=FALSE);  N2 <- N1     
    for (i in 1:n) {
        N2[i] <- sadmvn(lower=rep(-Inf,p), upper=Q[,i], mean=matrix(rep(0,p),p), varcov=Lambda) 
    }
    SN <- 2^p * N1 * N2
    return(SN)
}


dmst <- function(dat, mu=NULL, sigma=NULL, delta=NULL, dof=1, known=NULL) {
    if(dof==Inf) return(dmsn(dat, mu, sigma, delta, known))
    if(missing(dat)) stop("dat must be provided")
    n <- nrow(dat);  p <- ncol(dat);                   
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
        dof <- known$dof
    } 
    if(is.null(mu)) mu <- rep(0,p)
    if(is.null(delta)) delta <- rep(0,p)
    if(is.null(sigma)) sigma<-diag(p)
    mu <- as.numeric(mu)
    delta <- as.numeric(delta)
    Delta <- diag(p); diag(Delta) <- delta            
    Ydiff <- t(dat) - mu%*% matrix(1,1,n)             
    Omega <- sigma + Delta %*% t(Delta)               
    InvOmega <- solve(Omega)                          
    Lambda <- diag(p) - Delta %*% InvOmega %*% Delta  
    Q <- Delta %*% InvOmega %*% Ydiff                 
    eta <- mahalanobis(dat, as.numeric(mu), InvOmega, T)  
    N <- sqrt((dof+p)/abs(dof+eta))                   
    T1 <- dmt(dat, as.numeric(mu), Omega, dof);  T2 <- T1     
    for (i in 1:n) T2[i] <- pmt(N[i]*Q[,i], rep(0,p), Lambda, round(dof+p))
    ST <- 2^p * T1 * T2
    return(ST)
}

dfmmst <- function(dat, mu=NULL, sigma=NULL, delta=NULL, dof=NULL, pro=NULL, known=NULL) {
    n <- nrow(dat);  p <- ncol(dat)
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)  
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
        dof <- known$dof
        pro <- known$pro
    }
    if(missing(dof)) stop("dof is missing.")
    g <- length(dof)
    if(g==1) return(dmst(dat, known=as.fmmst(1, mu, sigma, delta, dof))) 
    if(is.null(dof)) dof <- rep(1,g)
    if(is.null(pro)) pro <- rep(1/g,g)
    ST <- rep(0, length=n)
    for (i in 1:g) ST <- ST + pro[i]*dmst(dat, as.numeric(mu[[i]]), sigma[[i]], as.numeric(delta[[i]]), round(dof[i]))
    return(ST)
}




