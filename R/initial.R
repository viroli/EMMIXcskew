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
#  initial.r
#             Initial values for fmcfust
#
################################################################################

init.cfust <- init.fmcfust <- function(g, dat, q=p, initial=NULL, known=NULL, clust=NULL, nkmeans=20, method=c("moments","transformation","EMMIXskew","EMMIXuskew")) {
    Y <- t(dat); P <- dim(Y);  k <- p <- P[1];  N <- P[2]; fixed <- known
    w <- options("warn"); on.exit(options(w)); options(warn=-1)
   
    if(method[1]=="EMMIXskew") {if(is.null(initial)) stop("initial is NULL. To use method=\"EMMIXskew\", please provide a fitted model obtained from EMMIXskew as initial.") else initial <- EMMIX2cfust(initial,q=q)}
    if(method[1]=="EMMIXuskew") {if(is.null(initial)) stop("initial is NULL. To use method=\"EMMIXuskew\", please provide a fitted model obtained from EMMIXuskew as initial.") else initial <- fmmst2cfust(initial)}  
    if(method[1]=="restricted") {if(is.null(initial)) stop("initial is NULL. To use method=\"restricted\", please provide a fitted model obtained from EMMIXcskew with q=1 as initial.") else initial <- restricted2cfust(initial,q=q)}  
    
    MU <- initial$mu
    SIGMA <- initial$sigma
    DELTA <- initial$delta
    DOF <- initial$dof
    PI <- initial$pro
    if (!is.null(fixed$mu)) MU <- fixed$mu
    if (!is.null(fixed$sigma)) SIGMA <- fixed$sigma
    if (!is.null(fixed$delta)) DELTA <- fixed$delta
    if (!is.null(fixed$dof)) DOF <- fixed$dof
    if (!is.null(fixed$pro)) PI <- fixed$pro
    MUflag <- is.null(MU); SIGMAflag <- is.null(SIGMA); DELTAflag <- is.null(DELTA); PIflag <- is.null(PI); DOFflag <- is.null(DOF)
    if(!PIflag) if(!checkPI(g, PI)) {
        cat("WARNING: in fmcfust initialisation, pro is not correctly specified.\n")
        PIflag <- TRUE
    }
    if(!DOFflag) if(!checkDOF(g, DOF)) {
        cat("WARNING: in fmcfust  initialisation, dof is not correctly specified.\n")
        DOFflag <- TRUE
    }
    if(!MUflag) if(!checkMU(g, k, MU)) {
        cat("WARNING: in fmcfust  initialisation, mu is not correctly specified.\n")
        MUflag <- TRUE
    }
    if(!SIGMAflag) if(!checkSIGMA(g, k, SIGMA)) {
        cat("WARNING: in fmcfust  initialisation, sigma is not correctly specified.\n")
        SIGMAflag <- TRUE
    }
    if(!DELTAflag) if(!checkDELTA(g, k, q, DELTA)) {
        cat("WARNING: in fmcfust  initialisation, delta is not correctly specified.\n")
        DELTAflag <- TRUE
    }
       
    if (MUflag || SIGMAflag || DELTAflag || PIflag || DOFflag) {    
        if(MUflag) MU <- list()
        if(SIGMAflag) SIGMA <- list()
        if(DELTAflag) DELTA <- list()
        if(PIflag) PI <- vector()
        if(DOFflag) DOF <- vector()    
        
        if(method[1]=="transformation") {if(is.null(clust)) clust<-kmeans(dat, g)$cluster; INIT <- fust_init_eigen(dat, g, clust, itmax=20)}
        else { 
            INIT <- fust_init(g,Y,q,initial,known,clust,nkmeans=nkmeans, k=k, N=N, MUflag=MUflag, SIGMAflag=SIGMAflag, DELTAflag=DELTAflag, DOFflag=DOFflag, PIflag=PIflag, MU=MU, SIGMA=SIGMA, DELTA=DELTA, DOF=DOF, PI=PI)               
        } 
    } else{  
        TT <- fust_computeTAU(g, t(Y), MU, SIGMA, DELTA, PI, DOF)
        logL <- TT$logL 
        TAU <- TT$TAU
        maxRESULT <- list("MU"=MU, "SIGMA"=SIGMA, "DELTA"=DELTA, "PI"=PI, "DOF"=DOF, "logL"=logL, "TAU"=TAU, "pro"=PI, "mu"=MU, "sigma"=SIGMA, "delta"=DELTA, "dof"=DOF)
        maxRESULT$fflag <- list("MU"=!is.null(fixed$mu), "SIGMA"=!is.null(fixed$sigma), "DELTA"=!is.null(fixed$delta), "PI"=!is.null(fixed$pro), "DOF"=!is.null(fixed$dof))
        INIT <- maxRESULT
    }
    INIT$loglik <- INIT$logL
    return(INIT)  

}


EMMIX2cfust <- function(Fit, i=NULL, q=p) {
    g <- length(Fit$pro)
    p <- nrow(Fit$mu)
    nFit <- Fit
    nFit$mu <- nFit$delta <- nFit$sigma <- list()
    if(all(Fit$dof==0)) Fit$dof <- rep(100,g)
    if(is.null(i) || !(i %in% 1:g)) {
        for(i in 1:g) {
            nFit$mu[[i]] <- matrix(Fit$mu[,i],p)
            nFit$delta[[i]] <- matrix(c(Fit$delta[,i],rep(0,(q-1)*p)), p, q)
            nFit$sigma[[i]] <- matrix(Fit$sigma[,,i],p,p)
        }
    } else{ 
        nFit$pro <- 1
        nFit$dof <- Fit$dof[i]
        nFit$mu[[1]] <- matrix(Fit$mu[,i],p)
        nFit$delta[[1]] <- matrix(c(Fit$delta[,i]), rep(0,(q-1)*p),p,q)
        nFit$sigma[[1]] <- matrix(Fit$sigma[,,i],p,p)
    }
    return(nFit)
}

fmmst2cfust <- function(Fit, i=NULL) {
    g <- length(Fit$pro)
    p <- nrow(Fit$mu[[1]])
    nFit <- Fit
    nFit$delta <- list()
    if(is.null(i) || !(i %in% 1:g)) {
        for(i in 1:g) {nFit$delta[[i]] <- diag(p); diag(nFit$delta[[i]])<-Fit$delta[[i]]}
    } else{
        nFit$pro <- 1
        nFit$dof <- Fit$dof[i]
        nFit$mu[[1]] <- Fit$mu[[i]]
        nFit$delta[[1]] <- diag(p); diag(nFit$delta[[1]])<-Fit$delta[[i]]
        nFit$sigma[[1]] <- Fit$sigma[[i]]
    }
    return(nFit)
}

restricted2cfust <- function(Fit, q=p, i=NULL) {
    g <- length(Fit$pro) 
    p <- nrow(Fit$mu[[1]])
    nFit <- Fit
    nFit$delta <- list()
    if(is.null(i) || !(i %in% 1:g)) {
        for(i in 1:g) {nFit$delta[[i]] <- matrix(0, p, q); nFit$delta[[i]][,1] <- Fit$delta[[i]]}
    } else{
        nFit$delta[[1]] <- matrix(0, p, q); nFit$delta[[i]][,1] <- Fit$delta[[i]]
    }
    return(nFit)
}      



fust_init_eigen <- function(Y, g, clust, itmax=20) {
    INIT <- list(); p <- ncol(Y)
    INIT$clust <- clust; INIT$pro <- as.numeric(table(clust)/nrow(Y))
    for(i in 1:g) {
        clusti <- c(1:nrow(Y))[clust==i]
        COV <- cov(Y[clusti,]) 
        CC <- eigen(COV, T)
        CCC <- CC$vectors
        Yst <- t(solve(CCC)%*%t(Y[clusti,]))
        Fiti <- fmmst(1,Yst, itmax=itmax, print=F)
        INIT$mu[[i]] <- CCC%*%Fiti$mu[[1]]
        INIT$sigma[[i]] <- CCC %*% Fiti$sigma[[1]] %*% t(CCC)
        INIT$delta[[i]] <- CCC %*% diag(as.numeric(Fiti$delta[[1]]),p)
        INIT$dof[i] <- Fiti$dof[1]
    }           
    TT <- fust_computeTAU(g, Y, INIT$mu, INIT$sigma, INIT$delta, INIT$pro, INIT$dof)
    INIT$logL <- TT$logL; INIT$TAU <- TT$TAU  
    return(INIT)
}

fust_init <- function(g, Y, q, initial=NULL, fixed=NULL, clust=NULL, ndelta=1, nkmeans=100, k=1, N=1, MUflag, SIGMAflag, DELTAflag, DOFflag, PIflag, MU, SIGMA, DELTA, DOF, PI) {           
    N <- ncol(Y); k <- p <- nrow(Y)   
    init_param <- function (Fit, g, Y, k, q, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag, ndelta=1){
        Fit$MU <- Fit$SIGMA <- Fit$DELTA <- list()
        if(is.null(Fit$size)) Fit$size <- as.numeric(table(Fit$cluster))
        for (i in 1:g) {
            selectY <- {if(k==1) t(matrix(Y[Fit$cluster==i],1)) else t(Y[,Fit$cluster==i])}            
            if(MUflag)  Fit$MU[[i]] <- {if(k==1) mean(Y[Fit$cluster==i]) else matrix(rowSums(Y[,Fit$cluster==i])/Fit$size[i],k,1) } else Fit$Mu[[i]] <- MU[[i]]
            if(SIGMAflag) Fit$SIGMA[[i]] <- cov(selectY) else Fit$SIGMA[[i]] <- SIGMA[[i]]
            skew <- skewness(selectY)
            if(DELTAflag) {  #1:q, 2:p, 3:diagonal, 4:r                                              
                if(ndelta==3 || ndelta==2) Fit$DELTA[[i]] <- diag(-5*(skew <(-0.1)) + 5*(skew>0.1))       
                else {Fit$DELTA[[i]] <- matrix(0, k, q); Fit$DELTA[[i]][,1] <- 3*sign(skew)}                      
            } else Fit$DELTA[[i]] <- DELTA[[i]]
        }
        if(DOFflag) Fit$DOF <- rep(4, g) else Fit$DOF <- DOF
        if(PIflag) Fit$PI <- Fit$size/N  else Fit$PI <- PI  
        TT <- fust_computeTAU(g, t(Y), Fit$MU, Fit$SIGMA, Fit$DELTA, Fit$PI, Fit$DOF)
        Fit$logL <- TT$logL; Fit$TAU <- TT$TAU        
        return(Fit)  
    }         
                
    if(is.null(clust)) {    
        Try <- kmeans(t(Y), g);   
        if(any(Try$size < k)) {Try$logL <- maxLL <- -Inf}    
        maxRESULT <- Try <- init_param(Try, g, Y, k, q, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag, ndelta)
        maxLL <- Try$logL
        if(g > 1 && nkmeans>1) {
            savecluster <- list(); savecluster[[1]] <- maxRESULT$cluster; savecounter <- 2
            for (nk in 1:nkmeans) {
                Try <- kmeans(t(Y), g); newclust <- T
                for (m in 1:(savecounter-1)) {
                    if (error.rate(savecluster[[m]], Try$cluster)<0.1) newclust <- F;
                }
                if (!newclust) next;
                savecluster[[savecounter]] <- Try$cluster; savecounter <- savecounter + 1
                if(any(Try$size < k)) Try$logL <- -Inf
                else {
                    Try <- init_param(Try, g, Y, k, q, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag, ndelta)
                }
                if(Try$logL > maxLL) {maxRESULT <- Try; maxLL <- Try$logL} 
                #cat("... Initialisation with trial", savecounter-1, " kmeans, logL = ", Try$logL, " ...\n")
            }
        } else #cat("... Initialisation logL = ", Try$logL, " ...\n")
        INITIAL <- maxRESULT
    } else { 
        maxRESULT <- list(); maxRESULT$centers <-matrix(0,g,k)
        maxRESULT$cluster <- clust    
        maxRESULT <- init_param(maxRESULT, g, Y, k, q, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag, ndelta)
        maxLL <- maxRESULT$logL
        #cat("... Initialisation with given clustering: logL = ", maxLL, " ...\n")
        savecluster <- list(); savecluster[[1]] <- maxRESULT$cluster; savecounter <- 2
        saveRESULT <- list(); saveRESULT[[1]] <- maxRESULT; saveRESULTcounter <- 2
    }
        INITIAL <- maxRESULT     
        CLUST <- maxRESULT$cluster
        if(maxLL == -Inf) stop("Initialization failed. Sample size is too small. \n")
        aVec <- seq(0.1, 0.9, by=0.1)
        for (it in 1:length(aVec)) {
            a <- aVec[it]
            if (g==1) {     
                problem <- F
                if(k==1) {
                    sampleMEAN <- mean(Y)
                    sampleCOV <- matrix(var(Y))
                    sampleSKEW <- diagSCOV <- matrix(skewness(Y))
                }else{      
                    sampleMEAN <- t(t(rowSums(Y)))/N
                    sampleCOV <- cov(t(Y)); diagSCOV <- diag(k); diag(diagSCOV) <- diag(sampleCOV)
                    sampleSKEW <- skewness(t(Y))
                } 
                if(PIflag)    PI <- 1
                if(DOFflag)   DOF <- 40    
                if(SIGMAflag) SIGMA[[1]] <- sampleCOV + (a-1)*diagSCOV
                if(DELTAflag) {
                    if(ndelta==4) {DELTA[[i]] <- matrix(0, k, q); DELTA[[i]][,1] <- sqrt((1-a)/(1-2/pi))*matrix(sqrt(diag(sampleCOV)),k,1)*sign(sampleSKEW)}  
                    else DELTA[[1]] <- diag(sqrt((1-a)/(1-2/pi))*matrix(sqrt(diag(sampleCOV)),k,1)*sign(sampleSKEW))  
                }   
                if(MUflag)    MU[[1]] <- sampleMEAN - sqrt(2/pi)*DELTA[[1]]
                if (det(SIGMA[[1]]) <= 0) {problem <- T}
            } else{
                MEAN <- INITIAL$centers
                CLUST <- INITIAL$cluster
                if(PIflag)  PI <- INITIAL$size/N
                if(DOFflag) DOF <- rep(40,g)
                problem <- F
                for (i in 1:g) {
                    selectY <- {if(k==1) t(matrix(Y[CLUST==i],1)) else t(Y[,CLUST==i])}
                    sampleMEAN <- matrix(MEAN[i,],k,1)
                    sampleCOV <- cov(selectY); diagSCOV <- diag(k); diag(diagSCOV)<-diag(sampleCOV)
                    sampleSKEW <- skewness(selectY)
                    if(SIGMAflag)  SIGMA[[i]] <- sampleCOV + (a-1)*diagSCOV
                    if(DELTAflag)  {
                        if(ndelta==4) {DELTA[[i]] <- matrix(0, k, q); DELTA[[i]][,1] <- sqrt((1-a)/(1-2/pi))*matrix(sqrt(diag(sampleCOV)),k,1)*sign(sampleSKEW)}  
                        else DELTA[[i]] <- diag(sqrt((1-a)/(1-2/pi))*matrix(sqrt(diag(sampleCOV)),k,1)*sign(sampleSKEW))                 
                    }
                    if(MUflag)     MU[[i]] <- sampleMEAN - sqrt(2/pi)*DELTA[[i]]
                    if (det(SIGMA[[i]]) <= 0) {problem <- T; break;}
                }
            }
            if(problem) {LL <- -Inf; TAU <- NA}
            else {
                TT <- try(fust_computeTAU(g, t(Y), MU, SIGMA, DELTA, PI, DOF), silent=T)
                if(!is.numeric(TT)) {LL <- -Inf; TAU <- NA}
                else {LL <- TT$logL; TAU <- TT$TAU}
            }
            result <- list("MU"=MU, "SIGMA"=SIGMA, "DELTA"=DELTA, "PI"=PI, "DOF"=DOF, "logL"=LL, "TAU"=TAU)
            if (LL > maxLL) {maxRESULT <- result; maxLL <- LL}  
            #cat("... Initialisation with a =", a, ", logL = ", LL, " ...\n")
        }        
    maxRESULT$mu <- maxRESULT$MU; maxRESULT$sigma <- maxRESULT$SIGMA; maxRESULT$delta <- maxRESULT$DELTA; maxRESULT$dof <- maxRESULT$DOF; maxRESULT$pro <- maxRESULT$PI;     maxRESULT$fflag <- list("MU"=!is.null(fixed$mu), "SIGMA"=!is.null(fixed$sigma), "DELTA"=!is.null(fixed$delta), "PI"=!is.null(fixed$pro), "DOF"=!is.null(fixed$dof))
    return(maxRESULT)
}



checkMU <- function(g, p, MU) {
    pass <- T
    if(length(MU) !=g)   pass <- F
    for(i in 1:g) if(length(as.numeric(MU[[i]])) != p) pass <- F
    return(pass)
}

checkDELTA <- function(g, p, q, DELTA) {
    pass <- T
    if(length(DELTA) !=g)   pass <- F
    for(i in 1:g) if(length(as.numeric(DELTA[[i]])) != p*q) pass <- F
    return(pass)
}

checkSIGMA <- function(g, p, SIGMA) {
    pass <- T
    if(length(SIGMA) !=g)   pass <- F
    for(i in 1:g) if(any(dim(SIGMA[[i]]) != c(p,p))) pass <- F
    for(i in 1:g) if(any(det(SIGMA[[i]]) == 0)) pass <- F
    for(i in 1:g) if(any(det(SIGMA[[i]]) < 0)) pass <- F
    return(pass)
}

checkDOF <- function(g, DOF) {
    pass <- T
    if(length(DOF) !=g)   pass <- F
    if(any(DOF < 0)) pass <- F
    return(pass)
}

checkPI <- function(g, PI, tol=1e-4) {
    pass <- T
    if(length(PI) !=g)   pass <- F
    if(sum(PI) < 1-tol)  pass <- F
    if(sum(PI) > 1+tol)  pass <- F
    return(pass)
}

