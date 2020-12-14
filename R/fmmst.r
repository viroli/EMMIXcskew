#
#  EM algorithm for Mixture of Unrestricted Multivariate Skew t-distributioins
#  Package: EMMIX-uskew
#  Version: 0.11-6
#
#  Code by S. X. Lee
#  Updated on 31 Jul 2014
#
# Lee S. and Mclachlan, G.J. (2010) On the fitting of finite mixtures of
#   multivariate skew t-distributions via the EM algorithm
#

################################################################################
#  SECTION 1
#             Fitting Mixtures of Multivariate Skew t Distributions
#
################################################################################


fmmst <- function(g=1, dat, initial=NULL, known=NULL, itmax=100, eps=1e-3, clust=NULL, nkmeans=20, print=T) {

    if(!is.matrix(dat)) Y <- as.matrix(dat)  else Y <- dat
    p <- ncol(Y); n <- nrow(Y); fulldebug=F
    if (itmax > 1000) itmax <- 1000   #do not allow more than 1000 iterations (too much)
    if(n > 5000 && p>=3) {
#        cat("  NOTE: For large datasets, please consider using EmSkew for faster computation.\n")
        fulldebug = T
    }
    if(n <= 20*g) stop("sample size is too small!")

    if(print) {
        cat("Finite Mixture of Multivariate Skew t-distribution\n")
        if(g<=1) cat("with 1 component\n")
        else cat("with ", g, "components\n")
        cat("  ----------------------------------------------------\n\n")
    }
    initial <- init(g, t(Y), initial, known, clust, nkmeans)
    return(fmmstt(g, p, n, Y, initial, known, nkmeans, itmax, eps, print, fulldebug)) 
    
}


init <- function(g, Y, initial=NULL, fixed=NULL, clust=NULL, nkmeans=100) {

    
    init_param1 <- function(Fit, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag) {
        Fit$MU <- Fit$mu <- Fit$SIGMA <- Fit$sigma <- Fit$delta <- Fit$DELTA <- list()
        for (i in 1:g) {  #calculate initial parameter values
            if(k==1){
                if(MUflag)  Fit$MU[[i]] <- mean(Y[Fit$cluster==i]) else Fit$Mu[[i]] <- MU[[i]]
                if(SIGMAflag) Fit$SIGMA[[i]] <- var(Y[Fit$cluster==i]) else Fit$SIGMA[[i]] <- SIGMA[[i]]
                skew <- skewness(Y[Fit$cluster==i]) 
            }else{            
                if(MUflag) Fit$MU[[i]] <- matrix(rowSums(Y[,Fit$cluster==i])/length(Fit$cluster[Fit$cluster==i]),k,1) else Fit$MU[[i]] <- MU[[i]]
                if(SIGMAflag) Fit$SIGMA[[i]] <- cov(t(Y[,Fit$cluster==i])) else Fit$SIGMA[[i]] <- SIGMA[[i]]
                skew <- skewness(t(Y[,Fit$cluster==i]))
            }
            if(DELTAflag) Fit$DELTA[[i]] <- -5*(skew <(-0.1)) + 5*(skew>0.1)  else Fit$DELTA[[i]] <- DELTA[[i]]
            Fit$mu[[i]] <- Fit$MU[[i]]
            Fit$sigma[[i]] <- Fit$SIGMA[[i]]
            Fit$delta[[i]] <- Fit$DELTA[[i]]
        }
        if(DOFflag)  Fit$DOF <- Fit$dof <- rep(4, g) else Fit$dof <- Fit$DOF <- DOF
        Fit$PI <- Fit$pro <- Fit$size/N
        tmp <- computet(g, t(Y), Fit$MU, Fit$SIGMA, Fit$DELTA, Fit$PI, Fit$DOF)
        Fit$TAU <- tmp$TAU; Fit$logL <- tmp$logL
        return(Fit)
    }
    
    
    init_param2 <- function(Y, k, N, g, CLUST, a, iDOF, MUflag, SIGMAflag, DELTAflag, DOFflag, PIflag){
        problem <- F
        MU <- SIGMA <- DELTA <- list()
       if (g==1) {   
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
            if(DOFflag)   DOF <- iDOF
            if(SIGMAflag) SIGMA[[1]] <- sampleCOV + (a-1)*diagSCOV
            if(DELTAflag) DELTA[[1]] <- sqrt((1-a)/(1-2/pi))*matrix(sqrt(diag(sampleCOV)),k,1)*sign(sampleSKEW)
            if(MUflag)    MU[[1]] <- sampleMEAN - sqrt(2/pi)*DELTA[[1]]
            if (det(SIGMA[[1]]) < 0) {problem <- T}
    
        } else{
            if(PIflag)  PI <- as.numeric(table(CLUST))/N
            if(DOFflag) DOF <- rep(iDOF,g)
            for (i in 1:g) {
                if(k==1) {
                    sampleMEAN <- mean(Y[CLUST==i])
                    sampleCOV <- matrix(var(Y[CLUST==i]))
                    sampleSKEW <- diagSCOV <- matrix(skewness(Y[CLUST==i]))
                } else{  
                    sampleMEAN <- matrix(rowSums(Y[,CLUST==i],k,1))/length(CLUST[CLUST==i])
                    sampleCOV <- cov(t(Y[,CLUST==i])); diagSCOV <- diag(diag(sampleCOV))
                    sampleSKEW <- skewness(t(Y[,CLUST==i]))
                }
                if(SIGMAflag)  SIGMA[[i]] <- sampleCOV + (a-1)*diagSCOV                                                   #kxk
                if(DELTAflag)  DELTA[[i]] <- sqrt((1-a)/(1-2/pi))*matrix(sqrt(diag(sampleCOV)),k,1)*sign(sampleSKEW)     #kx1
                if(MUflag)     MU[[i]] <- sampleMEAN - sqrt(2/pi)*DELTA[[i]]                                              #kx1
                if (det(SIGMA[[i]]) < 0) {problem <- T; break;}
            }
        }
        tmp <- try(computet(g, t(Y), MU, SIGMA, DELTA, PI, DOF),silent=T)
        if(!is.list(tmp)) {tmp$logL <- -Inf; problem <- T}
        return(list("problem"=problem, "MU"=MU, "SIGMA"=SIGMA, "DELTA"=DELTA, "DOF"=DOF, "PI"=PI, "logL"=tmp$logL, "TAU"=tmp$TAU, "mu"=MU, "sigma"=SIGMA, "delta"=DELTA, "dof"=DOF, "pro"=PI))
    }    
    
    checkDELTA <- function(g, p, DELTA) {
        pass <- T
        if(length(DELTA) !=g)   pass <- F
        for(i in 1:g) if(length(as.numeric(DELTA[[i]])) != p) pass <- F
        return(pass)
    }    

    P <- dim(Y);  k <- P[1];  N <- P[2]
    w <- options("warn"); on.exit(options(w)); options(warn=-1)
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
        cat("WARNING: in fmmst initialisation, pro is not correctly specified.\n")
        PIflag <- TRUE
    }
    if(!DOFflag) if(!checkDOF(g, DOF)) {
        cat("WARNING: in fmmst initialisation, dof is not correctly specified.\n")
        DOFflag <- TRUE
    }
    if(!MUflag) if(!checkMU(g, k, MU)) {
        cat("WARNING: in fmmst initialisation, mu is not correctly specified.\n")
        MUflag <- TRUE
    }
    if(!SIGMAflag) if(!checkSIGMA(g, k, SIGMA)) {
        cat("WARNING: in fmmst initialisation, sigma is not correctly specified.\n")
        SIGMAflag <- TRUE
    }
    if(!DELTAflag) if(!checkDELTA(g, k, DELTA)) {
        cat("WARNING: in fmmst initialisation, delta is not correctly specified.\n")
        DELTAflag <- TRUE
    }
    if (MUflag || SIGMAflag || DELTAflag || PIflag || DOFflag) {
        if(MUflag) MU <- list()
        if(SIGMAflag) SIGMA <- list()
        if(DELTAflag) DELTA <- list()
        if(PIflag) PI <- vector()
        if(DOFflag) DOF <- vector()
        if(is.null(clust)) {
            Try <- kmeans(t(Y), g);
            Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
            maxLL <- Try$logL
            maxRESULT <- Try                      
            clusters <- list(); clusters[[1]] <- maxRESULT$cluster;
            if(g > 1 && nkmeans>1) {
                for (nk in 1:nkmeans) {
                    Try <- kmeans(t(Y), g); newclust <- T
                    for (m in 1:length(clusters)) {
                        if (error.rate(clusters[[m]], Try$cluster)<0.1) newclust <- F;
                    }
                    if (!newclust) next;
                    if(any(Try$size < k)) Try$logL <- -Inf
                    else {
                        clusters[[length(clusters)+1]] <- Try$cluster;
                        Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
                    }
                    if(Try$logL > maxLL) {maxRESULT <- Try; maxLL <- Try$logL}
                }
            }

        } else{
            Try <- list();  #Try$MU <- Try$DELTA <- Try$SIGMA <- list()
            if(any(Try$size < k)) stop("Some components have less than",k,"observations.\n")
            Try$cluster <- clust; Try$size <- as.numeric(table(clust))
            Try <- init_param1(Try, g, Y, k, N, MU, SIGMA, DELTA, DOF, MUflag, SIGMAflag, DELTAflag, DOFflag)
            maxLL <- Try$logL
            maxRESULT <- Try
        }

        INITIAL <- maxRESULT
        if(maxLL == -Inf) stop("Initialization failed. Sample size is too small. \n")
        aVec <- seq(0.1, 0.9, by=0.1)
        for (it in 1:length(aVec)) {
            a <- aVec[it]
            Try <- init_param2(Y, k, N, g, INITIAL$cluster, a, 40, MUflag, SIGMAflag, DELTAflag, DOFflag, PIflag)
            if(Try$problem) LL <- -Inf else LL <- Try$logL
            if (LL > maxLL) {maxRESULT <- Try; maxLL <- LL}
        }
    } else {
        logL <- sum(log(dfmmst(t(Y), MU, SIGMA, DELTA, DOF, PI)))
        maxRESULT <- list("MU"=MU, "SIGMA"=SIGMA, "DELTA"=DELTA, "PI"=PI, "DOF"=DOF, "logL"=logL)
    }
    maxRESULT$fflag <- list("MU"=!is.null(fixed$mu), "SIGMA"=!is.null(fixed$sigma), "DELTA"=!is.null(fixed$delta), "PI"=!is.null(fixed$pro), "DOF"=!is.null(fixed$dof))
    return(maxRESULT)
    
}
    



fmmstt <- function(g=1, p=1, n=1, Y, initial=NULL, known=NULL, nkmeans=20, itmax=100, eps=1e-6, debug=T, fulldebug=F) {
    N <- n
    MU <- initial$MU
    SIGMA <- initial$SIGMA
    DELTA <- initial$DELTA
    PI <- initial$PI
    DOF <- initial$DOF
    fflag <- initial$fflag   
    if(fflag$MU && fflag$SIGMA && fflag$DELTA && fflag$DOF && fflag$PI) {        
        known <- list("mu"=MU, "sigma"=SIGMA, "delta"=DELTA, "pro"=PI, "dof"=DOF) 
        tmp <- computet(g, Y, MU, SIGMA, DELTA, PI, DOF)
        known$tau <- tmp$TAU 
        known$clusters <- apply(known$tau,2,which.max)
        known$loglik <- known$lk <- tmp$logL  
        m <- g*(p + p + 0.5*p*(p+1)) + (g-1) + g   
        known$aic <- 2*m - 2*known$loglik
        known$bic <- m*log(N) - 2*known$loglik
        cat("NOTE: All parameters are known. uskew will terminate.\n");
        return(known)
    }
    if(fulldebug) cat("  ... Initialisation completed ...\n")
    TAU <- eta <- E1 <- E2 <- matrix(0, g, N); E3 <- E4 <- list();
    k <- 1; epsilon <- Inf; 
    m <- g*(p + p + 0.5*p*(p+1)) + (g-1) + g
    TAU <- computet(g, Y, MU, SIGMA, DELTA, PI, DOF)$TAU
    LL <- lk <- initial$logL
    aic <- 2*m - 2*LL; bic <- m*log(n) - 2*LL;
    singular <- F;  tauCutOff <- 5e-8; 
    sumTAU <- rowSums(TAU)                                         
    if(any(sumTAU < p)) stop("one or more of the cluster is too small. Please consider reducing g.\n")  #full options not available for GPL-ed version  
    if(debug) cat("  Iteration  0 : loglik = ", LL, "\n")
 
    while((k <= itmax) && (epsilon > eps)) {            
        saveDOF <- DOF;
        for(i in 1:g) {
            E3[[i]] <- E4[[i]] <- list()
            DD <- diag(p); diag(DD)<-as.numeric(DELTA[[i]])                 
            OMEGA <- SIGMA[[i]] + DD %*% DD                               
            invOMEGA <- solve(OMEGA)                                      
            LAMBDA <- diag(p) - DD%*%invOMEGA%*%DD                 
            Q <- DD%*%invOMEGA%*%(t(Y)-MU[[i]]%*%matrix(1,1,N))           
            eta[i,] <- mahalanobis(Y, as.numeric(MU[[i]]), invOMEGA, T)   
            Q1 <- Q*(matrix(1,p,1)%*%sqrt((DOF[i]+p+2)/(DOF[i]+eta[i,]))) 
            Q2 <- Q*(matrix(1,p,1)%*%sqrt((DOF[i]+p)/(DOF[i]+eta[i,])))   
            for(j in 1:n) {
                if (TAU[i,j] < tauCutOff) {E2[i,j] <- 0; E3[[i]][[j]] <- matrix(0,p,1); E4[[i]][[j]] <- matrix(0,p,p); next;}
                T1 <- pmt(Q1[,j], rep(0,p), LAMBDA, round(DOF[i])+p+2)
                T2 <- pmt(Q2[,j], rep(0,p), LAMBDA, round(DOF[i])+p)                                   
                if (T2==0) {cat("T2=0 at ",i,j,"; tau =",TAU[i,j],"\n"); E2[i,j] <- 0; E3[[i]][[j]] <- matrix(0,p,1); E4[[i]][[j]] <- matrix(0,p,p); next;}
                E2[i,j] <- (DOF[i]+p)/(DOF[i]+eta[i,j]) * T1 / T2         
                TT <- .dds(rep(0,p), Q[,j], (DOF[i]+eta[i,j])/(DOF[i]+p+2)*LAMBDA, round(DOF[i])+p+2)
                E3[[i]][[j]] <- TT$EX * E2[i,j]                           
                E4[[i]][[j]] <- TT$EXX * E2[i,j]                                   
                if(any(is.nan(E3[[i]][[j]]))||any(E3[[i]][[j]]==Inf)||any(is.nan(E4[[i]][[j]]))||any((E4[[i]][[j]]==Inf))) {E2[i,j] <- 0; E3[[i]][[j]] <- matrix(0,p,1); E4[[i]][[j]] <- matrix(0,p,p)}   
            }              
            if(fulldebug) cat("  ... E-step for component ", i, " completed ...\n",sep="")
        }    
        sumTAU <- rowSums(TAU)                                         
        if(any(sumTAU < p)) {
            warning("  uskew have detected one or more of the cluster is getting too small.\n")
            break;  #full options not available for GPL-ed version
        }
        for (i in 1:g) {
            DD <- diag(p); diag(DD)<-as.numeric(DELTA[[i]])              
            invSIGMA <- solve(SIGMA[[i]])
            if(!fflag$PI) PI[i] <- sumTAU[i]/N 
            M1 <- colSums(E2[i,]*TAU[i,]*Y)       
            M2 <- sum(E2[i,]*TAU[i,])            
            M3 <- matrix(0,p,1)
            M4 <- M5 <- M6 <- matrix(0,p,p)
            for (j in 1:n) {
              M3 <- M3 + TAU[i,j]*E3[[i]][[j]]    
              M4 <- M4 + TAU[i,j]*E4[[i]][[j]]    
            }
            if(!fflag$MU) MU[[i]] <- (M1 - DD%*%M3) / M2       
            for(j in 1:n) {
                M5 <- M5 + TAU[i,j]*(matrix(Y[j,],p)-MU[[i]])%*%t(E3[[i]][[j]])                  
                M6 <- M6 + TAU[i,j]*E2[i,j]*(matrix(Y[j,],p)-MU[[i]])%*%(Y[j,]-matrix(MU[[i]],1)) 
            }
            if(!fflag$DELTA) {
                A <- try(solve(invSIGMA * M4) %*% diag(invSIGMA%*%M5),silent=T)                   
                if(is.numeric(A)) DELTA[[i]] <- A  else singular <- T
            }
            DD <- diag(p); diag(DD)<-as.numeric(DELTA[[i]])                                       
            Den <- sumTAU[i]
            if(!fflag$SIGMA) SIGMA[[i]] <- (DD%*%M4%*%DD - DD%*%t(M5) - M5%*%DD + M6)/Den        
            if(!fflag$DOF) {
                Num <- sum(TAU[i,]*(digamma(0.5*(DOF[i]+p))-log(0.5*(DOF[i]+eta[i,]))-(DOF[i]+p)/(DOF[i]+eta[i,])))          
                DOFfun <- function(v) {log(v/2)-digamma(v/2)+1 + Num/Den}
                if(sign(DOFfun(1)) == sign(DOFfun(200)))  DOF[i] <- {if(abs(DOFfun(1)) < abs(DOFfun(400))) 1 else 200}  
                else DOF[i] <- uniroot(DOFfun, c(1,200))$root                   
            }     
        }  
        if(fulldebug) cat("  ... M-step completed ...\n") 
        if(singular) {k <- k+1; break;}
        else tmp <- computet(g, Y, MU, SIGMA, DELTA, PI, DOF)
        if(tmp$logL < LL) {DOF <- saveDOF; tmp <- computet(g, Y, MU, SIGMA, DELTA, PI, DOF)}                
        TAU <- tmp$TAU; newLL <- tmp$logL; lk <- c(lk,newLL)    
        if(debug) cat("  Iteration ",k,": loglik = ",newLL, "\n")
        if (k < 2) epsilon <- abs(LL-newLL)/abs(newLL)
        else {
            tmp <- (newLL - LL)/(LL-lk[length(lk)-1])
            tmp2 <- LL + (newLL-LL)/(1-tmp)
            epsilon <- abs(tmp2-newLL)
        }  
        LL <- newLL
        k <- k+1         
    }                         
    m <- g*(p + p + 0.5*p*(p+1)) + (g-1) + g
    aic <- 2*m - 2*LL; bic <- m*log(n) - 2*LL;    
    clusters <- apply(TAU,2,which.max)
    results <- list("pro"=PI, "mu"=MU, "sigma"=SIGMA, "delta"=DELTA, "dof"=DOF, 
            "tau"=TAU, "clusters"=clusters, "loglik"=LL, "lk"=lk, 
            "iter"=(k-1), "eps"=epsilon, "aic"=aic, "bic"=bic)
    attr(results, "class") <- "fmmst" 
    if(debug) {
        cat("  ----------------------------------------------------\n")
        cat("  Iteration ", k-1, ": loglik = ", LL,"\n\n",sep="")
        if(singular) cat("\nNOTE: Sigma is computationally singular at iteration",k-1,"\n\n")
        summary2(results)
    }
    return(results)     
    
}               


computet <- function(g, Y, MU, SIGMA, DELTA, PI, DOF) {
    n <- nrow(Y); p <- ncol(Y); 
    if(is.null(n)) {n <- 1; p <- length(Y)}   
    if (g == 1) {
        TAU <- matrix(1, g, n)
        logL <- sum(log(dmst(Y, MU[[1]], SIGMA[[1]], DELTA[[1]], DOF[1])))
    } else {
        TAU <- matrix(0,g, n)
        for (i in 1:g) TAU[i,] <- PI[i]*dmst(Y, as.numeric(MU[[i]]), SIGMA[[i]], DELTA[[i]], DOF[i]) 
        logL <- sum(log(colSums(TAU)))
        sumTAU <- matrix(1, g, 1) %*% matrix(colSums(TAU),1,n) 
        TAU <- TAU/sumTAU
    }
    return(list(TAU=TAU,logL=logL))
}


.dds <- function(a, mu, sigma, dof) {
    p <- nrow(sigma)
    Tp <- 1
    ET <- .ds2(sigma, dof, mu-a)
    mu <- matrix(mu,p)
    EX <- mu - ET$EX
    EXX <- ET$EXX - ET$EX%*%t(mu) - mu%*%t(ET$EX) + mu%*%t(mu)
    return(list(EX=EX/Tp, EXX=EXX/Tp))
}

.ds2 <- function(sigma, dof, a=NULL) {
    p <- dim(sigma)[1]; if(is.null(p)) {p<-1}
    if (is.null(a)) a<- rep(0,p); a <- as.numeric(a)
    if (p==1) {
        Tp <- pt(a/sqrt(sigma),df=dof)
        EX <- -0.5*sqrt(sigma*dof/pi)*(1+a^2/dof/sigma)^(-0.5*dof+0.5)*gamma(0.5*dof-0.5)/gamma(0.5*dof)
        EXX <- sigma*dof/(dof-2)*Tp + (dof-1)/(dof-2)*a*EX
        return(list("EX"=EX/Tp, "EXX"=EXX/Tp)); return(list("EX"=EX/Tp))
    }else{
        Tp <- pmt(a, rep(0,p), sigma, dof)
        xi <- ai <- rep(0,p); theta <- matrix(0, p, p); G1 <- gamma(0.5*dof-0.5)/gamma(0.5*dof)
        V1 <- sqrt(0.5*dof); V2 <- dof^(0.5*dof)/(dof-2)
        for (i in 1:p){
            Sii <- sigma[i,i]                                                                                            
            S1 <- (dof+a[i]^2/Sii)/(dof-1)*(sigma[-i,-i]-sigma[-i,i]%*%t(sigma[-i,i])/Sii)                               
            A1 <- a[-i] - a[i]/Sii*sigma[-i,i]                                                                        
            xi[i] <- (dof/(dof+a[i]^2/Sii))^(0.5*dof-0.5)/sqrt(2*pi*Sii)*G1*V1*pmt(A1, rep(0,p-1), S1, dof-1)[1]       
            if(p==2) {                                                     
                if (i == 1) {                                              
                    S2 <- solve(sigma)                                     
                    V3 <- dof + a %*% S2 %*% a                             
                    theta[1,2] <- -V2/(V3^(0.5*dof-1))/(2*pi*sqrt(Sii*sigma[2,2]-sigma[1,2]^2))
                    theta[2,1] <- theta[1,2]                               
                }                                                          
            }else{                                                         
                if (i != p) {                                              
                    for (j in (i+1):p) {                                   
                        S2 <- solve(matrix(c(Sii, sigma[i,j], sigma[i,j], sigma[j,j]),2,2))                                 
                        V3 <- dof + a[c(i,j)] %*% S2 %*% a[c(i,j)]                                                            
                        S3 <- as.numeric(V3/(dof-2))*(sigma[-c(i,j),-c(i,j)]-matrix(sigma[-c(i,j),c(i,j)],p-2,2)%*%S2%*%t(matrix(sigma[-c(i,j),c(i,j)], p-2, 2))) 
                        A2 <- a[-c(i,j)] - sigma[-c(i,j),c(i,j)]%*%S2%*%a[c(i,j)]                                          
                        theta[i,j] <- -V2/(V3^(0.5*dof-1))/(2*pi*sqrt(Sii*sigma[j,j]-sigma[i,j]^2))*pmt(as.numeric(A2), rep(0,p-2),S3, dof-2)
                        theta[j,i] <- theta[i,j]                           
                    }                                                      
                }                                                          
            }                                                              
            theta[i,i] <- a[i]*xi[i]/Sii - sum(sigma[i,-i]*theta[i,-i])/Sii
        } 
    } 
    EX <- -sigma%*%xi/Tp       
    EXX <- dof/(dof-2)*pmt(a,rep(0,p),dof/(dof-2)*sigma,dof-2)*sigma - sigma%*%theta%*%sigma
    EXX <- EXX/Tp           #pxp
    return(list("EX"=EX, "EXX"=EXX))
}


summary.fmmst <- function(object, ...) {
    g <- length(object$dof)
    p <- nrow(object$sigma[[1]]) 
    if(is.null(p)) p <- 1
    if(g<1) stop("invalid fmmst object: g < 1")
    cat("Finite Mixture of Multivarate Skew t-distributions\n")
    if(g==1) cat("with ", g, " component\n\n")
    else cat("with ", g, " components\n\n")
    if(g==1) cat("Mean:\n")
    else cat("Component means:\n")
    means <- deltas <- matrix(0,p,g)
    for (i in 1:g) {
        means[,i] <- matrix(object$mu[[i]],p,1)
        deltas[,i] <- matrix(object$delta[[i]],p,1)
    }
    print(means)
    cat("\n\n")
    if(g==1) cat("Scale matrix:\n")
    else cat("Component scale matrices:\n")
    print(object$sigma)
    if(g==1) cat("\nSkewness parameter:\n")
    else cat("\nComponent skewness parameters:\n")
    print(deltas)
    cat("\n\n")
    if(g==1) cat("Degrees of freedom:\n")
    else cat("Component degrees of freedom:\n")
    cat(object$dof, "\n\n")
    if(g>1) {
      cat("Component mixing proportions:\n")
      cat(object$pro, "\n\n")
    }
}

summary2 <- function(x) {
    g <- length(x$pro)
    p <- nrow(x$sigma[[1]])       
    if(is.null(p)) p <- 1
    if(g==1) cat("Mean:\n")
    else cat("Component means:\n")
    means <- deltas <- matrix(0,p,g)
    for (i in 1:g) {
        means[,i] <- matrix(x$mu[[i]],p,1)
        deltas[,i] <- matrix(x$delta[[i]],p,1)
    }
    print(means)
    cat("\n\n")
    if(g==1) cat("Scale matrix:\n")
    else cat("Component scale matrices:\n")
    print(x$sigma)
    if(g==1) cat("\nSkewness parameter:\n")
    else cat("\nComponent skewness parameters:\n")
    print(deltas)
    cat("\n\n")
    if(g==1) cat("Degrees of freedom:\n")
    else cat("Component degrees of freedom:\n")
    cat(x$dof, "\n\n")
    if(g>1) {
      cat("Component mixing proportions:\n")
      cat(x$pro, "\n\n")
    }
}

print.fmmst <- function(x, ...) {
    g <- length(x$dof)
    cat("Finite Mixture of Multivarate Skew t-distribution\n")
    if(g==1) cat("with ", g, " component\n\n")
    else cat("with ", g, " components\n\n")
    obj <- as.fmmst(g, x$mu, x$sigma, x$delta, x$dof, x$pro)
    obj$tau <- x$tau
    obj$clusters <- x$clusters
    obj$loglik <- x$loglik
    obj$aic <- x$aic
    obj$bic <- x$bic
    print(obj)
}

as.fmmst <- function(g, mu, sigma, delta, dof=rep(1,g), pro=rep(1/g,g)) {
    obj <- list()
    if(missing(mu)) stop("mu must be specified!") 
    if(g==1){
        obj$mu <- if(is.list(mu)) mu[[1]] else mu
        p <- length(as.numeric(obj$mu))
        if(missing(sigma)) sigma <- diag(p)
        if(missing(delta)) delta <- rep(0,p)
        obj$sigma <- if(is.list(sigma)) sigma[[1]] else sigma 
        obj$delta <- if(is.list(delta)) delta[[1]] else delta 
        obj$dof <- ifelse(length(dof)>1, dof[1], dof)
        obj$pro <- 1 
    }
    else {
        obj$mu <- mu  
        obj$sigma <- sigma
        obj$delta <- delta
        obj$dof <- dof
        obj$pro <- pro
    }
    return(obj)
}
