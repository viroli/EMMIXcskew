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
#  fmcfust.r
#             EM algorithm for fitting FM-CFUST
#
################################################################################

fmcfust <- function(g=1, dat, q, initial=NULL, known=NULL, clust=NULL, itmax=100, eps=1e-6, nkmeans=20, verbose=T, method=c("moments","transformation","EMMIXskew","EMMIXuskew"), convergence=c("Aitken","likelihood","parameters")) {

    if(!is.matrix(dat)) Y <- as.matrix(dat)  else Y <- dat
    p <- ncol(Y); n <- nrow(Y); fulldebug=F; 
    if(missing(q)) q <- p
    if (itmax > 1000) itmax <- 1000   #do not allow more than 1000 iterations (too much)
    if(n > 5000 && p>=3) {
#        cat("  NOTE: For large datasets, please consider using EmSkew for faster computation.\n")
        fulldebug = T
    }
    if(n <= 20*g) stop("sample size is too small!")      

    if(verbose) {
        cat("Finite Mixture of Multivariate CFUST Distributions\n")
        if(g<=1) cat("with 1 component\n")
        else cat("with ", g, "components\n")
        cat("  ----------------------------------------------------\n\n")
    }  
    initial <- init.cfust(g, Y, q, initial, known, clust, nkmeans, method=method)
    ndelta <- 1; if(q==1) ndelta <- 4
    res <- fmfustt(g, p, q, n, Y, initial, known, ndelta, itmax, eps, verbose, fulldebug, convergence=convergence)
    return(res)
}


fmfustt <- function(g=1, p=1, q, n=1, Y, initial=NULL, known=NULL, ndelta, itmax=100, eps=1e-6, debug=T, fulldebug=F, convergence=convergence) {
    N <- n
    MU <- initial$MU
    SIGMA <- initial$SIGMA
    DELTA <- initial$DELTA
    PI <- initial$PI
    DOF <- initial$DOF
    fflag <- initial$fflag
    if(fflag$MU && fflag$SIGMA && fflag$DELTA && fflag$DOF && fflag$PI) {
        known <- list("mu"=MU, "sigma"=SIGMA, "delta"=DELTA, "pro"=PI, "dof"=DOF)
        tmp <- fust_computeTAU(g, Y, MU, SIGMA, DELTA, PI, DOF)
        known$tau <- tmp$TAU
        known$clusters <- apply(known$tau,2,which.max)
        known$loglik <- known$lk <- tmp$logL
        m <- g*(p + (p*q) + 0.5*p*(p+1)) + (g-1) + g
        known$aic <- 2*m - 2*known$loglik
        known$bic <- m*log(N) - 2*known$loglik
        cat("NOTE: All parameters are known. fmcfust will terminate.\n");
        return(known)
    }
    if(fulldebug) cat("  ... Initialisation completed ...\n")
    TAU <- eta <- E1 <- E2 <- matrix(0, g, N); E3 <- E4 <- list();
    k <- 1; epsilon <- Inf;      
    m <- g*(p*(!fflag$MU) + (p*q)*(!fflag$DELTA) + 0.5*p*(p+1)*(!fflag$SIGMA)) + g*(!fflag$DOF) + (g-1)*(!fflag$PI)
    TAU <- initial$TAU
    LL <- lk <- initial$logL
    aic <- aicVec <- 2*m - 2*LL; bic <- bicVec <- m*log(n) - 2*LL;
    singular <- F;  tauCutOff <- 5e-8;  cvg <- convergence[1]
    sumTAU <- rowSums(TAU)
    if(any(sumTAU < p)) cat("one or more of the cluster is too small. Please consider reducing g.\n")  #full options not available for GPL-ed version
    if(debug) cat("  Iteration  0 : loglik = ", LL, "\n")
    E1s <- E1 <-  matrix(0, g, N)

    while((k <= itmax) && (epsilon > eps)) {
        saveDOF <- DOF;  
        for(i in 1:g) {
            E3[[i]] <- E4[[i]] <- list()
            DD <- DELTA[[i]]                                      
            OMEGA <- SIGMA[[i]] + DD %*% t(DD)                    
            invOMEGA <- solve(OMEGA)                              
            LAMBDA <- diag(q) - t(DD)%*%invOMEGA%*%DD             
            Q <- t(DD)%*%invOMEGA%*%(t(Y)-MU[[i]]%*%matrix(1,1,N))
            eta[i,] <- mahalanobis(Y, as.numeric(MU[[i]]), invOMEGA, T)   
            Q1 <- Q*(matrix(1,q,1)%*%sqrt((DOF[i]+p+2)/(DOF[i]+eta[i,]))) 
            Q2 <- Q*(matrix(1,q,1)%*%sqrt((DOF[i]+p)/(DOF[i]+eta[i,])))   
            for(j in 1:n) {
                if (TAU[i,j] < tauCutOff) {E2[i,j] <- 0; E3[[i]][[j]] <- matrix(0,q,1); E4[[i]][[j]] <- matrix(0,q,q); next;}
                T1 <- pmt(Q1[,j], rep(0,q), LAMBDA, round(DOF[i])+p+2)
                T2 <- pmt(Q2[,j], rep(0,q), LAMBDA, round(DOF[i])+p)
                if (T2==0) {cat("T2=0 at ",i,j,"; tau =",TAU[i,j],"\n"); E2[i,j] <- 0; E3[[i]][[j]] <- matrix(0,q,1); E4[[i]][[j]] <- matrix(0,q,q); next;}
                E2[i,j] <- (DOF[i]+p)/(DOF[i]+eta[i,j]) * T1 / T2
                TT <- .dds(rep(0,q), Q[,j], (DOF[i]+eta[i,j])/(DOF[i]+p+2)*LAMBDA, round(DOF[i])+p+2)
                E3[[i]][[j]] <- TT$EX * E2[i,j]
                E4[[i]][[j]] <- TT$EXX * E2[i,j]
                if(any(is.nan(E3[[i]][[j]]))||any(E3[[i]][[j]]==Inf)||any(is.nan(E4[[i]][[j]]))||any((E4[[i]][[j]]==Inf))) {E2[i,j] <- 0; E3[[i]][[j]] <- matrix(0,q,1); E4[[i]][[j]] <- matrix(0,q,q)}
            }
            if(fulldebug) cat("  ... E-step for component ", i, " completed ...\n",sep="")
        }
        if(fulldebug) cat("  ... Iteration ", k, " E-step completed ...\n",sep="") 
        sumTAU <- rowSums(TAU)
        if(cvg=="parameters") par.prev <- list(MU=MU, SIGMA=SIGMA, DELTA=DELTA, DOF=DOF, PI=PI)
        for (i in 1:g) {
            if(sumTAU[i] < p) cat("  NOTE: cluster", i, "is getting too small.\n")
            DD <- DELTA[[i]]
            invSIGMA <- solve(SIGMA[[i]])
            if(!fflag$PI) PI[i] <- sumTAU[i]/N
            M1 <- colSums(E2[i,]*TAU[i,]*Y)
            M2 <- sum(E2[i,]*TAU[i,])
            M3 <- matrix(0,q,1)
            M4 <- matrix(0, q, q)
            M5 <- matrix(0, p, q)
            M6 <- matrix(0, p,p)
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
                if(ndelta==4) {A <- matrix(0, p, q); A[,1] <- M5[,1]/M4[1,1]}     #rMST
                else if(ndelta==3) A <- try(solve(invSIGMA * M4) %*% diag(invSIGMA%*%M5),silent=T)  #uMST
                else A <- try(M5 %*% solve(M4),silent=T)    #general                
                if(is.numeric(A)) DELTA[[i]] <- A  else {singular <- T; browser()}
            }
            DD <- DELTA[[i]]
            Den <- sumTAU[i]
            if(!fflag$SIGMA) SIGMA[[i]] <- (DD%*%t(M4)%*%t(DD) - DD%*%t(M5) - M5%*%t(DD) + M6)/Den
            if(!fflag$DOF) {
                Num <- sum(TAU[i,]*(digamma(0.5*(DOF[i]+p))-log(0.5*(DOF[i]+eta[i,]))-(DOF[i]+p)/(DOF[i]+eta[i,])))
                DOFfun <- function(v) {log(v/2)-digamma(v/2)+1 + Num/Den}
                if(sign(DOFfun(1)) == sign(DOFfun(200)))  DOF[i] <- {if(abs(DOFfun(1)) < abs(DOFfun(400))) 1 else 200}
                else DOF[i] <- uniroot(DOFfun, c(1,200))$root
            }
        }
        if(fulldebug) cat("  ... Iteration ", k, " M-step completed ...\n",sep="") 
        if(singular) {k <- k+1; break;}
        else tmp <- fust_computeTAU(g, Y, MU, SIGMA, DELTA, PI, DOF)
        if(tmp$logL < LL) {DOF <- saveDOF; tmp <- fust_computeTAU(g, Y, MU, SIGMA, DELTA, PI, DOF)}
        TAU <- tmp$TAU; newLL <- tmp$logL; lk <- c(lk,newLL)
        if(debug) cat("  Iteration ",k,": loglik = ",newLL, "\n")
        if(cvg=="Aitken" && k > 2) {
            tmp <- (newLL - LL)/(LL-lk[length(lk)-1])
            tmp2 <- LL + (newLL-LL)/(1-tmp)
            epsilon <- abs(tmp2-newLL)}
        else {if(cvg=="parameters"){
          epsilon <- F
          for(i in 1:g){  
              if(abs(PI[i]-par.prev$PI[i])/abs(PI[i])>eps) {epsilon<-T; break} 
              if(abs(DOF[i]-par.prev$DOF[i])/abs(DOF[i])>eps) {epsilon<-T; break}
              if(all(abs(as.numeric(MU[[i]])-as.numeric(par.prev$MU[[i]]))/abs(as.numeric(MU[[i]]))>eps)) {epsilon<-T; break} 
              if(all(abs(as.numeric(SIGMA[[i]])-as.numeric(par.prev$SIGMA[[i]]))/abs(as.numeric(SIGMA[[i]]))>eps)) {epsilon<-T; break}
              if(all(abs(as.numeric(DELTA[[i]])-as.numeric(par.prev$DELTA[[i]]))/abs(as.numeric(DELTA[[i]]))>eps)) {epsilon<-T; break}   
          }    
          if(epsilon) epsilon<-Inf else epsilon<-0 
          }
          else  epsilon <- abs(LL-newLL)/abs(newLL)
        }
        aic <- 2*m - 2*newLL; bic <- m*log(n) - 2*newLL; aicVec <- c(aicVec, aic); bicVec <- c(bicVec, bic)
        LL <- newLL
        k <- k+1
    }
    aic <- 2*m - 2*LL; bic <- m*log(n) - 2*LL;
    clusters <- apply(TAU,2,which.max)
    results <- list("pro"=PI, "mu"=MU, "sigma"=SIGMA, "delta"=DELTA, "dof"=DOF,
            "tau"=TAU, "clusters"=clusters, "loglik"=LL, "lk"=lk, "iter"=(k-1), 
            "eps"=epsilon, "aic"=aic, "bic"=bic)
    attr(results, "class") <- "fmcfust"
    if(debug) {
        cat("  ----------------------------------------------------\n")
        cat("  Iteration ", k-1, ": loglik = ", LL,"\n\n",sep="")
        #if(singular) cat("\nNOTE: Sigma is computationally singular at iteration",k-1,"\n\n")
    }
    return(results)
}


fust_computeTAU <- function(g, Y, MU, SIGMA, DELTA, PI, DOF) {
    n <- nrow(Y); p <- ncol(Y)          
    if(is.null(n)) {n <- 1; p <- length(Y)}   
    if (g == 1) {
        TAU <- matrix(1, g, n)
        logL <- sum(log(dcfust(Y, MU[[1]], SIGMA[[1]], DELTA[[1]], DOF[1])))
    } else {
        TAU <- matrix(0,g, n)
        for (i in 1:g) TAU[i,] <- PI[i]*dcfust(Y, as.numeric(MU[[i]]), SIGMA[[i]], DELTA[[i]], DOF[i])
        logL <- sum(log(colSums(TAU)))
        sumTAU <- matrix(1, g, 1) %*% matrix(colSums(TAU),1,n)
        TAU <- TAU/sumTAU
    }
    return(list(TAU=TAU,logL=logL))
}

ETT.left <- function (a, mu, sigma, dof, integral = F, computeVar = T) {
    p <- nrow(sigma)
    Tp <- {
        if (integral) 
            pmt(mu - a, rep(0, p), sigma, dof)
        else 1
    }
    ET <- ett(sigma, dof, mu - a, integral, computeVar)
    mu <- matrix(mu, p)
    EX <- mu - ET$EX
    if (computeVar) 
        EXX <- ET$EXX - ET$EX %*% t(mu) - mu %*% t(ET$EX) + mu %*% 
            t(mu)
    return(if (computeVar) {
        list(EX = EX/Tp, EXX = EXX/Tp)
    } else EX/Tp)
}

ett <- function (sigma, dof, a = NULL, integral = F, computeVar = T) {
    p <- dim(sigma)[1]
    if (is.null(p)) {
        p <- 1
    }
    if (is.null(a)) 
        a <- rep(0, p)
    a <- as.numeric(a)
    if (p == 1) {
        Tp <- pt(a/sqrt(sigma), df = dof)
        EX <- -0.5 * sqrt(sigma * dof/pi) * (1 + a^2/dof/sigma)^(-0.5 * 
            dof + 0.5) * gamma(0.5 * dof - 0.5)/gamma(0.5 * dof)
        EXX <- sigma * dof/(dof - 2) * Tp + (dof - 1)/(dof - 
            2) * a * EX
        if (integral) {
            Tp <- 1
        }
        if (computeVar) 
            return(list(EX = EX/Tp, EXX = EXX/Tp))
        return(list(EX = EX/Tp))
    }
    else {
        if (integral) {
            Tp <- 1
        }
        else {
            Tp <- pmt(a, rep(0, p), sigma, dof)
        }
        xi <- ai <- rep(0, p)
        theta <- matrix(0, p, p)
        G1 <- gamma(0.5 * dof - 0.5)/gamma(0.5 * dof)
        V1 <- sqrt(0.5 * dof)
        V2 <- dof^(0.5 * dof)/(dof - 2)
        for (i in 1:p) {
            Sii <- sigma[i, i]
            S1 <- (dof + a[i]^2/Sii)/(dof - 1) * (sigma[-i, -i] - 
                sigma[-i, i] %*% t(sigma[-i, i])/Sii)
            A1 <- a[-i] - a[i]/Sii * sigma[-i, i]
            xi[i] <- (dof/(dof + a[i]^2/Sii))^(0.5 * dof - 0.5)/sqrt(2 * 
                pi * Sii) * G1 * V1 * pmt(A1, rep(0, p - 1), 
                S1, dof - 1)[1]
            if (computeVar) {
                if (p == 2) {
                  if (i == 1) {
                    S2 <- solve(sigma)
                    V3 <- dof + a %*% S2 %*% a
                    theta[1, 2] <- -V2/(V3^(0.5 * dof - 1))/(2 * 
                      pi * sqrt(Sii * sigma[2, 2] - sigma[1, 
                      2]^2))
                    theta[2, 1] <- theta[1, 2]
                  }
                }
                else {
                  if (i != p) {
                    for (j in (i + 1):p) {
                      S2 <- solve(matrix(c(Sii, sigma[i, j], 
                        sigma[i, j], sigma[j, j]), 2, 2))
                      V3 <- dof + a[c(i, j)] %*% S2 %*% a[c(i, 
                        j)]
                      S3 <- as.numeric(V3/(dof - 2)) * (sigma[-c(i, 
                        j), -c(i, j)] - matrix(sigma[-c(i, j), 
                        c(i, j)], p - 2, 2) %*% S2 %*% t(matrix(sigma[-c(i, 
                        j), c(i, j)], p - 2, 2)))
                      A2 <- a[-c(i, j)] - sigma[-c(i, j), c(i, 
                        j)] %*% S2 %*% a[c(i, j)]
                      theta[i, j] <- -V2/(V3^(0.5 * dof - 1))/(2 * 
                        pi * sqrt(Sii * sigma[j, j] - sigma[i, 
                        j]^2)) * pmt(as.numeric(A2), rep(0, p - 2), S3, dof - 2)
                      theta[j, i] <- theta[i, j]
                    }
                  }
                }
                theta[i, i] <- a[i] * xi[i]/Sii - sum(sigma[i, 
                  -i] * theta[i, -i])/Sii
            }
        }
    }
    EX <- -sigma %*% xi/Tp
    if (computeVar) {
        EXX <- dof/(dof - 2) * pmt(a, rep(0, p), dof/(dof - 2) * 
            sigma, dof - 2) * sigma - sigma %*% 
            theta %*% sigma
        EXX <- EXX/Tp
        return(list(EX = EX, EXX = EXX))
    }
    else {
        return(EX = EX)
    }
}


as.fmcfust <- function(g, mu, sigma, delta, dof=rep(1,g), pro=rep(1/g,g)) {
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

summary.fmcfust <- function(object, ...) {
    g <- length(object$dof)
    p <- nrow(object$sigma[[1]]) 
    if(is.null(p)) p <- 1
    if(g<1) stop("invalid fmcfust object: g < 1")
    cat("Finite Mixture of Multivarate CFUST distributions\n")
    if(g==1) cat("with ", g, " component\n\n")
    else cat("with ", g, " components\n\n")
    if(g==1) cat("Mean:\n")
    else cat("Component means:\n")
    means <- matrix(0,p,g) 
    for (i in 1:g) means[,i] <- matrix(object$mu[[i]],p,1)
    print(means)
    cat("\n\n")
    if(g==1) cat("Scale matrix:\n")
    else cat("Component scale matrices:\n")
    print(object$sigma)
    if(g==1) cat("\nSkewness matrix:\n")
    else cat("\nComponent skewness matrices:\n")
    print(object$delta)
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
    means <- matrix(0,p,g)
    for (i in 1:g) means[,i] <- matrix(x$mu[[i]],p,1)
    print(means)
    cat("\n\n")
    if(g==1) cat("Scale matrix:\n")
    else cat("Component scale matrices:\n")
    print(x$sigma)
    if(g==1) cat("\nSkewness matrix:\n")
    else cat("\nComponent skewness matrices:\n")
    print(x$delta)
    cat("\n\n")
    if(g==1) cat("Degrees of freedom:\n")
    else cat("Component degrees of freedom:\n")
    cat(x$dof, "\n\n")
    if(g>1) {
      cat("Component mixing proportions:\n")
      cat(x$pro, "\n\n")
    }
}

print.fmcfust <- function(x, ...) {
    g <- length(x$dof)
    cat("Finite Mixture of Multivarate CFUST Distribution\n")
    if(g==1) cat("with ", g, " component\n\n")
    else cat("with ", g, " components\n\n")
    obj <- as.fmcfust(g, x$mu, x$sigma, x$delta, x$dof, x$pro)
    obj$tau <- x$tau
    obj$clusters <- x$clusters
    obj$loglik <- x$loglik
    obj$aic <- x$aic
    obj$bic <- x$bic
    print(obj)
}
