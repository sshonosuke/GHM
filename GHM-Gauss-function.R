library(flexmix)

# function of GHM 
# delta: convergence criterion
# init: initial grouping number
# X should not include an intercept term
# G and L should be larger than 1

GHM.Gauss <- function(Y, X, ID, G=4, L=3){
  maxitr <- 500
  delta <- 0.01
  m <- length( unique(ID) )
  N <- length(Y)
  p <- dim(X)[2]
  
  # initial values 
  hmix <- matrix(0, m, L)
  hBeta <- array(0, c(m, p+1, L))
  hSig <- matrix(0.1, m, L)
  for(i in 1:m){
    init.fit <- flexmix(Y[ID==i]~X[ID==i,], k=L, model=FLXMRglm(family="gaussian"))
    dd <- dim(parameters(init.fit))[[2]]
    hBeta[i,,1:dd] <- parameters(init.fit)[1:(p+1),]
    hSig[i,1:dd] <- parameters(init.fit)[p+2,]
    hmix[i,1:dd] <- summary(init.fit)@comptab$ratio
  }
  
  num <- 10   # numner of initial points in k-means
  KM <- list()
  rate <- c()
  for(j in 1:num){
    KM[[j]] <- kmeans(hmix, G)
    rate[j] <- KM[[j]]$betweenss / KM[[j]]$totss
  }
  KM <- KM[[which.max(rate)]]
  hg <- KM$cluster
  hPi <- KM$centers
  hBeta <- apply(hBeta, c(2,3), median)
  hSig <- apply(hSig, 2, median)
  
  dimnames(hBeta)[[1]] <- paste0("v",0:p)
  dimnames(hBeta)[[2]] <- paste0("mix",1:L)
  dimnames(hPi)[[1]] <- paste0("g=", 1:G)
  dimnames(hPi)[[2]] <- paste0("mix", 1:L)
  
  # Iterations
  for(k in 1:maxitr){
    # save current values
    hBeta0 <- hBeta
    hSig0 <- hSig
    hPi0 <- hPi
    hg0 <- hg
    mu <- cbind(1, X)%*%hBeta 
    Sig.mat <- t( matrix(rep(hSig, N), L, N) )
    
    # E-step
    th <- exp(-700)  # thresholding value to avoid exact zero 
    dens <- dnorm(Y, mu, Sig.mat)
    weight <- matrix(NA, N, L)
    for(i in 1:m){
      pp <- hPi[hg[i], ] 
      dd <- dens[ID==i, ]
      wdd <- t( t(dd) * pp )
      sdd <- apply(wdd, 1, sum)
      sdd[ sdd<th ] <- th
      weight[ID==i, ] <- wdd/sdd
    }
    
    # M-step (Beta and Sigma)
    for(k in 1:L){
      ww <- weight[,k]
      if( sum(ww)>5 ){
        fit <- lm(Y~X, weights=ww)
        hBeta[,k] <- coef(fit)
        hSig[k] <- sqrt( sum(ww*fit$residuals^2)/sum(ww) )
      }
    }
    
    # M-step (mixing proportion)
    Ind <- hg[ID]
    for(g in 1:G){
      ww <- weight[Ind==g,]
      if( dim(ww)[1]>0 ){ hPi[g,] <- apply(ww,2,sum)/dim(ww)[1] }
    }
    
    # M-step (grouping variable)
    val <- matrix(NA,m,G)
    for(i in 1:m){
      for(g in 1:G){
        val[i,g] <- sum( log(hPi[g,]) * apply(weight[ID==i,], 2, sum) )
      }
    }
    hg <- apply(val, 1, which.max)
    
    # convergence check
    conv <- sum( abs(hBeta-hBeta0) )
    if( conv<delta ){  break  }
  }
  
  # Maximum likelihood value
  hmu <- cbind(1, X)%*%hBeta 
  Sig.mat <- t( matrix(rep(hSig, N), L, N) )
  dens <- dnorm(Y, hmu, Sig.mat)
  ML <- matrix(NA, N, L)
  for(i in 1:m){
    pp <- hPi[hg[i],] 
    dd <- dens[ID==i,]
    ML[ID==i,] <- t(t(dd)*pp)
  }
  ML <- sum( log( apply(ML, 1, sum) ) )
  BIC <- -2*ML + log(N)*( (p+1)*(L-1) + G*(L-1) + m )
  Result <- list(ML=ML, Beta=hBeta, Sig=hSig, Pi=hPi, g=hg, BIC=BIC)
  return(Result)
}





###     adaptive GHM     ###
aGHM.Gauss <- function(Y, X, ID, rG=c(2,6), rL=c(2,5), print=F){
  # Candidates of models
  sG <- rG[1]:rG[2]
  lG <- length(sG)
  sL <- rL[1]:rL[2]
  lL <- length(sL)
  DD=list()
  
  # Calculation of BIC
  BIC <- matrix(NA,lG,lL)
  dimnames(BIC)[[1]] <- paste0("G=", sG)
  dimnames(BIC)[[2]] <- paste0("L=", sL)
  for(g in 1:lG){
    dd <- list()
    for(l in 1:lL){
      dd[[l]] <- GHM.Gauss(Y, X, ID, G=sG[g], L=sL[l])
      BIC[g, l] <- dd[[l]]$BIC
      if(print){ 
        print( paste0("G=",sG[g],", L=",sL[l]) )
        print( BIC[g,l] ) 
      }
    }
    DD[[g]] <- dd
  }
  
  # optimal selection
  opt <- which.min(BIC)
  opt.G <- opt - lG*floor(opt/lG)
  if(opt.G==0){ opt.G <- lG }
  opt.L <- ceiling(opt/lG)
  fit.best <- DD[[opt.G]][[opt.L]]
  
  # Summary
  Res <- list(fit.best=fit.best, fit.all=DD, BIC=BIC)
  return(Res)
}




