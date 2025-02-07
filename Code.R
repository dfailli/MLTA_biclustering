# N sending nodes
# R receiving nodes
# J covariates
# G components (units clusters)
# D segments (variables clusters)

# beta (multinomial logit coefficients) Jx(G-1)
# eta (prior) NxG
# z (posterior) NxG
# a (column membership matrix) DxMxG
# b (logit intercept) 1xG
# mu.d (latent trait mean) Dx1
# u (latent trait) DxN

library(Rcpp)
library(devtools)
library(processx)


## Functions ####

Rcpp::sourceCpp("logistic.cpp") # logistic for b
Rcpp::sourceCpp("logistic2.cpp") # logistic for mu

randgenu=function(g, k, d)
{
  #randomly generate a membership matrix U
  U=matrix(0,k,d)
  U[,1]=1;
  U[1:d,]=diag(diag(matrix(1,d,d)))
  for(i in (d+1):k)
  {
    U[i,]=U[i, sample(d)];
    
  }
  return(U)
}


f_mlta <- function(X, DM, G, D, tol, maxiter, pdGH, beta0)
{
  
  
  N <- nrow(X)          
  M <- ncol(X)           
  J <- ncol(DM)       
  
  
  # Gauss-Hermite quadrature
  
  npoints <- 7           
  Q <- npoints ^ D       
  GaHer <- glmmML::ghq(npoints, FALSE)                     
  ugh <- as.matrix(expand.grid(rep(list(GaHer$zeros), D))) 
  ugh.star <- sqrt(2)*ugh                                  
  p.gh = as.matrix(expand.grid(rep(list(GaHer$weights), times = D)))	
  W =  (2)^(D/2) * exp(apply(ugh, 1, crossprod)) * apply(p.gh, 1, prod)
  W.prod =  apply(p.gh, 1, prod)
  Phi <- apply(ugh.star, 1, mvtnorm::dmvnorm)               
  
  
  # Matrix expansion
  
  onesnJ=matrix(1,1,N*M) 
  onesj=matrix(1,1,M)
  onesn=matrix(1,1,N)
  onesq=matrix(1,1,Q)   
  onesnJq=matrix(1,1,N*M*Q)
  onesnq=matrix(1,1,N*Q)
  
  
  # Incidence matrix expansion (QGMNx1)
  ym=rep(rep(c(X),G),Q) 
  
  
  # EM initialization
  
  # beta
  if(!is.null(beta0)) beta = beta0 else{
    beta <- rep(0,J*(G-1))
    beta <- matrix(beta,ncol=G-1)
  }
  
  
  # eta (priors)
  exb <- exp(DM %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  
  
  # z (posteriors)
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N)   z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,])) 
  
  
  # a_gk (kmeans)
  a=array(0,c(D,M,G))
  for (g in 1:G){
    X[z[,g]==1,]
    result.a = kmeans(t(X[z[,g]==1,]), D, nstart=50, iter.max=1000) # kmeans
    a.temp=matrix(0,D, ncol=M)
    for (k in 1:M){
      a.temp[result.a$cluster[k],k]=1
    }
    a[,,g]=a.temp
  }
  
  
  # b_g, mu_d
  
  b <- c(rnorm(G))
  mu.d=c(rnorm(D))
  
  
  # Iterative process
  
  v <- matrix(0, N, G)
  fy.uz <- array(0, c(N, Q, G))
  fy.u=matrix(0,N,Q)
  obj=array(0,c(N,Q,G))
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  tol <- 10^-4
  print(c(beta))
  
  
  while (diff > tol & iter < maxiter)
  {
    
    iter <- iter + 1
    beta.old=beta
    ll_old <- ll
    
    
    # E-step
    
    if(D==1){
      for (g in 1:G)
      {
        Agh <- t(tcrossprod(matrix(a[, , g],M,D), t(t(ugh.star)+mu.d)) + b[g]) 
        pgh <- 1 / (1 + exp(-Agh)) 
        fy.uz[, , g] <-
          exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
        # fy.uz[is.nan(fy.uz[, , g])] <- 0
      }
    }else{
      for (g in 1:G)
      {
        Agh <- t(tcrossprod(t(a[, , g]), t(t(ugh.star)+mu.d)) + b[g]) 
        pgh <- 1 / (1 + exp(-Agh)) 
        fy.uz[, , g] <-
          exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
        # fy.uz[is.nan(fy.uz[, , g])] <- 0
      }
    }
    
    for(q in 1:Q){
      fy.u[,q]=rowSums(eta*(fy.uz[,q,])) 
    }
    fyu = (t(t(fy.u)*W*Phi))
    
    fy.z <- apply(aperm(apply(fy.uz,c(1,3),function(x) x*W*Phi),c(2,1,3)),c(1,3),sum)
    fyz <- eta * (fy.z)
    fy <- apply(fyz, 1, sum)
    z <- fyz / fy

    
    for(q in 1:Q){
      obj[,q,]=eta*(fy.uz[,q,]) 
    }
    fyuz=aperm(apply(obj,c(1,3),function(x) x*W*Phi),c(2,1,3)) 
    fuz.y=apply(fyuz,c(2,3),function(x) x/fy) # NxQxG
    obj3.new=apply(fuz.y,c(2,3),function(x) t(onesj)%x%x) # MNxQxG
    ww2=c(aperm(obj3.new,c(1,3,2)) ) # QGMN
    
    
    # M- step
    
    U=aperm(a,c(2,3,1))   # GxMxD
    U=matrix(U,M*G,D)     # GMxD
    Um=U%x%t(onesn)       # GMNxD
    Umm=t(onesq)%x%Um     # QGMNxD
    zz=t(onesq)%x%(diag(rep(1,G))%x%t(onesnJ)) # QGMNxG
    
    yy <- matrix(ym, Q*G*M*N, 1)
    xxb=zz
    xxm=Umm
    ww22 <- matrix(ww2, Q*G*M*N, 1)
    offsb=Umm %*% mu.d + matrix(c(Um %*% t(ugh.star)),ncol=1)
    offsm=(zz%*%b + matrix(c(Um %*% t(ugh.star)),ncol=1))
    coeffb=b
    coeffm=mu.d
    
    
    # Update b_g and mu_d

    bcoeff <- callLogisticRegressionb(yy, xxb, offsb, ww22, coeffb)
    mucoeff <- callLogisticRegressionm(yy, xxm, offsm, ww22, coeffm)
    
    
    b=bcoeff[1:G]
    mu.d=mucoeff[1:D]
    
    
    # Update beta
    
    lk = sum(z*log(eta))
    it = 0; lko = lk
    XXdis = array(0,c(G,(G-1)*ncol(DM),N))
    for(i in 1:N){
      XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(DM[i,]))
    }
    while((lk-lko>10^-6 & it<100) | it==0){
      it = it+1; lko = lk 
      sc = 0; Fi = 0
      for(i in 1:N){
        pdis = eta[i,]
        sc = sc+t(XXdis[,,i])%*%(z[i,]-pdis)
        Fi = Fi+t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
      }
      
      dbe = as.vector(ginv(Fi)%*%sc)
      mdbe = max(abs(dbe))
      if(mdbe>0.5) dbe = dbe/mdbe*0.5
      be0 = c(beta)
      flag = TRUE
      while(flag){
        beta = be0+dbe
        Eta = matrix(0,N,G)
        for(i in 1:N){
          
          if(ncol(DM)==1) Eta[i,] = XXdis[,,i]*beta
          else Eta[i,] = XXdis[,,i]%*%beta
        }	
        if(max(abs(Eta))>100){
          dbe = dbe/2
          flag = TRUE	
        }else{
          flag = FALSE
        }	        	
      }
      if(iter/10 == floor(iter/10))       print(beta)
      
      beta = matrix(beta, J, G-1)
      exb <- exp(DM %*% beta) # update priors
      eta <- cbind(1,exb)/(rowSums(exb)+1)
      
      lk = sum(z*log(eta))
    }
    
    
    # Update a_gk 
    
    if(D==1){
      for (g in 1:G)
      {
        Agh <- t(tcrossprod(matrix(a[, , g],M,D), t(t(ugh.star)+mu.d)) + b[g]) 
        pgh <- 1 / (1 + exp(-Agh)) 
        fy.uz[, , g] <-
          exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
        # fy.uz[is.nan(fy.uz[, , g])] <- 0
      }
    }else{
      for (g in 1:G)
      {
        Agh <- t(tcrossprod(t(a[, , g]), t(t(ugh.star)+mu.d)) + b[g]) 
        pgh <- 1 / (1 + exp(-Agh)) 
        fy.uz[, , g] <-
          exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
        # fy.uz[is.nan(fy.uz[, , g])] <- 0
      }
    }
    
    for(q in 1:Q){
      obj[,q,]=eta*(fy.uz[,q,]) 
    }
    fyuz=aperm(apply(obj,c(1,3),function(x) x*W*Phi),c(2,1,3)) 
    fuz.y=apply(fyuz,c(2,3),function(x) x/fy) # NxQxG
    
    IDj=1:M                           # unique index for receiving node
    IDj=rep(IDj, each=N)              # repeated index of receiving, each N times   NMx1
    objective=array(NA, c(M, G, D))   # objective function MxG D times
    
    for (g in 1:G){
      
      for(d in 1:D){
        
        lp = b[g] +  mu.d[d] + ugh.star[,d]                     # linear predictor
        pgh=1 / (1 + exp(-lp))                                  # prob
        fy2=sapply(pgh,function(x) dbinom(X,1,x, log=T))        # bernoulli density (log) MxN
        
        pesi=t(onesj)%x%fuz.y[,,g]                              # weight
        
        likelihard=rowSums(fy2*(pesi))                          # sum over Q
        likejkdhard=tapply(likelihard, IDj, sum)                # Mx1 (sum over N)
        objective[,g,d]=likejkdhard                             # log-likelihood for every G, K, D  # MxG D times
      }
    }
    
    maxobjective=apply(objective, c(1,2), max)
    Maxobjective=array(maxobjective, c(M, G, D))
    Uest=objective==Maxobjective
    Uest=Uest+0
    

    for(g in 1:G){
      niter=0
      while(min(apply(matrix(Uest[,g,],M,D), c(2), sum))==0)
      {
        niter=niter+1
        idxf=which(apply(Uest[,g,], c(2), sum)==0)   
        for (h in 1:length(idxf))
        {
          idxr=order(objective[,  g, idxf[h]], decreasing=T)  
          Uest[idxr[niter],g, idxf[h]]=1
          Uest[idxr[niter],g,-idxf[h]]=0
        }
      }
    }
    
    a=aperm(Uest, c(3,1,2))
    
    
    # Log-likelihood
    
    if(D==1){
      for (g in 1:G)
      {
        Agh <- t(tcrossprod(matrix(a[, , g],M,D), t(t(ugh.star)+mu.d)) + b[g]) 
        pgh <- 1 / (1 + exp(-Agh)) 
        fy.uz[, , g] <-
          exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
        # fy.uz[is.nan(fy.uz[, , g])] <- 0
      }
    }else{
      for (g in 1:G)
      {
        Agh <- t(tcrossprod(t(a[, , g]), t(t(ugh.star)+mu.d)) + b[g]) 
        pgh <- 1 / (1 + exp(-Agh)) 
        fy.uz[, , g] <-
          exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
        # fy.uz[is.nan(fy.uz[, , g])] <- 0
      }
    }
    fy.z <- apply(aperm(apply(fy.uz,c(1,3),function(x) x*W*Phi),c(2,1,3)),c(1,3),sum)
    fyz <- eta * (fy.z)
    fy <- apply(fyz, 1, sum)
    ll <- sum(log(fy))
    
    
    # Stopping Criterion
    
    diff <- sum(abs(ll - ll_old))
    
    if(sum((ll - ll_old))<0) print(paste(iter, ll, sum((ll - ll_old)), "aaaaaa"))
    else print(paste(iter, ll, sum((ll - ll_old))))
    
  } # end while(diff>tol)
  
  LL <- ll
  BIC <-
    -2 * LL + (G * (M * (D) - D * (D - 1) / 2) + J*(G - 1) + G + D) * log(N)
  
  rownames(a) <- rownames(a, do.NULL = FALSE, prefix = "D=")
  colnames(a) <- colnames(a, do.NULL = FALSE, prefix = "Item ")
  dimnames(a)[[3]] <- paste("Group ", seq(1, G) , sep = '')
  
  colnames(beta) <- 2:G
  rownames(beta) <- c(colnames(DM))
  names(LL) <- c("Log-Likelihood:")
  names(BIC) <- c("BIC:")
  
  ord = order(b)
  be = cbind(0, beta)
  be = be[,ord]
  beta = be[,2:G]-be[,1]
  
  b = b[ord]
  bb = b-b[2]
  b = bb
  
  ord.mu = order(mu.d)
  mu.d = mu.d[ord.mu]
  mu.d=mu.d-mu.d[2]
  
  log_likelihood <- function(bb) {
    maxb=G
    maxmu=maxb+D
    maxbeta=maxmu+(J*(G-1))
    fy.uz=array(0,c(N,Q,G))
    if(D==1){
      for (g in 1:G)
      {
        Agh <- t(tcrossprod(matrix(a[, , g],M,D), t(t(ugh.star)+bb[(maxb+1):(maxmu)])) + bb[1:maxb][g])
        pgh <- 1 / (1 + exp(-Agh))
        fy.uz[, , g] <-
          exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh)))
        # fy.uz[is.nan(fy.uz[, , g])] <- 0
      }
    }else{
      for (g in 1:G)
      {
        Agh <- t(tcrossprod(t(a[, , g]), t(t(ugh.star)+bb[(maxb+1):(maxmu)])) + bb[1:maxb][g])
        pgh <- 1 / (1 + exp(-Agh))
        fy.uz[, , g] <-
          exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh)))
        # fy.uz[is.nan(fy.uz[, , g])] <- 0
      }
    }
    fy.z <- apply(aperm(apply(fy.uz,c(1,3),function(x) x*W*Phi),c(2,1,3)),c(1,3),sum)
    
    betaaa=matrix(bb[(maxmu+1):maxbeta],J,G-1)
    exb <- exp(DM %*% betaaa)
    eta <- cbind(1,exb)/(rowSums(exb)+1)
    
    fyz <- eta * (fy.z)
    fy <- apply(fyz, 1, sum)
    loglik <- sum(log(fy))
    
    return(loglik)
  }
  
  result <- optim(c(b,mu.d,c(beta)), log_likelihood, control = list(fnscale = -1, maxit=100), hessian = T)
  score = numDeriv::grad(log_likelihood, x = result$par)
  hessiana = result$hessian
  se=sqrt(diag(MASS::ginv(hessiana) %*% score %*% t(score) %*% MASS::ginv(hessiana)))
  
  out <- list(
    a = a,
    b = b,
    mu.d = mu.d,
    eta = eta,
    se=se,
    z = z,
    LL = LL,
    BIC = BIC,
    beta = beta, 
    rrr=exp(beta)
  )
  
  out
}


f_mlta_nstarts <- function(X, DM, G, D, nstarts, tol, maxiter, pdGH, beta0)
{
  out <- try(f_mlta(X, DM, G, D, tol, maxiter, pdGH, beta0))
  
  if(nstarts > 1){
    for(i in 2:nstarts) {
      
      out1 <- try(f_mlta(X, DM, G, D, tol, maxiter, pdGH, beta0))
      if(out1$LL > out$LL) out <- out1
      
    }
  }
  return(out)
}


ResTable <- function(bicG, restype)
{
  if (restype == 'll') {
    resBIC <- vector('list', 1)
    names(resBIC) <- 'Table of LL'
    resBIC[[1]] <- bicG
  }
  
  if (restype == 'llva') {
    resBIC <- vector('list', 1)
    names(resBIC) <- 'Table of LL'
    resBIC[[1]] <- bicG
  }
  
  if (restype == 'BIC') {
    resBIC <- vector('list', 2)
    names(resBIC) <- c('Table of BIC Results', 'Model Selection')
    resBIC[[1]] <- bicG
    resBIC[[2]] <-
      paste(colnames(bicG)[t(bicG == min(bicG)) * seq(1:ncol(bicG))], rownames(bicG)[(bicG ==
                                                                                        min(bicG)) * seq(1:nrow(bicG))], sep = ', ')
    names(resBIC[[2]]) <- 'Model with lower BIC:'
  }
  
  resBIC
}


tableBIC <- function(out){
  lout <- length(out) - 1
  bicG <- numeric(lout)
  names(bicG) <- names(out[- (lout + 1)])
  
  for(i in 1:lout) bicG[i] <- out[[i]]$BIC
  
  resBIC <- vector('list', 2)
  names(resBIC) <- c('Model Selection', 'Table of BIC Results')
  resBIC[[1]] <- names(bicG)[bicG == min(bicG)]
  names(resBIC[[1]]) <- 'Model with lower BIC:'
  resBIC[[2]] <- bicG
  resBIC
}



f_mlta_methods <- function(X, DM, G, D, nstarts, tol, maxiter, pdGH, beta0=NULL)
{
  if (D > 0 && G > 1) {
    out <- f_mlta_nstarts(X, DM, G, D, nstarts, tol, maxiter, pdGH, beta0=NULL)
  }
  class(out) <- c("mlta")
  return(out)
}

f_lca <-function(X, DM, G, tol, maxiter, beta0){ 
  
  N <- nrow(X)
  M <- ncol(X) 
  J <- ncol(DM)
  
  # Initialize EM Algorithm
  
  if(!is.null(beta0)) beta = beta0 else{
    beta <- rep(0,J*(G-1))
    beta <- matrix(beta,ncol=G-1)
  }
  
  exb <- exp(DM %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  # priors
  
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N)   z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,]))
  
  # Set up
  
  p <- (t(z) %*% X) / colSums(z) # proportion of success in group g
  
  # Initialize Aitken acceleration
  
  l <- numeric(3)
  lA <- numeric(2)
  
  # Iterative process
  
  v <- matrix(0, N, G)
  W <- list()
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  tol <- 0.1 ^ 6
  print(c(beta))
  
  while(diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    beta.old=beta
    ll_old <- ll
    #diff.ll <- 0
    
    
    # M-step 
    
    lk = sum(z*log(eta))
    it = 0; lko = lk
    XXdis = array(0,c(G,(G-1)*ncol(DM),N))
    for(i in 1:N){
      XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(DM[i,]))
    }
    while((lk-lko>10^-6 & it<100) | it==0){
      it = it+1; lko = lk 
      sc = 0; Fi = 0
      for(i in 1:N){
        pdis = eta[i,]
        sc = sc+t(XXdis[,,i])%*%(z[i,]-pdis)
        Fi = Fi+t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
      }
      
      dbe = as.vector(ginv(Fi)%*%sc)
      mdbe = max(abs(dbe))
      if(mdbe>0.5) dbe = dbe/mdbe*0.5
      be0 = c(beta)
      flag = TRUE
      while(flag){
        beta = be0+dbe
        Eta = matrix(0,N,G)
        for(i in 1:N){
          
          if(ncol(DM)==1) Eta[i,] = XXdis[,,i]*beta
          else Eta[i,] = XXdis[,,i]%*%beta
        }	
        if(max(abs(Eta))>100){
          dbe = dbe/2
          flag = TRUE	
        }else{
          flag = FALSE
        }	        	
      }
      if(iter/10 == floor(iter/10))       print(beta)
      
      beta = matrix(beta, J, G-1)    
      exb <- exp(DM %*% beta) # updfe priors
      eta <- cbind(1,exb)/(rowSums(exb)+1)
      
      lk = sum(z*log(eta))
    }
    
    
    # E-step
    PS = array(apply(p, 1, function(xx) dbinom(c(X), 1, xx)), c(N, M, G))
    v <- eta *  apply(PS, c(1,3), prod)
    # v[is.nan(v)] <- 0
    # v[v < 0] <- 0
    vsum <- rowSums(v)
    z <- v / vsum      # posterior
    ll <- rep(1, N) %*% log(vsum) # log-likelihood
    
    
    diff <- sum(abs(ll - ll_old))
    
    if(sum((ll - ll_old))<0) print(paste(iter, ll, sum((ll - ll_old)), "aaaaaa"))
    else print(paste(iter, ll, sum((ll - ll_old))))
    
  } # end EM
  
  p <- (t(z) %*% X) / colSums(z)	
  
  BIC <- -2*ll + (J*(G-1)+G*M) * log(N)
  
  expected <- vsum * N
  
  
  colnames(p) <- NULL
  rownames(p) <- rownames(p, do.NULL = FALSE, prefix = "Group ")
  colnames(p) <- colnames(p, do.NULL = FALSE, prefix = "p_g")
  colnames(beta) <- 2:G
  rownames(beta) <- c(colnames(DM))
  names(ll)<- "Log-Likelihood:"
  names(BIC)<-"BIC:"
  
  out <- list(p = p, eta = eta, LL = ll, BIC = BIC)
  list(p = p, eta = eta, LL = ll, BIC = BIC, z = z, expected = expected,
       beta = beta) 
}



f_lca_nstarts <- function(X, DM, G, nstarts, tol, maxiter, beta0) {
  
  out <- f_lca(X, DM, G, tol, maxiter, beta0)
  
  foreach(i=2:nstarts, .packages = c("MASS","igraph"), .export = "f_lca") %dopar% {
    out1 <- f_lca(X, DM, G, tol, maxiter, beta0)
    if(out1$LL > out$LL) out <- out1
  }
  out
}



lca <- function(X, DM, G, nstarts = 3, tol = 0.1^2, maxiter = 250, beta0=NULL) {
  
  if (any(G < 1)) {
    print("Specify G > 0!")
    return("Specify G > 0!")
  }
  
  if (any(G == 1)) {
    out <- f_lca_nstarts(X, DM, G, nstarts, tol, maxiter, beta0 = beta0)
  } else{
    if (length(G) == 1) {
      out <- f_lca_nstarts(X, DM, G, nstarts, tol, maxiter, beta0 = beta0)
    } else{
      out <- vector("list", length(G) + 1)
      names(out) <- c(paste('G', G, sep = '='), 'BIC')
      i <- 0
      for (g in G) {
        i <- i + 1
        out[[i]] <- f_lca_nstarts(X, DM, g, nstarts, tol, maxiter, beta0 = beta0)
      }
      out[[length(G) + 1]] <- tableBIC(out)
    }
  }
  out
}


randgenug1=function(k, d)
{
  #randomly generate a membership matrix U
  U=matrix(0,k,d)
  U[,1]=1;
  U[1:d,]=diag(diag(matrix(1,d,d)))
  for(i in (d+1):k)
  {
    U[i,]=U[i, sample(d)];
    
  }
  return(U)
}


f_lta <- function(X, D, tol, maxiter, pdGH) {
  
  
  N <- nrow(X)          
  M <- ncol(X)           
  
  
  # Gauss-Hermite quadrature
  
  npoints <- 7           
  Q <- npoints ^ D       
  GaHer <- glmmML::ghq(npoints, FALSE)                     
  ugh <- as.matrix(expand.grid(rep(list(GaHer$zeros), D))) 
  ugh.star <- sqrt(2)*ugh                                  
  p.gh = as.matrix(expand.grid(rep(list(GaHer$weights), times = D)))	
  W =  (2)^(D/2) * exp(apply(ugh, 1, crossprod)) * apply(p.gh, 1, prod)
  W.prod =  apply(p.gh, 1, prod)
  Phi <- apply(ugh.star, 1, mvtnorm::dmvnorm)                
  
  onesnJ=matrix(1,1,N*M) 
  onesj=matrix(1,1,M)
  onesn=matrix(1,1,N)
  onesq=matrix(1,1,Q)   
  onesnJq=matrix(1,1,N*M*Q)
  onesnq=matrix(1,1,N*Q)
  
  ym=rep(c(X),Q) 
  
  
  # a_gk (kmeans)
  a=array(0,c(D,M))
  result.a = kmeans(t(X), D, nstart=50, iter.max=1000) # kmeans
  a.temp=matrix(0,D, ncol=M)
  for (k in 1:M){
    a.temp[result.a$cluster[k],k]=1
  }
  a=a.temp
  
  
  # b_g, mu_d
  
  b <- c(rnorm(1))
  mu.d=c(rnorm(D))
  
  
  # Iterative process
  
  fy.u=matrix(0,N,Q)
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  tol <- 10^-4
  
  
  while (diff > tol & iter < maxiter)
  {
    
    iter <- iter + 1
    ll_old <- ll
    
    
    # E-step
    
    Agh <- t(tcrossprod(t(a), t(t(ugh.star)+mu.d)) + b) 
    pgh <- 1 / (1 + exp(-Agh)) 
    fy.u <- exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
    # fy.uz[is.nan(fy.uz[, , g])] <- 0
    
    fyu = (t(t(fy.u)*W*Phi))
    
    fy <- apply(fyu, 1, sum)
    
    fu.y=apply(fyu,c(2),function(x) x/fy) # NxQ
    prova3.new=apply(fu.y,c(2),function(x) t(onesj)%x%x) # MNxQ
    ww2=c(aperm(prova3.new,c(1,2)) ) # QMN
    
    
    # M- step
    
    U=t(a)   # MxD
    U=matrix(U,M,D)     # MxD
    Um=U%x%t(onesn)       # MNxD
    Umm=t(onesq)%x%Um     # QMNxD
    zz=t(onesq)%x%(1%x%t(onesnJ)) # QMNx1
    
    # Update b_g and mu_d
    
    mod2=glm(cbind(ym,1-ym) ~ zz - 1, offset=Umm %*% mu.d + matrix(c(Um %*% t(ugh.star)),ncol=1),
             family = binomial(link = logit), weights=ww2) #control = list(maxit = 2)
    mod1=glm(cbind(ym,1-ym) ~ Umm  - 1, offset=(zz%*%b + matrix(c(Um %*% t(ugh.star)),ncol=1)),
             family = binomial(link = logit), weights=ww2)#control = list(maxit = 2)
    
    b=mod2$coefficients
    mu.d=mod1$coefficients
    
    
    # Update a_gk 
    
    Agh <- t(tcrossprod(t(a), t(t(ugh.star)+mu.d)) + b) 
    pgh <- 1 / (1 + exp(-Agh)) 
    fy.u <- exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
    
    fyu = (t(t(fy.u)*W*Phi))
    
    fy <- apply(fyu, 1, sum)
    
    fu.y=apply(fyu,c(2),function(x) x/fy) # NxQ
    
    IDj=1:M                           # unique index for receiving node
    IDj=rep(IDj, each=N)              # repeated index of receiving, each N times   NMx1
    objective=array(NA, c(M, D))      # objective function MxD 
    
    for(d in 1:D){
      lp = b +  mu.d[d] + ugh.star[,d]                         # linear predictor
      pgh=1 / (1 + exp(-lp))                                  # prob
      fy2=sapply(pgh,function(x) dbinom(X,1,x, log=T))        # bernoulli density (log) MxN
      
      pesi=t(onesj)%x%fu.y                                    # weight
      
      likelihard=rowSums(fy2*(pesi))                          # sum over Q
      likejkdhard=tapply(likelihard, IDj, sum)                # Mx1 (sum over N)
      objective[,d]=likejkdhard                               # log-likelihood for every K, D  # MxD 
    }
    
    maxobjective=apply(objective, c(1), max)
    Maxobjective=array(maxobjective, c(M, D))
    Uest=objective==Maxobjective
    Uest=Uest+0
    
    
    niter=0
    while(min(apply(Uest, c(2), sum))==0)
    {
      niter=niter+1
      idxf=which(apply(Uest, c(2), sum)==0)
      
      for (h in 1:length(idxf))
      {
        idxr=order(objective[,  idxf[h]], decreasing=T)
        Uest[idxr[niter], idxf[h]]=1
        Uest[idxr[niter],-idxf[h]]=0
      }
    }
    a=t(Uest)
    
    
    # Log-likelihood
    
    Agh <- t(tcrossprod(t(a), t(t(ugh.star)+mu.d)) + b) 
    pgh <- 1 / (1 + exp(-Agh)) 
    fy.u <- exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh))) 
    
    fyu = (t(t(fy.u)*W*Phi))
    
    fy <- apply(fyu, 1, sum)
    
    ll <- sum(log(fy))
    
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
    if(sum((ll - ll_old))<0) print(paste(iter, ll, sum((ll - ll_old)), "aaaaaa"))
    else print(paste(iter, ll, sum((ll - ll_old))))
    
  } # end while(diff>tol)
  
  LL <- ll
  BIC <-
    -2 * LL + ((M * (D) - D * (D - 1) / 2) + 1 + D) * log(N)
  
  names(LL) <- c("Log-Likelihood:")
  names(BIC) <- c("BIC:")
  
  out <- list(
    a = a,
    b = b,
    mu.d = mu.d,
    LL = LL,
    BIC = BIC
  )
  
  out
}


f_lta_nstarts <- function(X, D, nstarts, tol, maxiter, pdGH)
{
  out <- f_lta(X, D, tol, maxiter, pdGH)
  if(nstarts > 1){ 
    for(i in 2:nstarts){
      out1 <- f_lta(X, D, tol, maxiter, pdGH)
      if(out1$LL > out$LL) out <- out1
    }
  }
  out
}


lta <- function(X, D, nstarts = 3, tol = 0.1^2, maxiter = 250, pdGH = 21) {
  
  if(any(D == 0)) stop("D must be > 0")
  
  if(length(D) == 1){ 
    out <- f_lta_nstarts(X, D, nstarts, tol, maxiter, pdGH)
  }else{
    out<-vector("list", length(D) + 1)
    names(out) <- c(paste('Dim y', D, sep = '='), 'BIC')
    i<-0
    for(diy in D){
      i <- i + 1
      out[[i]] <- f_lta_nstarts(X, diy, nstarts, tol, maxiter, pdGH)
    }
    out[[length(D) + 1]]<-tableBIC(out)
  }
  out
}


mlta <- function(X, DM, G, D, nstarts = 3, tol = 10^-5, maxiter = 250, pdGH = 50, 
                 beta0 = NULL){
  if (any(G < 1)) {
    print("Specify G > 0!")
    return("Specify G > 0!")
  }
  
  if (any(D < 0)) {
    print("Specify D >= 0!")
    return("Specify D >= 0!")
  }
  
  if (length(D) == 1 && length(G) == 1) {
    out <- f_mlta_methods(X, DM, G, D, nstarts, tol, maxiter, pdGH, beta0 = beta0)
  } else{
    out <- vector("list", length(D) * length(G) + 3)
    names(out) <- c(t(outer(
      paste('G=', G, sep = ''),
      paste('dim y=', D, sep = ''),
      paste,
      sep = ','
    )),
    'BIC', 'LL', 'LLva')
    
    bictab <- matrix(0, length(D), length(G))
    lltab <- matrix(0, length(D), length(G))
    
    rownames(bictab) <- paste('D=', D, sep = '')
    colnames(bictab) <- paste('G=', G, sep = '')
    
    rownames(lltab) <- paste('D=', D, sep = '')
    colnames(lltab) <- paste('G=', G, sep = '')
    
    llvatab <- lltab
    
    i <- 0
    for (g in G){
      for (diy in D) {
        i <- i + 1
        
        out[[i]] <- f_mlta_methods(X, DM, g, diy, nstarts, tol, maxiter, pdGH, beta0=beta0)
        print("ciao")
        bictab[i] <- out[[i]]$BIC
        lltab[i] <- out[[i]]$LL
        
        if (diy == 0) {
          llvatab[i] <- out[[i]]$LL
        }
        # } else{
        #   llvatab[i] <- out[[i]]$LLva
        # }
        
        out[[length(G) * length(D) + 1]] <- ResTable(bictab, restype = 'BIC')
        out[[length(G) * length(D) + 2]] <- ResTable(lltab, restype = 'll')
        out[[length(G) * length(D) + 3]] <- ResTable(llvatab, restype = 'llva')
      }
    }
    class(out) <- "mmlta"
  }
  
  out
}



print.mlta <- function(x){
  stopifnot(inherits(x, 'mlta'))
  cat("b:\n")
  print(x$b) 
  cat("\nw:\n")  
  print(x$a)
  cat("\neta:\n")
  print(x$eta)
  cat("\n")
  print(x$LL)
  print(x$BIC) 
}



print.mmlta <- function(x){
  cat("Log-Likelihood:\n")
  print(x$LL$`Table of LL (G-H Quadrature correction)`)
  cat("BIC:\n")
  print(x$BIC$`Table of BIC Results`)
  print(x$BIC$`Model Selection`)
}


## Simulation study ####

library(MASS)
library(mclust)
library(fastglm)
library(foreach)
library(parallel)
library(doParallel)
registerDoParallel(cores=50)
getDoParWorkers()

a.sim=list()
b.sim <- list()
mu.sim=list()
beta.sim <- list()

rand.g <- c()
rand.d <- c()

S <- 100
nstarts <- 100
G <- 3
D <- 2
N <- 100
M <- 20


# G=3 D=2
beta = matrix(c(1,-0.4,1.5,-0.9),2,2)
b0=c(-1.7, 0, 1.7)
mu0=c(-2,0.5)

bic.sim <- array(,c(3,3,S)) # (for BIC simulation)
icl.sim <- array(,c(3,3,S)) # (for ICL simulation)

# # G=4 D=3
# beta = matrix(c(1,-0.4, 1.5,-0.9, 2, -1.3), 2, 3, byrow=F)
# b0=c(-1.7, 0, 1.7, 0.7)
# mu0=c(-2,0.5,1.5)

res <- foreach (s=1:S, .packages = c("MASS","igraph")) %dopar% {
  
  set.seed(s)
  
  DM = cbind(rep(1, N), rnorm(N, 1, 1))
  exb <- exp(DM %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N){
    z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,]))
  }
  
  ord.g <- order(b0)
  cl.g = apply(z[,ord.g], 1, which.max) 
  
  Ulist=lapply(1:G,randgenu,M,D)
  Utrue=array(unlist(Ulist),c(M,D,G))    
  a=array(0,c(D,M,G))
  a=aperm(Utrue,c(2,1,3))
  
  ord.d=order(mu0)
  if(D==1){
    a.true=a
    a0=a.true
  }else{
    a.true=a[ord.d, , ]
    a.true=a0=a.true[ , ,ord.g]
  }
  
  cl.d=matrix(, M,G)
  if(D==1){
    cl.d=a.true
  }else{
    for (g in 1:G){
      cl.d[,g]=apply(a.true[,,g], 2, which.max)
    }
  }
  
  C <- array(0, c(D, D, N))
  mu.u <- matrix(0, D, N)
  u <- matrix(NA, D, N)
  for (i in 1:N) {
    C[,,i] <- diag(D)
    u[,i] <- mvrnorm(1, mu.u[,i], C[,,i])
  }
  
  p <- matrix(NA, nrow=N, ncol=M)
  for (i in 1:N){
    for (k in 1:M){
      lpred=b0[which.max(z[i,])] + t(a[, k, which.max(z[i,])]) %*% (mu0 + u[ ,i])
      p[i,k] <- 1/(1+exp(-lpred))
    }
  }
  
  y <- matrix(NA, nrow=N, ncol=M)
  for (i in 1:N){
    for (k in 1:M){
      y[i,k] <- rbinom(1, size = 1, prob = p[i,k])
    }
  }
  
  mod.sim <- mlta(X=y, DM=DM, G = G, D = D, beta0=NULL, nstarts=nstarts)
  
  list(mod.sim,cl.g,cl.d)
}


for(s in 1:S){
  
  a.sim[[s]] <- res[[s]][[1]]$a
  b.sim[[s]] <- res[[s]][[1]]$b
  mu.sim[[s]] <- res[[s]][[1]]$mu.d
  beta.sim[[s]] <- res[[s]][[1]]$beta
  se.sim[[s]]=res[[s]][[1]]$se
  
  # bic.sim[,,s] <- t(res[[s]][[1]]$BIC$`Table of BIC Results`) # for BIC simulation
  
  ord.g.sim = order(b.sim[[s]])
  cl.g.sim <- apply(res[[s]][[1]]$z[,ord.g.sim],1,which.max)
  rand.g[s] <- adjustedRandIndex(res[[s]][[2]],cl.g.sim)
  
  ord.d.sim = order(mu.sim[[s]])
  a.est=res[[s]][[1]]$a 
  a.est=a.est[ord.d.sim, , ]
  a.est=a.est[ , ,ord.g.sim]
  cl.d.sim=matrix(, M, G)
  for (g in 1:G){
    cl.d.sim[,g]=apply(a.est[,,g], 2, which.max)
  }
  
  rand.d[s]=adjustedRandIndex(res[[s]][[3]],cl.d.sim)
  print(rand.d)
}


# ARI
summary(rand.g)
summary(rand.d) 


# Beta

b.sim <- b.sim[!sapply(b.sim, function(x) all(is.na(x)))]
beta.sim <- beta.sim[!sapply(beta.sim, function(x) all(is.na(x)))]
mu.sim <- mu.sim[!sapply(mu.sim, function(x) all(is.na(x)))]

S=length(b.sim)

beta.sim2 = array(, c(2,G-1,S))
for (s in 1:S){
  print(s)
  if(!any(is.na(b.sim[[s]]))){
    ord = order(b.sim[[s]])
    be = cbind(0, beta.sim[[s]])
    be = be[,ord]
    beta.sim2[,,s] = be[,2:G]-be[,1]
  }
}

betao = cbind(0, beta)
ord = order(b0)
betao = betao[,ord]
print(betao[,2:G]-betao[,1])
print(apply(beta.sim2, c(1,2), mean))
print(apply(beta.sim2, c(1,2), median))

betao=betao[,2:G]-betao[,1]
par(mar=c(5,3.3,2,2))
par(mfrow=c(2,G-1))
par(cex.lab=2, cex.axis=1.5)
for (h in 1:2) {
  for (k in 1:(G-1)) {
    boxplot(beta.sim2[h,k,], #xlab=bquote(bold(beta)[paste(.(h-1),.(k+1))]),
            col="lightcyan3", border=NA, frame=F, 
            las=2, ylim=c(betao[h,k]-10,betao[h,k]+10))
    grid(lty="solid")
    boxplot(beta.sim2[h,k,], #xlab=bquote(bold(beta)[paste(.(h),.(k+1))]), 
            col="lightcyan3", border="black", frame=F,
            add=TRUE, las=2, ylim=c(betao[h,k]-1.5,betao[h,k]+1.5))
    abline(h=betao[h,k], col="red3",lwd=2, lty="longdash")
    mtext(expression(paste(italic("N="),2000," ")),side=3,outer=TRUE,padj=1.3, cex=1.3)
  }
}

for(h in 1:2){
  for (k in 1:(G-1)) {
    print(mean((beta.sim2[h,k,]-betao[h,k])^2))
  }
}


for (h in 1:2) {
  for (k in 1:(G-1)) {
    print( (mean(beta.sim2[h,k,])-betao[h,k])/betao[h,k])
  }
}


# b_g

b.sim2 = matrix(, G,S)
for (s in 1:S){
  ord = order(b.sim[[s]])
  b.sim2[,s] = b.sim[[s]][ord]
  bb= b.sim2[,s]-b.sim2[2,s]
  b.sim2[,s] = bb
}
ord = order(b0)
b2 = b0[ord]
b2 = b2-b2[2]
print(b2)
print(apply(b.sim2, 1, mean))
print(apply(b.sim2, 1, sd))
print(b2-apply(b.sim2, 1, mean))

par(mfrow=c(1,G))
for (k in 1:(G)) {
  boxplot(b.sim2[k,], main=paste("b",k), ylim=c(b2[k]-1.5,b2[k]+1.5))
  abline(h=b2[k], col=2)
}

for (k in 1:G) {
  print(mean((b.sim2[k,]-b2[k])^2))
}


for (k in 1:(G)) {
  print( (mean(b.sim2[k,])-b2[k])/b2[k])
}


# mu_d

mu.sim2 = matrix(, D,S)
for (s in 1:S){
  ord = order(mu.sim[[s]])
  mu.sim2[,s] = mu.sim[[s]][ord]
  mu.sim2[,s]=mu.sim2[,s]-mu.sim2[2,s]
}
ord = order(mu0)
mu2 = mu0[ord]
mu2 = mu2-mu2[2]
print(mu2)
print(apply(mu.sim2, 1, mean))
print(mu2-apply(mu.sim2, 1, mean))

par(mfrow=c(1,D))
for (k in 1:(D)) {
  boxplot(mu.sim2[k,], main=paste("mu",k), ylim=c(mu2[k]-1.5,mu2[k]+1.5))
  abline(h=mu2[k], col=2)
}

for (k in 1:D) {
  print(mean((mu.sim2[k,]-mu2[k])^2))
}


for (d in 1:D) {
  print( (mean(mu.sim2[d,])-mu2[d])/mu2[d])
}


# BIC
for(dg in 1:9){
  print(sum(apply(bic.sim[,,1:S],c(3),which.min)==dg)/S) 
}


# COVERAGE 

true_param=c(b,mu.d,beta)

coverage=list()
for(s in 1:S){
  confint=cbind( c((b.sim[[s]]-qnorm(0.975)*(se.sim[[s]][1:G])), 
                   (mu.sim[[s]]-qnorm(0.975)*se.sim[[s]][(G+1):(G+D)]),
                   (beta.sim[[s]]-qnorm(0.975)*se.sim[[s]][(G+D+1):(G+D+J*(G-1))])),
                 c((b.sim[[s]]+qnorm(0.975)*se.sim[[s]][1:G]), 
                   (mu.sim[[s]]+qnorm(0.975)*se.sim[[s]][(G+1):(G+D)]),
                   (beta.sim[[s]]+qnorm(0.975)*se.sim[[s]][(G+D+1):(G+D+J*(G-1))])) )
  coverage[[s]] <- (true_param >= confint[,1]) & (true_param <= confint[,2])
}


b.cov.sim=matrix(NA,G,S)
mu.cov.sim=matrix(NA,D,S)
beta.cov.sim=matrix(NA,J*(G-1),S)
for(s in 1:S){
  b.cov.sim[,s]=coverage[[s]][1:G]
  mu.cov.sim[,s]=coverage[[s]][(G+1):(G+D)]
  beta.cov.sim[,s]=coverage[[s]][(G+D+1):(G+D+J*(G-1))]
  
}

coverage_rate_b <- rowMeans(b.cov.sim) 
coverage_rate_mu <- rowMeans(mu.cov.sim) 
coverage_rate_beta <- rowMeans(beta.cov.sim) 

c(coverage_rate_b,coverage_rate_mu,coverage_rate_beta)

mean(coverage_rate_b)
mean(coverage_rate_mu)
mean(coverage_rate_beta)


# Standard errors

result <- optim(c(b,mu.d,c(beta)), log_likelihood, control = list(fnscale = -1, maxit=100), hessian = T)
library(numDeriv)
library(MASS)
score = grad(log_likelihood, x = result$par)
hessiana = result$hessian
se.true=sqrt(diag(ginv(hessiana) %*% score %*% t(score) %*% ginv(hessiana)))

se.true[(G+D+1):(G+D+(G-1)*J)]

b.se.sim=matrix(NA,G,S)
mu.se.sim=matrix(NA,D,S)
beta.se.sim=matrix(NA,J*(G-1),S)

for(s in 1:S){
  b.se.sim[,s]=se.sim[[s]][1:G]
  mu.se.sim[,s]=se.sim[[s]][(G+1):(G+D)]
  beta.se.sim[,s]=se.sim[[s]][(G+D+1):(G+D+J*(G-1))]
  
}

for(g in 1:G){
  print(mean(b.se.sim[g,])-se.true[1:G][g])
}

for(d in 1:D){
  print(mean(mu.se.sim[d,])-se.true[(G+1):(G+D)][d])
}

for(gj in 1:((G-1)*J)){
  print((mean(beta.se.sim[gj,])-se.true[(G+D+1):(G+D+(G-1)*J)][gj]))
}

## Application ####

### Data ####

load("appl_pediatric2.RData")

library(dplyr)

# summary(app_data$Alvarado_Score)
# app_data <- mutate(app_data, Alvarado_Score = ifelse(Alvarado_Score > 3, 1, 0))
# 
# summary(app_data$Paedriatic_Appendicitis_Score)
# app_data <- mutate(app_data, Paedriatic_Appendicitis_Score = ifelse(Paedriatic_Appendicitis_Score >= 6, 1, 0))
# 
# table(app_data$Appendix_on_US)
# app_data <- mutate(app_data, Appendix_on_US = ifelse(Appendix_on_US == "yes", 1, 0))
# 
# table(app_data$Migratory_Pain)
# app_data <- mutate(app_data, Migratory_Pain = ifelse(Migratory_Pain == "yes", 1, 0))
# 
# table(app_data$Lower_Right_Abd_Pain)
# app_data <- mutate(app_data, Lower_Right_Abd_Pain = ifelse(Lower_Right_Abd_Pain == "yes", 1, 0))
# 
# table(app_data$Contralateral_Rebound_Tenderness)
# app_data <- mutate(app_data, Contralateral_Rebound_Tenderness = ifelse(Contralateral_Rebound_Tenderness == "yes", 1, 0))
# 
# table(app_data$Coughing_Pain)
# app_data <- mutate(app_data, Coughing_Pain = ifelse(Coughing_Pain == "yes", 1, 0))
# 
# table(app_data$Nausea)
# app_data <- mutate(app_data, Nausea = ifelse(Nausea == "yes", 1, 0))
# 
# table(app_data$Loss_of_Appetite)
# app_data <- mutate(app_data, Loss_of_Appetite = ifelse(Loss_of_Appetite == "yes", 1, 0))
# 
# summary(app_data$Body_Temperature)
# app_data <- mutate(app_data, Body_Temperature = ifelse(Body_Temperature > 37.36, 1, 0))
# 
# summary(app_data$WBC_Count)
# app_data <- mutate(app_data, WBC_Count = ifelse(WBC_Count > 11 | WBC_Count < 4.5, 1, 0))
# 
# summary(app_data$Neutrophil_Percentage)
# app_data <- mutate(app_data, Neutrophil_Percentage = ifelse(Neutrophil_Percentage > 60 | Neutrophil_Percentage < 40, 1, 0))
# 
# table(app_data$Neutrophilia)
# app_data <- mutate(app_data, Neutrophilia = ifelse(Neutrophilia == "yes", 1, 0))
# 
# summary(app_data$RBC_Count)
# app_data <- mutate(app_data, RBC_Count = ifelse(RBC_Count > 5.5 | RBC_Count < 4, 1, 0))
# 
# summary(app_data$Hemoglobin)
# app_data <- mutate(app_data, Hemoglobin = ifelse(Hemoglobin < 10 | Hemoglobin > 15, 1, 0))
# 
# summary(app_data$RDW)
# app_data <- mutate(app_data, RDW = ifelse(RDW < 12.3 | RDW > 14.1, 1, 0))
# 
# summary(app_data$Thrombocyte_Count)
# app_data <- mutate(app_data, Thrombocyte_Count = ifelse(Thrombocyte_Count < 250 | 
#                                                           Thrombocyte_Count > 450, 1, 0))
# 
# summary(app_data$CRP)
# app_data <- mutate(app_data, CRP = ifelse(CRP > 10, 1, 0))
# 
# table(app_data$Dysuria)
# app_data <- mutate(app_data, Dysuria = ifelse(Dysuria == "yes", 1, 0))
# 
# table(app_data$Stool)
# app_data <- mutate(app_data, Stool = ifelse(Stool != "normal", 1, 0))
# 
# table(app_data$Peritonitis)
# app_data <- mutate(app_data, Peritonitis = ifelse(Peritonitis != "no", 1, 0))
# 
# table(app_data$Psoas_Sign)
# app_data <- mutate(app_data, Psoas_Sign = ifelse(Psoas_Sign == "yes", 1, 0))
# 
# table(app_data$US_Performed)
# app_data <- mutate(app_data, US_Performed = ifelse(US_Performed == "yes", 1, 0))

# table(app_data$Free_Fluids)
# app_data <- mutate(app_data, Free_Fluids = ifelse(Free_Fluids == "yes", 1, 0))

sum(app_data)/(nrow(app_data)*ncol(app_data)) # degree

### Descriptive plot ####

library(RColorBrewer)
par(mar=c(2,2,2,2))
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(800)
heatmap(m, Colv = NA, Rowv = NA, scale="none",col=colMain,margins=c(13,7),
        labCol=c("Alvarado Score","Paedriatic Appendicitis Score","Appendix on US",                  
                 "Migratory Pain","Lower Right Abd Pain","Contralateral Rebound Tenderness",
                 "Coughing Pain","Nausea","Loss of Appetite",                
                 "Body Temperature","WBC Count","Neutrophil Percentage",           
                 "Neutrophilia","RBC Count","Hemoglobin",                      
                 "RDW","Thrombocyte Count","CRP",                             
                 "Dysuria","Stool","Peritonitis",                     
                 "Psoas Sign","US Performed","Free Fluids"))
legend(x="topright", legend=c("0", "1"),fill=c("#F7FBFF","#08156B"),box.col="white",
       title="Values")

# eps 1000 1000
library(ggplot2)
library(reshape2)
library(RColorBrewer)
m_long <- melt(m)
m_long$value=as.factor(m_long$value)
fill_colors <- c("0" = "#F7FBFF", "1" = "#08156B")
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(800)
ggplot(data = m_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = fill_colors, labels = c("0", "1"), drop = FALSE) +
  labs(x = NULL, y = NULL, fill = "Values") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),  
        legend.title = element_text(size = 20)) +
  scale_x_discrete(labels = c("Alvarado Score","Paedriatic Appendicitis Score","Appendix on US",                  
                              "Migratory Pain","Lower Right Abd Pain","Contralateral Rebound Tenderness",
                              "Coughing Pain","Nausea","Loss of Appetite",                
                              "Body Temperature","WBC Count","Neutrophil Percentage",           
                              "Neutrophilia","RBC Count","Hemoglobin",                      
                              "RDW","Thrombocyte Count","CRP",                             
                              "Dysuria","Stool","Peritonitis",                     
                              "Psoas Sign","US Performed","Free Fluids")) +
  guides(fill = guide_legend(title = "Values"))


### Variables plots ####

library(ggplot2)
library(gridExtra)
library(patchwork)

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Alvarado_Score)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot1=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "Proportion") +
  ggtitle("Alvarado Score") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Paedriatic_Appendicitis_Score)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot2=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Paediatr Append Score") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Appendix_on_US)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot3=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Appendix on US") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Migratory_Pain)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot4=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Migratory Pain") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Lower_Right_Abd_Pain)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot5=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("LowerRight Abd Pain") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Contralateral_Rebound_Tenderness)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot6=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Contr Rebound \n Tenderness") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Coughing_Pain)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot7=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "Proportion") +
  ggtitle("Coughing Pain") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Nausea)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot8=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Nausea") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Loss_of_Appetite)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot9=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Loss of Appetite") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Body_Temperature)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot10=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Body Temperature") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$WBC_Count)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot11=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("WBC Count") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Neutrophil_Percentage)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot12=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Neutrophil %") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Neutrophilia)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot13=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "Proportion") +
  ggtitle("Neutrophilia") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$RBC_Count)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot14=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("RBC Count") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Hemoglobin)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot15=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Hemoglobin") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$RDW)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot16=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("RDW") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Thrombocyte_Count)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot17=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Thrombocyte Count") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$CRP)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot18=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("CRP") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Dysuria)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot19=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "Proportion") +
  ggtitle("Dysuria") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Stool)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot20=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Stool") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Peritonitis)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot21=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Peritonitis") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Psoas_Sign)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot22=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Psoas Sign") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(cbind(0,table(app_data$US_Performed))/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot23=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("US Performed") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Index = c("0", "1"),
  Value = as.vector(table(app_data$Free_Fluids)/nrow(app_data))
)
myorder <- c("0", "1")
data$Index <- factor(data$Index, levels = myorder)
plot24=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "", y = "") +
  ggtitle("Free Fluids") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

plot = plot1+plot2+plot3+plot4+plot5+plot6+plot7+plot8+
  plot9+plot10+plot11+plot12+plot13+plot14+plot15+
  plot16+plot17+plot18+plot19+plot20+plot21+
  plot22+plot24 + plot_layout(ncol = 6)
print(plot) # eps: 1200, 900


### Covariates plots ####
library(ggplot2)
library(gridExtra)

summary(DM$Age)
plotc1=ggplot(DM, aes(x = Age, y = ..count../sum(..count..))) +
  geom_histogram(binwidth = 0.5, fill = "gray70", color = "gray60") +
  labs(title = "Age") +
  labs(x = "", y = "Proportion") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5))


summary(DM$BMI)
plotc2=ggplot(DM, aes(x = BMI, y = ..count../sum(..count..))) +
  geom_histogram(binwidth = 0.5, fill = "gray70", color = "gray60") +
  labs(title = "BMI") +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5))


data <- data.frame(
  Index = c("female", "male"),
  Value = as.vector(table(DM$Sex)/nrow(DM))
)
myorder <- c("female", "male")
data$Index <- factor(data$Index, levels = myorder)
plotc3=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70", color = "gray60") +
  labs(x = "", y = "Proportion") +
  ggtitle("Gender") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


summary(DM$Height)
plotc4=ggplot(DM, aes(x = Height, y = ..count../sum(..count..))) +
  geom_histogram(binwidth = 3, fill = "gray70", color = "gray60") +
  labs(title = "Height") +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5))


summary(DM$Weight)
plotc5=ggplot(DM, aes(x = Weight, y = ..count../sum(..count..))) +
  geom_histogram(binwidth = 3, fill = "gray70", color = "gray60") +
  labs(title = "Weight") +
  labs(x = "", y = "Proportion") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5))


summary(DM$Length_of_Stay)
plotc6=ggplot(DM, aes(x = Length_of_Stay, y = ..count../sum(..count..))) +
  geom_histogram(binwidth = 1, fill = "gray70", color = "gray60") +
  labs(title = "Length of Stay") +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5))


data <- data.frame(
  Index = c("conservative","primary surgical","secondary surgical"),
  Value = as.vector(table(DM$Management)/nrow(DM))
)
myorder <- c("conservative","primary surgical","secondary surgical")
data$Index <- factor(data$Index, levels = myorder)
plotc7=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70", color = "gray60") +
  labs(x = "", y = "Proportion") +
  ggtitle("Management") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Index = c("complicated","uncomplicated"),
  Value = as.vector(table(DM$Severity)/nrow(DM))
)
myorder <- c("uncomplicated","complicated")
data$Index <- factor(data$Index, levels = myorder)
plotc8=ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70", color = "gray60") +
  labs(x = "", y = "") +
  ggtitle("Severity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

plot = plotc1+plotc2+plotc3+plotc4+plotc5+plotc6+plotc7+plotc8 + plot_layout(ncol = 2)
print(plot) # eps: 1200, 900

### Model ####

library(MASS)

# head(DM)
# str(DM)
# DM$BMI=as.numeric(DM$BMI)
# DM$Sex=as.factor(DM$Sex)
# DM$Management=as.factor((DM$Management))
# DM$Severity=as.factor(DM$Severity)
# str(DM)
# DesMat <- model.matrix(~ Age + BMI + Sex + Height + 
#                          Weight + Length_of_Stay + Management + 
#                          Severity, data=DM)
# dim(DesMat)

m=as.matrix(app_data)
modg1d1=lta(X=m, D=1, nstarts=10)
modg1d2=lta(X=m, D=2, nstarts=10)
modg1d3=lta(X=m, D=3, nstarts=10)
modg1d4=lta(X=m, D=4, nstarts=10)

modg2d1=mlta(X=m, DM=DesMat, G = 2, D = 1, beta0=NULL, nstarts=10)
modg2d2=mlta(X=m, DM=DesMat, G = 2, D = 2, beta0=NULL, nstarts=10)
modg2d3=mlta(X=m, DM=DesMat, G = 2, D = 3, beta0=NULL, nstarts=10)
modg2d4=mlta(X=m, DM=DesMat, G = 2, D = 4, beta0=NULL, nstarts=10)

modg3d1=mlta(X=m, DM=DesMat, G = 3, D = 1, beta0=NULL, nstarts=10)
modg3d2=mlta(X=m, DM=DesMat, G = 3, D = 2, beta0=NULL, nstarts=10)
modg3d3=mlta(X=m, DM=DesMat, G = 3, D = 3, beta0=NULL, nstarts=10)
modg3d4=mlta(X=m, DM=DesMat, G = 3, D = 4, beta0=NULL, nstarts=10)

modg4d1=mlta(X=m, DM=DesMat, G = 4, D = 1, beta0=NULL, nstarts=10)
modg4d2=mlta(X=m, DM=DesMat, G = 4, D = 2, beta0=NULL, nstarts=10)
modg4d3=mlta(X=m, DM=DesMat, G = 4, D = 3, beta0=NULL, nstarts=10)
modg4d4=mlta(X=m, DM=DesMat, G = 4, D = 4, beta0=NULL, nstarts=10)

mod=modg2d3
G=2
D=3
a = mod$a
npoints <- 7           
Q <- npoints ^ D       
GaHer <- glmmML::ghq(npoints, FALSE)                     
ugh <- as.matrix(expand.grid(rep(list(GaHer$zeros), D))) 
ugh.star <- sqrt(2)*ugh                                  
p.gh = as.matrix(expand.grid(rep(list(GaHer$weights), times = D)))	
W =  (2)^(D/2) * exp(apply(ugh, 1, crossprod)) * apply(p.gh, 1, prod)
W.prod =  apply(p.gh, 1, prod)
Phi <- apply(ugh.star, 1, mvtnorm::dmvnorm)                
X = m
DM = DesMat

log_likelihood <- function(bb) {
  maxb=G
  maxmu=maxb+D
  maxbeta=maxmu+(J*(G-1))
  fy.uz=array(0,c(N,Q,G))
  if(D==1){
    for (g in 1:G)
    {
      Agh <- t(tcrossprod(matrix(a[, , g],M,D), t(t(ugh.star)+bb[(maxb+1):(maxmu)])) + bb[1:maxb][g])
      pgh <- 1 / (1 + exp(-Agh))
      fy.uz[, , g] <-
        exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh)))
      # fy.uz[is.nan(fy.uz[, , g])] <- 0
    }
  }else{
    for (g in 1:G)
    {
      Agh <- t(tcrossprod(t(a[, , g]), t(t(ugh.star)+bb[(maxb+1):(maxmu)])) + bb[1:maxb][g])
      pgh <- 1 / (1 + exp(-Agh))
      fy.uz[, , g] <-
        exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh)))
      # fy.uz[is.nan(fy.uz[, , g])] <- 0
    }
  }
  fy.z <- apply(aperm(apply(fy.uz,c(1,3),function(x) x*W*Phi),c(2,1,3)),c(1,3),sum)
  
  betaaa=matrix(bb[(maxmu+1):maxbeta],J,G-1)
  exb <- exp(DM %*% betaaa)
  eta <- cbind(1,exb)/(rowSums(exb)+1)
  
  fyz <- eta * (fy.z)
  fy <- apply(fyz, 1, sum)
  loglik <- sum(log(fy))
  return(loglik)
}

b = mod$b
mu.d = mod$mu.d
beta = mod$beta

J=ncol(DesMat)
N=nrow(m)
npoints <- 7           
Q <- npoints ^ D
result <- optim(c(b,mu.d,c(beta)), log_likelihood, control = list(fnscale = -1, maxit=100), hessian = T)
library(numDeriv)
library(MASS)
score = grad(log_likelihood, x = result$par)
hessiana = result$hessian
se=sqrt(diag(ginv(hessiana) %*% score %*% t(score) %*% ginv(hessiana)))
se


### Results ####

mod=modg2d3

mod$b
mod$mu.d

# z
library(ggplot2)
table(apply(mod$z,1,which.max))
ord.g=order(mod$b)
data <- data.frame(
  Index = c("Component 1", "Component 2"),
  Value = as.vector(table(apply(mod$z,1,which.max))[ord.g]/nrow(m))
)
myorder <- c("Component 1", "Component 2")
data$Index <- factor(data$Index, levels = myorder)
ggplot(data, aes(x = factor(Index), y = Value)) +
  geom_bar(stat = "identity", fill = "gray70", color = "gray60") +
  labs(x = "", y = "Proportion") +
  ggtitle("") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


### Ordered data matrix ####
cl.g=apply(mod$z,1,which.max)
tabu = c(0, cumsum(table(cl.g)))
clU = cl.g
ordu = order(clU)
D=3
G=2
N=nrow(m)
C <- array(0, c(D, D, N))
mu.u <- matrix(0, D, N)
u <- matrix(NA, D, N)
for (i in 1:N) {
  C[,,i] <- diag(D)
  u[,i] <- MASS::mvrnorm(1, mu.u[,i], C[,,i])
}
M=ncol(m)
p <- matrix(NA, nrow=N, ncol=M)
for (i in 1:N){
  for (k in 1:M){
    lpred=mod$b[which.max(mod$z[i,])] + t(mod$a[, k, which.max(mod$z[i,])]) %*% (mod$mu.d + u[ ,i])
    p[i,k] <- 1/(1+exp(-lpred))
  }
}
rownames(p)=rownames(m)
colnames(p)=colnames(m)
y=p # y=m 
yy = y[ordu,]
ord.m=order(mod$mu.d)
M=ncol(m)
cl.d=matrix(, M,G)
a=mod$a[ord.m,,]
for (gg in 1:G){
  cl.d[,gg]=apply(a[,,gg], 2, which.max)
}
clV = cl.d
idxU = rep(1:G, each=M)
ordV = c()
yy2 = c()
clVk = c()
aaa = c()
for(gg in 1:G){
  whV = clV[idxU==gg]
  tabv = c(0, cumsum(table(whV)))
  #if(t==1) print(tabv)
  ordV = order(whV)
  clVk = c(clVk, whV)
  #aaa = c(aaa, paste(k, whV[ordV], sep = "_"))
  yy2 = rbind(yy2, yy[(tabu[gg]+1):tabu[gg+1],ordV])
}

library(RColorBrewer)
par(mar=c(2,2,2,2))
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(800)
heatmap(yy2, Colv = NA, Rowv = NA, scale="none",col=colMain,margins=c(11,11))
colnames(yy2) # per riordinare le etichette del plot
heatmap(yy2, Colv = NA, Rowv = NA, scale="none",col=colMain,margins=c(13,7),
        labCol=c("Migratory Pain","Contralateral Rebound Tenderness","Coughing Pain",
                 "RBC Count","Hemoglobin","RDW","Thrombocyte Count","Dysuria",
                 "Stool","Peritonitis","Psoas Sign","Paedriatic Appendicitis Score",   
                 "Nausea","Loss of Appetite","Body Temperature","WBC Count","Neutrophilia",
                 "CRP","Free Fluids","Alvarado Score","Appendix on US","Lower Right Abd Pain",
                 "Neutrophil Percentage","US Performed"))

# eps 1000 1000
library(ggplot2)
library(reshape2)
library(RColorBrewer)
yy2_long <- melt(yy2)
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(800)
col_labels <- c("Migratory Pain","Contralateral Rebound Tenderness","Coughing Pain",
                "RBC Count","Hemoglobin","RDW","Thrombocyte Count","Dysuria",
                "Stool","Peritonitis","Psoas Sign","Paedriatic Appendicitis Score",   
                "Nausea","Loss of Appetite","Body Temperature","WBC Count","Neutrophilia",
                "CRP","Free Fluids","Alvarado Score","Appendix on US","Lower Right Abd Pain",
                "Neutrophil Percentage","US Performed")
ggplot(data = yy2_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = colMain) +
  labs(x = NULL, y = NULL, fill = "Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = col_labels) +
  scale_y_discrete(labels = rownames(yy2))

# NEW
#519 519 size

cl.g=apply(mod$z,1,which.max)
tabu = c(0, cumsum(table(cl.g)))
clU = cl.g
ordu = order(clU)
D=3
G=2
N=nrow(m)
M=ncol(m)
p <- matrix(NA, nrow=N, ncol=M)
rownames(p)=rownames(m)
colnames(p)=colnames(m)
y=p # y=m 
yy = y[ordu,]
ord.m=order(mod$mu.d)
M=ncol(m)
cl.d=matrix(, M,G)
a=mod$a[ord.m,,]
for (gg in 1:G){
  cl.d[,gg]=apply(a[,,gg], 2, which.max)
}
clV = cl.d
idxU = rep(1:G, each=M)
ordV = c()
yy2 = c()
clVk = c()
aaa = c()

# Order row clusters
for(gg in rev(1:G)){  
  whV = clV[idxU==gg]
  tabv = c(0, cumsum(table(whV)))
  ordV = order(whV)
  
  rows_in_group = (tabu[gg]+1):tabu[gg+1]
  yy2 = rbind(yy2, yy[rows_in_group, ordV])  
}

library(RColorBrewer)
par(mar=c(2,2,2,2))
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(800)
heatmap(yy2, Colv = NA, Rowv = NA, scale="none",col=colMain,margins=c(11,11))
colnames(yy2) 
heatmap(yy2, Colv = NA, Rowv = NA, scale="none",col=colMain,margins=c(13,7),
        labCol=c("Migratory Pain","Contralateral Rebound Tenderness","Coughing Pain",
                 "RBC Count","Hemoglobin","RDW","Thrombocyte Count","Dysuria",
                 "Stool","Peritonitis","Psoas Sign","Paedriatic Appendicitis Score",   
                 "Nausea","Loss of Appetite","Body Temperature","WBC Count","Neutrophilia",
                 "CRP","Free Fluids","Alvarado Score","Appendix on US","Lower Right Abd Pain",
                 "Neutrophil Percentage","US Performed"))


### Items classification ####
ord.g=order(mod$b)
ord.d=order(mod$mu.d)
a.new=mod$a[ord.d,,ord.g]
class=apply(a.new,c(2,3),which.max)
rownames(class)=colnames(m)
class
library(knitr)
kable(class, format = "latex", col.names = c("Component 1", "Component 2"))

### beta ####
beta2=mod$beta
se2=sqrt(se)
lb=beta2-1.96*se2[6:15]
ub=beta2+1.96*se2[6:15]

dat=cbind(beta2,lb,ub)
rownames(dat) <- sub("class", "", rownames(dat))
dat=cbind(c(rownames(dat)),dat)
colnames(dat)=c("category","estimate","lower","upper")
dat=as.data.frame(dat)
dat$estimate=as.numeric(dat$estimate)
dat$lower=as.numeric(dat$lower)
dat$upper=as.numeric(dat$upper)
dat$category=as.character(dat$category)

library(knitr)
dat$formatted_interval <- sprintf("(%.2f, %.2f)", dat$lower, dat$upper)
dat <- dat[, c("estimate", "formatted_interval")]
kable(dat, format = "latex", col.names = c("Estimates", "Confidence interval"))

library(plotrix)
beta_2 <- mod$beta
u.beta2 <- ub
l.beta2 <- lb
y <- 1:length(mod$beta)
par(mfrow=c(1,1))
par(mar=c(4.3,14,2.1,0), xpd = T)
plotCI(y=y, x=beta_2, ui=u.beta2, li=l.beta2, xlab = "", ylab = "", pch = 19, col="gray40", main = "Component 2 vs 1",
       cex.lab=1, cex.axis=1, cex.main=1.5, cex.sub=1, cex=1, xlim=c(-2,2), lwd=2, yaxt="n", err="x", frame.plot=F)
par(mar=c(4.3,17,2.1,0), xpd = F)
grid(lty="solid")
plotCI(y=y, x=beta_2, ui=u.beta2, li=l.beta2, xlab = "", ylab = "", pch = 19, col="gray40", main = "Component 2 vs 1",
       cex.lab=1, cex.axis=1, cex.main=1.5, cex.sub=1, cex=1, xlim=c(-2,2), lwd=2, add=TRUE, err="x")
axis(2, at=1:25, labels=c(rownames(dat)), las=1, cex.axis=1)
abline(v=0, lty=2)
mtext(expression(beta[2]),side=1,las=1,line=3.5, cex=1.5)


### Predicted probabilities ####
PI=array(0,c(G,D,N))
for (gg in 1:G) {
  for(d in 1:D) {
    class=apply(mod$z,1,which.max)
    id=which(class==gg)
    id.k=which(mod$a[d,,gg]==1)
    PI[gg,d,]=1/(1+exp(-(mod$b[gg]+mod$a[d,id.k,gg][1]*(mod$mu.d[d]+u[d,]))))
  }
}
order(mod$b)
order(mod$mu.d)

PI11=PI[1,3,]
PI12=PI[1,1,]
PI13=PI[1,2,]
PI21=PI[2,3,]
PI22=PI[2,1,]
PI23=PI[2,2,]

library(ggplot2)

data <- data.frame(
  Value = c(PI11, PI12, PI13, PI21, PI22, PI23),
  Group = rep(c("g=1 d=1", "g=1 d=2", "g=1 d=3", "g=2 d=1", "g=2 d=2", "g=2 d=3"), each = length(PI11))
)
means <- aggregate(data$Value, list(Group = data$Group), mean)
data <- merge(data, means, by = "Group")
colnames(data)[ncol(data)] <- "Mean"
ggplot(data, aes(x = Value, fill = Group)) +
  geom_histogram(binwidth = 0.05, color = "gray40") +
  geom_vline(data = data, aes(xintercept = Mean), color = "firebrick", linetype = "dashed", size = 1) +
  facet_wrap(~Group, nrow = 2, ncol = 3, scales = "free") +
  scale_fill_manual(values = rep("gray85", 6)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(strip.text = element_text(face = "italic", size = 15),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") + 
  coord_cartesian(xlim = c(0, 1))

ggplot(data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(color = "gray40", width = 0.5) +
  geom_point(stat = "summary", fun = mean, color = "firebrick", size = 3, shape = 18) +
  facet_wrap(~Group, nrow = 2, ncol = 3, scales = "free") +
  scale_fill_manual(values = rep("gray85", 6)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(face = "italic", size = 16, color = "gray40"),  
        axis.title = element_blank(),  
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1))


## Original MLTA ####

# # Remove the biclustering functions:
# rm(f_lca); rm(f_lca_nstarts); rm(f_lta); rm(f_lta_nstarts); rm(f_mlta); rm(f_mlta_methods);
# rm(f_mlta_nstarts); rm(lca); rm(lta); rm(mlta); rm(print.mlta); rm(print.mmlta); rm(ResTable);
# rm(tableBIC)

library(lvm4net)
G <- 1:4
D <- 1:4

mod.mlta <- mlta(m, G = G, D = D, wfix = FALSE, nstarts=10) 
mod.mlta.wfix <- mlta(m, G = G, D = 1:4, wfix = TRUE, nstarts=10) 

which.min(mod.mlta$BIC$`Table of BIC Results`)
mod.mlta$BIC
which.min(mod.mlta.wfix$BIC$`Table of BIC Results`)
mod.mlta.wfix$BIC # G=2, D=1

MLTA <- mod.mlta.wfix[[5]] 

### Predicted probabilities ####

D=1
G=2
N=nrow(m)
C <- array(0, c(D, D, N))
mu.u <- matrix(0, D, N)
u <- matrix(NA, D, N)
for (i in 1:N) {
  C[,,i] <- diag(D)
  u[,i] <- MASS::mvrnorm(1, mu.u[,i], C[,,i])
}
M=ncol(m)
PI=array(0,c(G,M,N))
for (gg in 1:G) {
  for( k in 1:M){
    class=apply(MLTA$z,1,which.max)
    id=which(class==gg)
    #id.k=which(mod$a[d,,gg]==1)
    PI[gg,k,]=1/(1+exp(-(MLTA$b[gg,k]+MLTA$w[1,k]*(u[1,])))) 
  }
}
order(apply(MLTA$b,1, which.max))
library(tidyverse)
df_long <- PI %>%
  as.data.frame.table() %>%
  rename(Group = Var1, Column = Var2, Obs = Var3, Value = Freq)
df_long <- df_long %>%
  mutate(Group = factor(Group, labels = paste("Component", 1:G)),
         Column = factor(Column, labels = paste("Column", 1:M)))
ggplot(df_long, aes(x = Column, y = Value, fill=Group)) +
  geom_boxplot(color = "gray40", width = 0.5) +
  geom_point(stat = "summary", fun = mean, color = "firebrick", size = 3, shape = 18) +
  facet_grid(Group ~ ., labeller = labeller(Group = c("Component 1" = "Component 1", 
                                                      "Component 2" = "Component 2"))) +
  scale_fill_manual(values = rep("gray85", 46)) +
  scale_x_discrete(labels = c("Column 1" = "Alvarado Score","Column 2" = "Pedriatic Appendicitis Score",
                              "Column 3" = "Appendix on US","Column 4" = "Migratory Pain",
                              "Column 5" = "Lower Right Abd Pain","Column 6" = "Contralateral Rebound Tenderness",
                              "Column 7" = "Coughing Pain","Column 8" = "Nausea",
                              "Column 9" = "Loss of Appetite","Column 10" = "Body Temperature",
                              "Column 11" = "WBC Count","Column 12" = "Neutrophil Percentage",
                              "Column 13" = "Neutrophilia","Column 14" = "RBC Count",
                              "Column 15" = "Hemoglobin","Column 16" = "RDW",
                              "Column 17" = "Thrombocyte Count","Column 18" = "CRP",
                              "Column 19" = "Dysuria","Column 20" = "Stool",
                              "Column 21" = "Peritonitis","Column 22" = "Psoas Sign",
                              "Column 23" = "Free Fluids")) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none", strip.text = element_text(size = 15)) + 
  coord_cartesian(ylim = c(0, 1))


### Covariates plots ####

library(ggplot2)
library(gridExtra)
DM$cluster=apply(MLTA$z,1,which.max)
DM$cluster=ifelse(DM$cluster==1,"Component 1",DM$cluster);DM$cluster=ifelse(DM$cluster==2,"Component 2",DM$cluster)

plotc1=ggplot(DM, aes(x = Age, y = ..count../sum(..count..), fill=cluster)) +
  geom_histogram(binwidth = 0.5, fill = "gray70", color = "gray60") +
  facet_wrap(~ cluster, scales = "fixed") +
  labs(title = "Age") +
  labs(x = "", y = "Proportion") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size = 15))


plotc2=ggplot(DM, aes(x = BMI,  y = after_stat(count)/sum(after_stat(count)), fill=cluster)) +
  geom_histogram(binwidth = 0.5, fill = "gray70", color = "gray60") +
  facet_wrap(~ cluster, scales = "fixed") +
  labs(title = "BMI") +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size = 15))


plotc3=ggplot(DM, aes(x = Sex, fill=cluster)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill = "gray70", color = "gray60") +
  facet_wrap(~ cluster, scales = "fixed") +
  labs(title = "Gender") +
  labs(x = "", y = "Proportion") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 1))

plotc4=ggplot(DM, aes(x = Height, y = ..count../sum(..count..), fill=cluster)) +
  geom_histogram(binwidth = 3, fill = "gray70", color = "gray60") +
  facet_wrap(~ cluster, scales = "fixed") +
  labs(title = "Height") +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size = 15))

plotc5=ggplot(DM, aes(x = Weight, y = ..count../sum(..count..), fill=cluster)) +
  geom_histogram(binwidth = 3, fill = "gray70", color = "gray60") +
  facet_wrap(~ cluster, scales = "fixed") +
  labs(title = "Weight") +
  labs(x = "", y = "Proportion") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size = 15))


plotc6=ggplot(DM, aes(x = Length_of_Stay, y = ..count../sum(..count..), fill=cluster)) +
  geom_histogram(binwidth = 3, fill = "gray70", color = "gray60") +
  facet_wrap(~ cluster, scales = "fixed") +
  labs(title = "Length of Stay") +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size = 15))


myorder <- c("conservative","primary surgical","secondary surgical")
DM$Management <- factor(DM$Management, levels = myorder)
plotc7=ggplot(DM, aes(x = Management, fill=cluster)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill = "gray70", color = "gray60") +
  facet_wrap(~ cluster, scales = "fixed") +
  labs(title = "Management") +
  labs(x = "", y = "Proportion") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15, angle=45, hjust=1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 1))


myorder <- c("uncomplicated","complicated")
DM$Severity <- factor(DM$Severity, levels = myorder)
plotc8=ggplot(DM, aes(x = Severity, fill=cluster)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), fill = "gray70", color = "gray60") +
  facet_wrap(~ cluster, scales = "fixed") +
  labs(title = "Severity") +
  labs(x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 15, angle=45, hjust=1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 1))

grid.arrange(plotc1, plotc2, plotc3, plotc4, plotc5, plotc6, plotc7, plotc8, ncol = 2)
# eps: 1690 1300


### Ordered data matrix ####

cl.g=apply(MLTA$z,1,which.max)
tabu = c(0, cumsum(table(cl.g)))
clU = cl.g
ordu = order(clU)
D=1
G=2
N=nrow(m)
M=ncol(m)
p <- matrix(NA, nrow=N, ncol=M)
for (i in 1:N){
  for (k in 1:M){
    lpred=MLTA$b[which.max(MLTA$z[i,]),k] + (MLTA$w[1, k] * (u[1,i]))
    p[i,k] <- 1/(1+exp(-lpred))
  }
}
rownames(p)=rownames(m)
colnames(p)=colnames(m)
y=p # y=m 
yy2 = y[ordu,]

library(RColorBrewer)
par(mar=c(2,2,2,2))
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(800)
heatmap(yy2, Colv = NA, Rowv = NA, scale="none",col=colMain,margins=c(11,11))
colnames(yy2) # per riordinare le etichette del plot
heatmap(yy2, Colv = NA, Rowv = NA, scale="none",col=colMain,margins=c(13,7),
        labCol=c("Migratory Pain","Contralateral Rebound Tenderness","Coughing Pain",
                 "RBC Count","Hemoglobin","RDW","Thrombocyte Count","Dysuria",
                 "Stool","Peritonitis","Psoas Sign","Appendicitis Score",   
                 "Nausea","Loss of Appetite","Body Temperature","WBC Count","Neutrophilia",
                 "CRP","Free Fluids","Alvarado Score","Appendix on US","Lower Right Abd Pain",
                 "Neutrophil Percentage","US Performed"),
        cexCol = 1.5, cexRow = 1.5)
segments(x0 = 0.1,
         y0 = 0.68,
         x1 = 0.82,
         y1 = 0.68,
         col = "red", lwd = 4)


# NEW

cl.g = apply(MLTA$z,1,which.max)
tabu = c(0, cumsum(table(cl.g)))
clU = cl.g
ordu = order(clU)
D = 1
G = 2
N = nrow(MLTA$z)
M = ncol(MLTA$w)
p <- matrix(NA, nrow = N, ncol = M)
for (i in 1:N) {
  for (k in 1:M) {
    lpred = MLTA$b[which.max(MLTA$z[i,]), k] + (MLTA$w[1, k] * (u[1,i]))
    p[i, k] <- 1 / (1 + exp(-lpred))
  }
}
rownames(p) = rownames(MLTA$z)
y = p # y = m
yy = y[ordu,]
yy2 = c()
for (gg in rev(1:G)) {  
  rows_in_group = (tabu[gg] + 1):tabu[gg + 1]  
  yy2 = rbind(yy2, yy[rows_in_group, ])  
}

library(RColorBrewer)
par(mar = c(2,2,2,2))
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(800)
heatmap(yy2, Colv = NA, Rowv = NA, scale = "none", col = colMain, margins = c(11,11))
heatmap(yy2, Colv = NA, Rowv = NA, scale = "none", col = colMain, margins = c(13,7),
        labCol = c("Migratory Pain","Contralateral Rebound Tenderness","Coughing Pain",
                   "RBC Count","Hemoglobin","RDW","Thrombocyte Count","Dysuria",
                   "Stool","Peritonitis","Psoas Sign","Appendicitis Score",   
                   "Nausea","Loss of Appetite","Body Temperature","WBC Count","Neutrophilia",
                   "CRP","Free Fluids","Alvarado Score","Appendix on US","Lower Right Abd Pain",
                   "Neutrophil Percentage","US Performed"),
        cexCol = 1.5, cexRow = 1.5)
