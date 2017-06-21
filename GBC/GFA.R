# This file contains two functions GFA_MCMC and GFA_EM
# written by Changgee Chang
# version: 20170606
#
#
# function GFA_EM: runs EM for GFA
#
# Parameters:
# X: p by n data matrix where p is the number of covariates and n is the number of samples
# type: p by 1 vector of data types (see below)
# E: e by 2 matrix for edges where e is the number of edges.
#    E must be sorted by the first column and then the second column.
#    Both (j,k) and (k,j) must be added and therefore e is twice the actual number of edges.
# L: number of factors
# v0,v1: variances of spike and slab for W
# lam: shrinkage parameter for W
# eta: smoothing parameter for W
# param: p by 1 vector of parameters nu, n, r, or N (see data types below)
# W.init: initial value of W
# GBC: if TRUE, run GBC instead of GFA
# u0,u1: if GBC, variances of spike and slab for Z
# zeta: shrinkage parameter for Z
# PX: If TRUE, use the PX-EM algorithm.
#
# Data types:
# 0: continuous Gaussian - nu is required as param
# 1: Binomial - n_j must be provided as param
# 2: Negative Binomial - r_j must be provided as param
# 3: Poisson - N must be provided as param

GFA_EM <- function(X,type,E,L,v0,v1,lam,eta,param,intercept=T,smoothing="MRF",W.init=NULL,GBC=F,u0=v0,u1=v1,zeta=lam,PX=F)
{
  p = nrow(X)
  n = ncol(X)
  
  # edge
  e = nrow(E)
  Eidx = rep(1,p+1)
  for ( j in 1:p )
  {
    Eidx[j+1] = Eidx[j]
    while ( Eidx[j+1] <= e )
    {
      if ( E[Eidx[j+1],1] == j )
        Eidx[j+1] = Eidx[j+1] + 1
      else
        break
    }
  }
  
  # normalization
  
  # initialization
  Y = matrix(0,p,n)
  psi = matrix(0,p,n)
  kappa = matrix(0,p,n)
  b = matrix(0,p,n)
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
    {
      Y[j,] = X[j,]
      psi[j,] = X[j,]
    }
    if ( type[j] == 1 )
    {
      pbar = (X[j,]+1)/(param[j]+2)
      Y[j,] = log(pbar/(1-pbar))
      kappa[j,] = X[j,]-param[j]/2
      b[j,] = param[j]
    }
    if ( type[j] == 2 )
    {
      pbar = param[j]/(param[j]+X[j,]+1)
      Y[j,] = log(pbar/(1-pbar))
      kappa[j,] = (X[j,]-param[j])/2
      b[j,] = param[j]+X[j,]
    }
    if ( type[j] == 3 )
    {
      Y[j,] = log(pmax(1,X[j,]))
      psi[j,] = log(param[j])
      kappa[j,] = X[j,]-param[j]/2
      b[j,] = param[j]
    }
  }
  
  if ( !intercept )
  {
    m = rep(0,p)
    for ( j in 1:p )
      if ( type[j] == 1 )
        m[j] = log(1/49)
  }
  else
  {
    m = apply(Y,1,median)
    for ( j in 1:p )
    {
      if ( type[j] == 1 )
      {
        if ( m[j] == 0 )
          pr = 1/50
        else if ( m[j] == param[j] )
          pr = 49/50
        else
          pr = m[j]/param[j]
        m[j] = log(pr/(1-pr))
      }
      else if ( type[j] == 3 )
        m[j] = log(max(1,m[j]))
    }
  }
  
  if ( !is.null(W.init) )
  {
    W = W.init
    WW = t(W)%*%W
    Z = chol2inv(chol(WW)) %*% t(W) %*% (Y-m)
  }
  else
  {
    init_svd = svd(Y-m,L,L)
    W = init_svd$u %*% diag(sqrt(init_svd$d[1:L]),L)    
    Z = diag(sqrt(init_svd$d[1:L]),L) %*% t(init_svd$v)
  }

  rho = matrix(0,p,n)
  alpha = matrix(0,p,L)
  thetaW = 1/(1+exp(-alpha))
  if ( GBC )
    thetaZ = matrix(0.5,L,n)
  
  
  iv0 = 1/v0
  iv1 = 1/v1
  div = iv1-iv0
  lv0 = log(v0)
  lv1 = log(v1)
  dlv = lv1-lv0

  iu0 = 1/u0
  iu1 = 1/u1
  diu = iu1-iu0
  lu0 = log(u0)
  lu1 = log(u1)
  dlu = lu1-lu0
  
  iter = 0
  while (T)
  {
    ptheta = thetaW
    iter = iter + 1
    
    # E-step for rho
    mu = m + W%*%Z
    for ( j in 1:p )
    {
      if ( type[j] == 0 )
      {
        a_rho = (param[j]+n)/2
        b_rho = param[j]/2 + sum((mu[j,]-psi[j,])^2)/2
        rho[j,] = a_rho/b_rho
      }
      else
      {
        tmp = abs(mu[j,]-psi[j,])
        for ( i in 1:n )
          if ( tmp[i] < 1e-5 )
            rho[j,i] = b[j,i]*(1-tmp[i]/2+tmp[i]^2/6-tmp[i]^3/24)/2/(1+exp(-tmp[i]))
          else
            rho[j,i] = b[j,i]*(1-exp(-tmp[i]))/2/tmp[i]/(1+exp(-tmp[i]))
      }
    }
    
    
    # E-step for gamma
    for ( l in 1:L )
    {
      dir = rep(0,p)
      for ( j in 1:p )
        if ( Eidx[j] < Eidx[j+1] )
          if ( smoothing == "MRF" )
            dir[j] = eta*sum(thetaW[E[Eidx[j]:(Eidx[j+1]-1),2],l])
          else
            dir[j] = eta*sum(2*thetaW[E[Eidx[j]:(Eidx[j+1]-1),2],l]-1)
      dir = dir - lam - dlv/2 - W[,l]^2*div/2 - alpha[,l]
      grad = thetaW[,l]*(1-thetaW[,l])*dir
      dg = sum(dir*grad)
      
      stheta = sum(thetaW[,l])
      if ( smoothing == "MRF" )
        Q = -stheta*dlv/2 - sum(W[,l]^2*thetaW[,l])*div/2 - lam*stheta + eta*sum(thetaW[E[,1],l]*thetaW[E[,2],l])/2
      else
        Q = -stheta*dlv/2 - sum(W[,l]^2*thetaW[,l])*div/2 - lam*stheta + eta*sum(2*thetaW[E[,1],l]*thetaW[E[,2],l]-thetaW[E[,1],l]-thetaW[E[,2],l]+1)/2
      s = 1
      while (T)
      {
        newalpha = alpha[,l] + s*dir
        newtheta = 1/(1+exp(-newalpha))
        newstheta = sum(newtheta)
        if ( smoothing == "MRF" )
          newQ = -newstheta*dlv/2 - sum(W[,l]^2*newtheta)*div/2 - lam*newstheta + eta*sum(newtheta[E[,1]]*newtheta[E[,2]])/2
        else
          newQ = -newstheta*dlv/2 - sum(W[,l]^2*newtheta)*div/2 - lam*newstheta + eta*sum(2*newtheta[E[,1]]*newtheta[E[,2]]-newtheta[E[,1]]-newtheta[E[,2]]+1)/2
        if ( newQ-Q >= s*dg/2 )
        {
          alpha[,l] = newalpha
          thetaW[,l] = newtheta
          break
        }
        s = s / 2
      }
    }

    if ( GBC )
    {
      # E-step for delta
      thetaZ = 1/(1+sqrt(u1/u0)*exp(Z^2*diu/2+zeta))
    }    
    
    # M-step for Z
    for ( i in 1:n )
    {
      D = rho[,i]
      c = kappa[,i] + D*(psi[,i]-m)
      if ( GBC )
        iU = iu0+thetaZ[,i]*diu
      else
        iU = rep(1,L)
      chizvar = chol( diag(iU,L) + t(W)%*%(W*D) )
      Z[,i] = backsolve(chizvar,forwardsolve(t(chizvar),t(W)%*%c))
    }
    
    
    # M-step for m
    if ( intercept )
      for ( j in 1:p )
        m[j] = (-sum(W[j,]*(Z%*%rho[j,])) + sum(psi[j,]*rho[j,]) + sum(kappa[j,])) / sum(rho[j,])

    if ( PX )
    {
      # M-step for A
      A = Z%*%t(Z)/n
      
      # Rotation
      cA = t(chol(A))
      Z = forwardsolve(cA,Z)
    }
    
    # M-step for W
    for ( j in 1:p )
    {
      G = rho[j,]
      f = kappa[j,] + G*(psi[j,]-m[j])
      iV = iv0 + thetaW[j,]*div
      chiwvar = chol( diag(iV,L) + Z%*%(t(Z)*G) )
      W[j,] = backsolve(chiwvar,forwardsolve(t(chiwvar),Z%*%f))
    }

    
    # W and Z adjustments
    for ( l in 1:L )
    {
      ww = sum(W[,l]^2)
      zz = sum(Z[l,]^2)
      W[,l] = sqrt(sqrt(zz/ww)) * W[,l]
      Z[l,] = sqrt(sqrt(ww/zz)) * Z[l,]
    }
    
#    print(paste("Iteration ",iter," completed."))

    if ( max(abs(ptheta-thetaW)) < 1e-3 )
      break
  }
  
  ans = list(p=p,n=n,L=L,e=e,v0=v0,v1=v1,lam=lam,eta=eta,param=param,m=m,W=W,Z=Z,rho=rho,alpha=alpha,thetaW=thetaW,iter=iter)
  if ( GBC )
    ans = c(ans,list(u0=u0,u1=u1,thetaZ=thetaZ))
  ans
}





# function GFA_MCMC: runs MCMC for GFA
#
# Parameters:
# s: number of MCMC samples
# X: p by n data matrix where p is the number of covariates and n is the number of samples
# type: p by 1 vector of data types (see below)
# E: e by 2 matrix for edges where e is the number of edges.
#    E must be sorted by the first column and then the second column.
#    Both (j,k) and (k,j) must be added and therefore e is twice the actual number of edges.
# L: number of factors
# v0,v1: variances of spike and slab
# lam: shrinkage parameter
# eta: smoothing parameter
# bi: number of burn-ins
# param: p by 1 vector of parameters nu, n, r, or N
#
# Data types:
# 0: continuous Gaussian - nu is required as param
# 1: Binomial - n_j must be provided as param
# 2: Negative Binomial - r_j must be provided as param
# 3: Poisson - N must be provided as param

#library(BayesLogit)

GFA_MCMC <- function(s,X,type,E,L,v0,v1,lam,eta,bi,param)
{
  p = nrow(X)
  n = ncol(X)
  
  # edge
  e = nrow(E)
  Eidx = rep(1,p+1)
  for ( j in 1:p )
  {
    Eidx[j+1] = Eidx[j]
    while ( Eidx[j+1] <= e )
    {
      if ( E[Eidx[j+1],1] == j )
        Eidx[j+1] = Eidx[j+1] + 1
      else
        break
    }
  }
  
  # normalization
  
  # initialization
  W = matrix(0,p,L)
  Z = matrix(rnorm(L*n),L,n)
  m = rep(0,p)
  rho = matrix(0,p,n)
  gam = matrix(0,p,L)
  
  psi = matrix(0,p,n)
  kappa = matrix(0,p,n)
  b = matrix(0,p,n)
  
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
    {
      psi[j,] = X[j,]
    }
    if ( type[j] == 1 )
    {
      kappa[j,] = X[j,]-param[j]/2
      b[j,] = param[j]
    }
    if ( type[j] == 2 )
    {
      kappa[j,] = (X[j,]-param[j])/2
      b[j,] = param[j]+X[j,]
    }
    if ( type[j] == 3 )
    {
      psi[j,] = log(param[j])
      kappa[j,] = X[j,]-param[j]/2
      b[j,] = param[j]
    }
  }
  
  # containers
  Ws = array(0,c(p,L,s))
  Zs = array(0,c(L,n,s))
  ms = matrix(0,p,s)
  Gams = array(0,c(p,L,s))
  rhos = array(0,c(p,n,s))
  
  cnt = 0
  while ( cnt < s+bi )
  {
    # Sample rho
    mu = m + W%*%Z
    for ( j in 1:p )
    {
      if ( type[j] == 0 )
      {
        a_rho = (param[j]+n)/2
        b_rho = param[j]/2 + sum((mu[j,]-psi[j,])^2)/2
        rho[j,] = rgamma(1,a_rho,b_rho)
      }
      else
      {
        rho[j,] = rpg(n,b[j,],abs(mu[j,]-psi[j,]))
      }
    }
    
    # Sample Z
    for ( i in 1:n )
    {
      D = rho[,i]
      c = kappa[,i] + D*(psi[,i]-m)
      chizvar = chol( diag(1,L) + t(W)%*%(W*D) )
      chzvar = backsolve(chizvar,diag(1,L))
      zmean = chzvar %*% (t(chzvar) %*% (t(W) %*% c))
      Z[,i] = zmean + chzvar %*% rnorm(L)
    }
    
    # Sample W
    for ( j in 1:p )
    {
      G = rho[j,]
      f = kappa[j,] + G*(psi[j,]-m[j])
      iV = 1/v0 + gam[j,]*(1/v1-1/v0)
      chiwvar = chol( diag(iV,L) + Z%*%(t(Z)*G) )
      chwvar = backsolve(chiwvar,diag(1,L))
      wmean = chwvar %*% (t(chwvar) %*% (Z %*% f))
      W[j,] = wmean + chwvar %*% rnorm(L)
    }
    
    # Sample m
    for ( j in 1:p )
      m[j] = rnorm(1,(-sum(W[j,]*Z%*%rho[j,])+sum(psi[j,]*rho[j,])+sum(kappa[j,]))/sum(rho[j,]),sqrt(1/sum(rho[j,])))
    
    # Sample gamma
    for ( j in 1:p )
    {
      lp0 = -log(v0)/2 - W[j,]^2/v0/2
      lp1 = -log(v1)/2 - W[j,]^2/v1/2 - lam
      if ( Eidx[j] < Eidx[j+1] )
        for ( k in Eidx[j]:(Eidx[j+1]-1) )
          lp1 = lp1 + eta*gam[E[k,2],]
      gam[j,] = rbinom(L,1,1/(1+exp(lp0-lp1)))
    }
    
    cnt = cnt + 1
    if ( cnt > bi )
    {
      Ws[,,cnt-bi] = W
      Zs[,,cnt-bi] = Z
      ms[,cnt-bi] = m
      Gams[,,cnt-bi] = gam
      rhos[,,cnt-bi] = rho
    }
  }
  
  list(p=p,n=n,L=L,e=e,v0=v0,v1=v1,lam=lam,eta=eta,param=param,W=Ws,Z=Zs,m=ms,gam=Gams,rho=rhos)
}