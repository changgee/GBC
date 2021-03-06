# GBC for Data Analysis


if ( file.exists("GBC/GFA.R") )
  source("GBC/GFA.R")
if ( file.exists("/home/changgee/project/GBC/GBC/GFA.R") )
  source("/home/changgee/project/GBC/GBC/GFA.R")
if ( file.exists("/home/cchan40/project/GBC/GBC/GFA.R") )
  source("/home/cchan40/project/GBC/GBC/GFA.R")

if ( file.exists("Eval/gbcmetric.R") )
  source("Eval/gbcmetric.R")
if ( file.exists("/home/changgee/project/GBC/Eval/gbcmetric.R") )
  source("/home/changgee/project/GBC/Eval/gbcmetric.R")
if ( file.exists("/home/cchan40/project/GBC/Eval/gbcmetric.R") )
  source("/home/cchan40/project/GBC/Eval/gbcmetric.R")


DataGBC_BCV <- function(X,E,L,k,v0,lam,eta,param,intercept=F,smoothing="MRF",thres=0.5,fold=3,seed=100,run=NULL)
{
  p = nrow(X)
  n = ncol(X)
  D1 = length(v0)
  D2 = length(lam)
  
  BCV = array(0,c(D1,D2,fold,fold))
  
  type = rep(0,p)
  
  set.seed(seed)
  Wfold = fold(p,fold)
  Zfold = fold(n,fold)
  
  for ( s in 1:fold )
    for ( t in 1:fold )
    {
      Widx = which(Wfold==s)
      Zidx = which(Zfold==t)
      XA = X[Widx,Zidx]
      XB = X[Widx,-Zidx]
      XC = X[-Widx,Zidx]
      XD = X[-Widx,-Zidx]
      typeD = type[-Widx]
      
      # creating new edge information for XD
      newidx = rep(0,p)
      pA = length(Widx)
      pD = p - pA
      nA = length(Zidx)
      nD = n - nA
      newidx[-Widx] = 1:(pD)
      EDidx = Wfold[E[,1]]!=s & Wfold[E[,2]]!=s
      eD = sum(EDidx)
      ED = matrix(newidx[E[EDidx,]],eD,2)
      
      for ( d1 in 1:D1 )
        for ( d2 in 1:D2 )
        {
          fit = GFA_EM(XD,typeD,ED,L,v0[d1],k*v0[d1],lam[d2],eta,param,intercept,smoothing,GBC=T)
          idx = apply(fit$thetaW>thres,2,sum)>0 & apply(fit$thetaZ>thres,1,sum)>0
          What = fit$W*(fit$thetaW>thres)
          What = matrix(What[,idx],nrow=pD)
          Zhat = fit$Z*(fit$thetaZ>thres)
          Zhat = matrix(Zhat[idx,],ncol=nD)
          if ( intercept )
          {
            What = cbind(fit$m,What)
            Zhat = rbind(1,Zhat)
          }
          WWhat = t(What)%*%What
          ZZhat = Zhat%*%t(Zhat)
          err = XA - (XB%*%t(Zhat))%*%(chol2inv(chol(ZZhat))%*%chol2inv(chol(WWhat)))%*%(t(What)%*%XC)
          BCV[d1,d2,s,t] = sum(err^2)
        }
    }

  list(L=L,k=k,v0=v0,lam=lam,eta=eta,fold=fold,BCV=BCV)
}  


DataGBC_CCV <- function(X,E,L,k,v0,lam,eta,param,intercept=F,smoothing="Ising",thres=0.5,fold=3,seed=100,run=NULL)
{
  p = nrow(X)
  n = ncol(X)
  D1 = length(v0)
  D2 = length(lam)
  
  CCV = array(0,c(D1,D2,fold,2))
  
  type = rep(0,p)
  
  set.seed(seed)
  Wfold = fold(p,fold)
  Zfold = fold(n,fold)
  
  for ( s in 1:fold )
  {
    if ( is.null(run) | (run[1]==1 & run[2]==s) )
    {
      Widx = which(Wfold==s)
      XA = X[Widx,]
      XB = X[-Widx,]
      typeB = type[-Widx]
      
      # creating new edge information for XB
      newidx = rep(0,p)
      pA = length(Widx)
      pB = p - pA
      newidx[-Widx] = 1:(pB)
      EBidx = Wfold[E[,1]]!=s & Wfold[E[,2]]!=s
      eB = sum(EBidx)
      EB = matrix(newidx[E[EBidx,]],eB,2)
      
      for ( d1 in 1:D1 )
        if ( is.null(run) | run[3]==d1 )
          for ( d2 in 1:D2 )
            if ( is.null(run) | run[4]==d2 )
            {
              fit = GFA_EM(XB,typeB,EB,L,v0[d1],k*v0[d1],lam[d2],eta,param,intercept,smoothing,GBC=T)
              idx = apply(fit$thetaW>thres,2,sum)>0 & apply(fit$thetaZ>thres,1,sum)>0
              Zhat = fit$Z*(fit$thetaZ>thres)
              Zhat = matrix(Zhat[idx,],ncol=n)
              if ( intercept )
                Zhat = rbind(1,Zhat)
              ZZhat = Zhat%*%t(Zhat)
              XAZhat = XA%*%t(Zhat)
              CCV[d1,d2,s,1] = sum(XA^2) - sum(XAZhat*(XAZhat%*%chol2inv(chol(ZZhat))))
            }
    }
      
    if ( is.null(run) | (run[1]==2 & run[2]==s) )
    {
      Zidx = which(Zfold==s)
      XA = X[,Zidx]
      XB = X[,-Zidx]
      for ( d1 in 1:D1 )
        if ( is.null(run) | run[3]==d1 )
          for ( d2 in 1:D2 )
            if ( is.null(run) | run[4]==d2 )
            {
              fit = GFA_EM(XB,type,E,L,v0[d1],k*v0[d1],lam[d2],eta,param,intercept,smoothing,GBC=T)
              idx = apply(fit$thetaW>thres,2,sum)>0 & apply(fit$thetaZ>thres,1,sum)>0
              What = fit$W*(fit$thetaW>thres)
              What = matrix(What[,idx],nrow=p)
              if ( intercept )
                XA = XA - fit$m
              WWhat = t(What)%*%What
              XAWhat = t(XA)%*%What
              CCV[d1,d2,s,2] = sum(XA^2) - sum(XAWhat*(XAWhat%*%chol2inv(chol(WWhat))))
            }
      
    }
  }
  
  if ( is.null(run) )
  {
    opt = which.min(apply(CCV,c(1,2),sum))
    opt_v0 = v0[(opt-1)%%D1+1]
    opt_lam = lam[(opt-1)%/%D1+1]
  }
  else
  {
    opt_v0 = v0[run[3]]
    opt_lam = lam[run[4]]
  }
  
  if ( is.null(run) | (run[1]==1 & run[2]==1) )
  {
    time = system.time(fit <- GFA_EM(X,type,E,L,opt_v0,k*opt_v0,opt_lam,eta,param,intercept,smoothing,GBC=T))
    Shat = list()
    for ( l in 1:L )
      Shat[[l]] = list(r=which(fit$thetaW[,l]>thres),c=which(fit$thetaZ[l,]>thres))
  }
  else
  {
    fit = NULL
    time = NULL
    Shat = list()
  }
  
  list(L=L,k=k,v0=v0,lam=lam,eta=eta,fold=fold,CCV=CCV,opt_v0=opt_v0,opt_lam=opt_lam,opt_fit=fit,opt_biclus=Shat,time=as.numeric(time[1]))
}  

DataGBC_BIC <- function(datapath,outpath,name,L,k,v0,lam,bias,eta,smoothing="Ising",thres=0.5)
{
  D1 = length(L)
  D2 = length(k)
  D3 = length(v0)
  D4 = length(lam)
  D5 = length(bias)
  
  BIC = array(0,c(D1,D2,D3,D4,D5))
  opt_BIC = Inf

  load(datapath)

  for ( d1 in 1:D1 )
    for ( d2 in 1:D2 )
      for ( d3 in 1:D3 )
        for ( d4 in 1:D4 )
          for ( d5 in 1:D5 )
          {
            fname = sprintf("res_%s_GBC_%02d_%02d_%.4f_%.4f_%.1f_%.2f",name,L[d1],k[d2],v0[d3],lam[d4],bias[d5],eta)
            fpath = paste(outpath,fname,sep="/")
            load(fpath)
            
            What = What*(thetaWhat>thres)
            Zhat = Zhat*(thetaZhat>thres)
            muhat = What %*% Zhat + mhat
            
            nParam = sum(thetaWhat>thres) + sum(thetaZhat>thres)
            nParam = nParam + p
            
            BIC[d1,d2,d3,d4,d5] = -2*llk(X,muhat,type,param) + nParam*log(n*p)
            
            if ( opt_BIC > BIC[d1,d2,d3,d4,d5] )
            {
              opt_BIC = BIC[d1,d2,d3,d4,d5]
              opt_L = L[d1]
              opt_k = k[d2]
              opt_v0 = v0[d3]
              opt_lam = lam[d4]
              opt_bias = bias[d5]
            }
          }

  list(name=name,L=L,k=k,v0=v0,lam=lam,eta=eta,BIC=BIC,opt_BIC=opt_BIC,opt_L=opt_L,opt_k=opt_k,opt_v0=opt_v0,opt_lam=opt_lam,opt_bias=opt_bias)
}  


DataGBC_RES <- function(datapath,outpath,name,L,k,v0,lam,bias,eta,smoothing="Ising",thres=0.5)
{
  D1 = length(L)
  D2 = length(k)
  D3 = length(v0)
  D4 = length(lam)
  D5 = length(bias)
  D6 = length(eta)
  
  CE = array(0,c(D1,D2,D3,D4,D5,D6))
  CS = array(0,c(D1,D2,D3,D4,D5,D6))
  nbc = array(0,c(D1,D2,D3,D4,D5,D6))
  ng = array(0,c(D1,D2,D3,D4,D5,D6))
  ns = array(0,c(D1,D2,D3,D4,D5,D6))
  BIC = array(0,c(D1,D2,D3,D4,D5,D6))

  load(datapath)
  S = list()
  for ( i in 1:9 )
    S[[i]] = list(r=1,c=which(org==i))

  for ( d1 in 1:D1 )
    for ( d2 in 1:D2 )
      for ( d3 in 1:D3 )
        for ( d4 in 1:D4 )
          for ( d5 in 1:D5 )
            for ( d6 in 1:D6 )
            {
              fname = sprintf("res_%s_GBC_%02d_%02d_%.5f_%.4f_%02d_%.1f",name,L[d1],k[d2],v0[d3],lam[d4],bias[d5],eta[d6])
              fpath = paste(outpath,fname,sep="/")
              load(fpath)
            
              Shat = list()
              for ( l in 1:L[d1] )
                Shat[[l]] = list(r=which(thetaWhat[,l]>thres),c=which(thetaZhat[l,]>thres))
    
              pf = gbcmetric(Shat,S,p,n,1)
              CE[d1,d2,d3,d4,d5,d6] = pf$CE
              CS[d1,d2,d3,d4,d5,d6] = pf$CS
              nbc[d1,d2,d3,d4,d5,d6] = pf$L1
              
              
              W = matrix(0,p,L[d1])
              Z = matrix(0,L[d1],n)
              
              nParam = p
              ng[d1,d2,d3,d4,d5,d6] = 0
              ns[d1,d2,d3,d4,d5,d6] = 0
              for ( l in 1:L[d1] )
              {
                Widx = which(thetaWhat[,l]>thres)
                Zidx = which(thetaZhat[l,]>thres)
                if ( length(Widx) != 0 & length(Zidx) != 0 )
                {
                  nParam = nParam + length(Widx) + length(Zidx)
                  ng[d1,d2,d3,d4,d5,d6] = ng[d1,d2,d3,d4,d5,d6] + length(Widx)
                  ns[d1,d2,d3,d4,d5,d6] = ns[d1,d2,d3,d4,d5,d6] + length(Zidx)
                  W[Widx,l] = What[Widx,l]
                  Z[l,Zidx] = Zhat[l,Zidx]
                }
              }
              ng[d1,d2,d3,d4,d5,d6] = ng[d1,d2,d3,d4,d5,d6] / pf$L1
              ns[d1,d2,d3,d4,d5,d6] = ns[d1,d2,d3,d4,d5,d6] / pf$L1
              muhat = W %*% Z + mhat
              
              BIC[d1,d2,d3,d4,d5,d6] = -2*llk(X,muhat,type,param) + nParam*log(n*p)
            }
  
  list(name=name,L=L,k=k,v0=v0,lam=lam,eta=eta,CE=CE,CS=CS,nbc=nbc,ng=ng,ns=ns,BIC=BIC)
}  

DataGBC_Plain <- function(datapath,outpath,name,L,k,v0,lam,bias,eta,center=1,smoothing="Ising")
{
  D1 = length(v0)
  D2 = length(lam)

  load(datapath)
  
  for ( d1 in 1:D1 )
    for ( d2 in 1:D2 )
    {
      fname = sprintf("res_%s_GBC_%02d_%02d_%.5f_%.4f_%02d_%.1f",name,L,k,v0[d1],lam[d2],bias,eta)
      fpath = paste(outpath,fname,sep="/")
      if ( !file.exists(fpath) )
      {
        time = system.time(fit <- GFA_EM(X,type,E,L,v0[d1],k*v0[d1],lam[d2],eta,param,center,smoothing=smoothing,GBC=T,zeta=lam[d2]+bias))
        
        What = fit$W
        Zhat = fit$Z
        mhat = fit$m
        thetaWhat = fit$thetaW
        thetaZhat = fit$thetaZ
        iter = fit$iter
        
        save(What,Zhat,mhat,thetaWhat,thetaZhat,iter,time,file=fpath)
      }
    }    

}  



# N: sample size
# K: # of folds
# r: repeats
fold <- function(N,K,r=1)
{
  Nk = rep(0,K)
  n = N
  for ( k in 1:K )
  {
    Nk[k] = ceiling(n/(K+1-k))
    n = n - Nk[k]
  }
  
  fold = matrix(0,N,r)
  for ( i in 1:r )
    fold[,i] = sample(rep(1:K,Nk),N)
  
  fold
}

