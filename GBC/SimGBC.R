# Simulation GBC
#


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

if ( file.exists("Sim/SimData.R") )
  source("Sim/SimData.R")
if ( file.exists("/home/changgee/project/GBC/Sim/SimData.R") )
  source("/home/changgee/project/GBC/Sim/SimData.R")
if ( file.exists("/home/cchan40/project/GBC/Sim/SimData.R") )
  source("/home/cchan40/project/GBC/Sim/SimData.R")

  

SimGBC_BCV <- function(R,seed,p,n,type,param,overlap,L,k,v0,lam,eta,intercept=T,smoothing="MRF",thres=0.5,fold=3,batch=0)
{
  D1 = length(v0)
  D2 = length(lam)
  
  BCV = array(0,c(D1,D2,fold,fold,R))
  opt_BCV = rep(Inf,R)
  opt_v0 = rep(0,R)
  opt_lam = rep(0,R)
  
  bclus = list()
  CE = rep(0,R)
  FP = rep(0,R)
  FN = rep(0,R)
  SEN = rep(0,R)
  SPE = rep(0,R)
  MCC = rep(0,R)

  S = list()
  
  for ( r in 1:R )
  {
    data = SimData(seed+batch+r,p,n,type,param,overlap)

    LL = nrow(data$Z)
    S[[r]] = list()
    for ( l in 1:LL )
      S[[r]][[l]] = list(r=which(data$W[,l]!=0),c=which(data$Z[l,]!=0))
    
    set.seed(seed+batch+r)
    Wfold = fold(p,fold)
    Zfold = fold(n,fold)
    
    for ( d1 in 1:D1 )
      for ( d2 in 1:D2 )
      {
        for ( s in 1:fold )
          for ( t in 1:fold )
          {
            Widx = which(Wfold==s)
            Zidx = which(Zfold==t)
            XA = data$X[Widx,Zidx]
            XB = data$X[Widx,-Zidx]
            XC = data$X[-Widx,Zidx]
            XD = data$X[-Widx,-Zidx]
            typeD = data$type[-Widx]
            paramD = data$param[-Widx]
            
            # creating new edge information for XD
            pA = length(Widx)
            pD = p - pA
            nA = length(Zidx)
            nD = n - nA
            newidx = rep(0,p)
            newidx[-Widx] = 1:(pD)
            EDidx = Wfold[data$E[,1]]!=s & Wfold[data$E[,2]]!=s
            eD = sum(EDidx)
            ED = matrix(newidx[data$E[EDidx,]],eD)
        
            fit = GFA_EM(XD,typeD,ED,L,v0[d1],k*v0[d1],lam[d2],eta,paramD,intercept,smoothing,GBC=T)

            WDhat = fit$W*(fit$thetaW>thres)
            ZDhat = fit$Z*(fit$thetaZ>thres)
            idx = apply(fit$thetaW>thres,2,sum)>0 & apply(fit$thetaZ>thres,1,sum)>0
            
            if ( sum(idx) > 0 )
            {
              What = matrix(What[,idx],nrow=pD)
              Zhat = matrix(Zhat[idx,],ncol=nD)
              WWhat = t(What)%*%What
              ZZhat = Zhat%*%t(Zhat)
              muAhat = (XB%*%t(Zhat))%*%(chol2inv(chol(ZZhat))%*%chol2inv(chol(WWhat)))%*%(t(What)%*%XC)
            }
            else
              muAhat = matrix(0,pA,nA)
            BCV[d1,d2,s,t,r] = sum(err^2)
          }
      }
    
  }
  
  list(p=p,n=n,type=type,param=param,overlap=overlap,L=L,k=k,v0=v0,lam=lam,eta=eta,BCV=BCV)
}



SimGBC_CCV <- function(R,seed,p,n,type,param,overlap,L,k,v0,lam,eta,intercept=T,smoothing="MRF",thres=0.5,fold=3,batch=0)
{
  D1 = length(v0)
  D2 = length(lam)
  
  CCV = array(0,c(D1,D2,fold,2,R))
  opt_v0 = rep(0,R)
  opt_lam = rep(0,R)
  
  bclus = list()
  CE = rep(0,R)
  FP = rep(0,R)
  FN = rep(0,R)
  SEN = rep(0,R)
  SPE = rep(0,R)
  MCC = rep(0,R)
  
  S = list()
  
  for ( r in 1:R )
  {
    data = SimData(seed+batch+r,p,n,type,param,overlap)
    
    set.seed(seed+batch+r)
    Wfold = fold(p,fold)
    Zfold = fold(n,fold)
    
    for ( s in 1:fold )
    {
      Widx = which(Wfold==s)
      XA = data$X[Widx,]
      XB = data$X[-Widx,]
      typeB = data$type[-Widx]
      paramB = data$param[-Widx]
      
      # creating new edge information for XB
      newidx = rep(0,p)
      pA = length(Widx)
      pB = p - pA
      newidx[-Widx] = 1:(pB)
      EBidx = Wfold[data$E[,1]]!=s & Wfold[data$E[,2]]!=s
      eB = sum(EBidx)
      EB = matrix(newidx[data$E[EBidx,]],eB,2)
      
      for ( d1 in 1:D1 )
        for ( d2 in 1:D2 )
        {
          fit = GFA_EM(XB,typeB,EB,L,v0[d1],k*v0[d1],lam[d2],eta,paramB,intercept,smoothing,GBC=T)
          idx = apply(fit$thetaW>thres,2,sum)>0 & apply(fit$thetaZ>thres,1,sum)>0
          Zhat = fit$Z*(fit$thetaZ>thres)
          Zhat = matrix(Zhat[idx,],ncol=n)
          if ( intercept )
            Zhat = rbind(1,Zhat)
          ZZhat = Zhat%*%t(Zhat)
          XAZhat = XA%*%t(Zhat)
          CCV[d1,d2,s,1,r] = sum(XA^2) - sum(XAZhat*(XAZhat%*%chol2inv(chol(ZZhat))))
        }

      Zidx = which(Zfold==s)
      XA = data$X[,Zidx]
      XB = data$X[,-Zidx]
      for ( d1 in 1:D1 )
        for ( d2 in 1:D2 )
        {
          fit = GFA_EM(XB,type,data$E,L,v0[d1],k*v0[d1],lam[d2],eta,data$param,intercept,smoothing,GBC=T)
          idx = apply(fit$thetaW>thres,2,sum)>0 & apply(fit$thetaZ>thres,1,sum)>0
          What = fit$W*(fit$thetaW>thres)
          What = matrix(What[,idx],nrow=p)
          if ( intercept )
            XA = XA - fit$m
          WWhat = t(What)%*%What
          XAWhat = t(XA)%*%What
          CCV[d1,d2,s,2,r] = sum(XA^2) - sum(XAWhat*(XAWhat%*%chol2inv(chol(WWhat))))
        }
      
    } 
  }
  
  list(p=p,n=n,type=type,param=param,overlap=overlap,L=L,k=k,v0=v0,lam=lam,eta=eta,CCV=CCV)
}


SimGBC_BIC <- function(R,seed,p,n,type,param,overlap,L,k,v0,lam,eta,center=1,smoothing="MRF",thres=0.5,batch=0)
{
  D1 = length(v0)
  D2 = length(lam)
  
  BIC = array(0,c(D1,D2,R))
  opt_BIC = rep(Inf,R)
  opt_v0 = rep(0,R)
  opt_lam = rep(0,R)
  
  bclus = list()
  CE = rep(0,R)
  FP = rep(0,R)
  FN = rep(0,R)
  SEN = rep(0,R)
  SPE = rep(0,R)
  MCC = rep(0,R)

  S = list()

  for ( r in 1:R )
  {
    data = SimData(seed+batch+r,p,n,type,param,overlap)
    
    LL = nrow(data$Z)
    S[[r]] = list()
    for ( l in 1:LL )
      S[[r]][[l]] = list(r=which(data$W[,l]!=0),c=which(data$Z[l,]!=0))

    for ( d1 in 1:D1 )
      for ( d2 in 1:D2 )
      {
        fit = GFA_EM(data$X,data$type,data$E,L,v0[d1],k*v0[d1],lam[d2],eta,data$param,intercept,smoothing,GBC=T)

        What = fit$W*(fit$thetaW>thres)
        Zhat = fit$Z*(fit$thetaZ>thres)
        muhat = What %*% Zhat + fit$m

        nParam = sum(fit$thetaW>thres) + sum(fit$thetaZ>thres)
        if ( center != 0 )
          nParam = nParam + p

        BIC[d1,d2,r] = -2*llk(data$X,muhat,data$type,data$param) + nParam*log(n*p)
        
        if ( opt_BIC[r] > BIC[d1,d2,r] )
        {
          opt_BIC[r] = BIC[d1,d2,r]
          opt_v0[r] = v0[d1]
          opt_lam[r] = lam[d2]
          opt_fit = fit
        }
      }

    Shat = list()
    for ( l in 1:L )
      Shat[[l]] = list(r=which(opt_fit$thetaW[,l]>thres),c=which(opt_fit$thetaZ[l,]>thres))
    
    eval = gbcmetric(Shat,S[[r]],p,n)
    
    bclus[[r]] = Shat
    CE[r] = eval$CE
    FP[r] = eval$FP_CE
    FN[r] = eval$FN_CE
    SEN[r] = eval$SEN_CE
    SPE[r] = eval$SPE_CE
    MCC[r] = eval$MCC_CE
  }
  
  list(p=p,n=n,type=type,param=param,overlap=overlap,L=L,k=k,v0=v0,lam=lam,eta=eta,S=S,Shat=bclus,CE=CE,FP=FP,FN=FN,SEN=SEN,SPE=SPE,MCC=MCC,BIC=BIC,opt_v0=opt_v0,opt_lam=opt_lam)
  
}


SimGBC_Plain <- function(R,seed,p,n,type,param,overlap,L,k,v0,lam,eta,center=1,smoothing="Ising",thres=0.5,batch=0)
{
  D1 = length(v0)
  D2 = length(lam)
  
  CE = array(0,c(D1,D2,R))
  FP = array(0,c(D1,D2,R))
  FN = array(0,c(D1,D2,R))
  SEN = array(0,c(D1,D2,R))
  SPE = array(0,c(D1,D2,R))
  MCC = array(0,c(D1,D2,R))
  
  S = list()
  fits = list()
  
  for ( r in 1:R )
  {
    print(r)
    data = SimData(seed+batch+r,p,n,type,param,overlap)
    
    LL = nrow(data$Z)
    S[[r]] = list()
    for ( l in 1:LL )
      S[[r]][[l]] = list(r=which(data$W[,l]!=0),c=which(data$Z[l,]!=0))
    
    fits[[r]] = list()
    for ( d1 in 1:D1 )
    {
      fits[[r]][[d1]] = list()
      for ( d2 in 1:D2 )
      {
        fit = GFA_EM(data$X,data$type,data$E,L,v0[d1],k*v0[d1],lam[d2],eta,data$param,center,smoothing,GBC=T)
        
        Shat = list()
        for ( l in 1:L )
          Shat[[l]] = list(r=which(fit$thetaW[,l]>thres),c=which(fit$thetaZ[l,]>thres))
        eval = gbcmetric(Shat,S[[r]],p,n)
        
        fits[[r]][[d1]][[d2]] = fit
        CE[d1,d2,r] = eval$CE
        FP[d1,d2,r] = eval$FP_CE
        FN[d1,d2,r] = eval$FN_CE
        SEN[d1,d2,r] = eval$SEN_CE
        SPE[d1,d2,r] = eval$SPE_CE
        MCC[d1,d2,r] = eval$MCC_CE
      }
    }
  }
  
  list(p=p,n=n,type=type,param=param,overlap=overlap,L=L,k=k,v0=v0,lam=lam,eta=eta,S=S,CE=CE,FP=FP,FN=FN,SEN=SEN,SPE=SPE,MCC=MCC)
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


tuneparam <- function(k,b,lam)
{
  a = b^2*(k-1)/2/k/sqrt(2*log(k))
  a/sqrt(lam)
}



beta_cross <- function(v0,v1,lam)
{
  sqrt(v0*v1*(log(v1/v0)+2*lam)/(v1-v0))
}



