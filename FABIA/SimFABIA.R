# Simulation GBC
#


if ( !require(fabia) )
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("fabia")
  library(fabia)
}

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


SimFABIA_BCV <- function(R,seed,p,n,type,param,overlap,L,thrW,thrZ,fold=3,batch=0)
{
  D1 = length(thrW)
  D2 = length(thrZ)
  
  BCV = array(0,c(D1,D2,fold,fold,R))
  opt_thrW = rep(0,R)
  opt_thrZ = rep(0,R)

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
    data = SimData(seed+batch+r,overlap,sigma2,p)
    
    n = ncol(data$Z)
    
    Wfold = fold(p,fold)
    Zfold = fold(n,fold)

    LL = nrow(data$Z)
    S[[r]] = list()
    for ( l in 1:LL )
      S[[r]][[l]] = list(r=which(data$W[,l]!=0),c=which(data$Z[l,]!=0))
    
    for ( s in 1:fold )
      for ( t in 1:fold )
      {
        Widx = which(Wfold==s)
        Zidx = which(Zfold==t)
        XA = data$X[Widx,Zidx]
        XB = data$X[Widx,-Zidx]
        XC = data$X[-Widx,Zidx]
        XD = data$X[-Widx,-Zidx]
        pD = p - length(Widx)
        nD = n - length(Zidx)

        fit = fabia(XD,L,random=0,center=2)
        
        for ( d1 in 1:D1 )
          for ( d2 in 1:D2 )
          {
            fabia_bic = extractBic(fit,thrW[d1],thrZ[d2])
            What = matrix(0,pD,L)
            Zhat = matrix(0,L,nD)
            idx = rep(T,L)
            for ( l in 1:L )
            {
              Widx = fabia_bic$numn[l,1]$numng
              Zidx = fabia_bic$numn[l,2]$numnp
              if ( length(Widx) != 0 & length(Zidx) != 0 )
              {
                What[Widx,l] = fit@L[Widx,l]
                Zhat[l,Zidx] = fit@Z[l,Zidx]
              }
              else
                idx[l] = F
            }
            if ( sum(idx) > 0 )
            {
              What = What[,idx]
              Zhat = Zhat[idx,]
              WWhat = t(What)%*%What
              ZZhat = Zhat%*%t(Zhat)
              err = XA - (XB%*%t(Zhat))%*%(chol2inv(chol(ZZhat))%*%chol2inv(chol(WWhat)))%*%(t(What)%*%XC)
            }
            else
              err = XA
            BCV[d1,d2,s,t,r] = sum(err^2)
          }
      }
    
    opt = which.min(apply(BCV[,,,,r],c(1,2),sum))
    opt_thrW[r] = thrW[(opt-1)%%D1+1]
    opt_thrZ[r] = thrZ[(opt-1)%/%D1+1]
    
    fit = fabia(data$X,L,random=0,center=2)
    fabia_bic = extractBic(fit,opt_thrW[r],opt_thrZ[r])
    
    Shat = list()
    for ( l in 1:L )
      Shat[[l]] =  list(r=fabia_bic$numn[l,1]$numng,c=fabia_bic$numn[l,2]$numnp )

    eval = gbcmetric(Shat,S[[r]],p,n)

    bclus[[r]] = Shat
    CE[r] = eval$CE
    FP[r] = eval$FP_CE
    FN[r] = eval$FN_CE
    SEN[r] = eval$SEN_CE
    SPE[r] = eval$SPE_CE
    MCC[r] = eval$MCC_CE
  }
  
  list(overlap=overlap,sigma2=sigma2,L=L,thrW=thrW,thrZ=thrZ,S=S,Shat=bclus,CE=CE,FP=FP,FN=FN,SEN=SEN,SPE=SPE,MCC=MCC,BCV=BCV,opt_thrW=opt_thrW,opt_thrZ=opt_thrZ)
}



SimFABIA_CCV <- function(R,seed,overlap,sigma2,L,thrW,thrZ,p=1000,fold=3,batch=0)
{
  D1 = length(thrW)
  D2 = length(thrZ)
  
  CCV = array(0,c(D1,D2,fold,fold,R))
  opt_thrW = rep(0,R)
  opt_thrZ = rep(0,R)
  
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
    data = SimData(seed+batch+r,overlap,sigma2,p)
    
    n = ncol(data$Z)
    
    Wfold = fold(p,fold)
    Zfold = fold(n,fold)
    
    LL = nrow(data$Z)
    S[[r]] = list()
    for ( l in 1:LL )
      S[[r]][[l]] = list(r=which(data$W[,l]!=0),c=which(data$Z[l,]!=0))
    
    for ( s in 1:fold )
    {
      Widx = which(Wfold==s)
      XA = data$X[Widx,]
      XB = data$X[-Widx,]
      pB = p - length(Widx)
      
      fit = fabia(XB,L,random=0,center=2)
      
      for ( d1 in 1:D1 )
        for ( d2 in 1:D2 )
        {
          fabia_bic = extractBic(fit,thrW[d1],thrZ[d2])
          What = matrix(0,pB,L)
          Zhat = matrix(0,L,n)
          idx = rep(T,L)
          for ( l in 1:L )
          {
            Widx = fabia_bic$numn[l,1]$numng
            Zidx = fabia_bic$numn[l,2]$numnp
            if ( length(Widx) != 0 & length(Zidx) != 0 )
            {
              What[Widx,l] = fit@L[Widx,l]
              Zhat[l,Zidx] = fit@Z[l,Zidx]
            }
            else
              idx[l] = F
          }
          if ( sum(idx) > 0 )
          {
            Zhat = Zhat[idx,]
            ZZhat = Zhat%*%t(Zhat)
            XAZhat = XA%*%t(Zhat)
            CCV[d1,d2,s,1,r] = sum(XA^2) - sum(XAZhat*(XAZhat%*%chol2inv(chol(ZZhat))))
          }
          else
            CCV[d1,d2,s,1,r] = sum(XA^2)
        }
      
      Zidx = which(Zfold==s)
      XA = data$X[,Zidx]
      XB = data$X[,-Zidx]
      nB = n - length(Zidx)
      
      fit = fabia(XB,L,random=0,center=2)
      
      for ( d1 in 1:D1 )
        for ( d2 in 1:D2 )
        {
          fabia_bic = extractBic(fit,thrW[d1],thrZ[d2])
          What = matrix(0,p,L)
          Zhat = matrix(0,L,nB)
          idx = rep(T,L)
          for ( l in 1:L )
          {
            Widx = fabia_bic$numn[l,1]$numng
            Zidx = fabia_bic$numn[l,2]$numnp
            if ( length(Widx) != 0 & length(Zidx) != 0 )
            {
              What[Widx,l] = fit@L[Widx,l]
              Zhat[l,Zidx] = fit@Z[l,Zidx]
            }
            else
              idx[l] = F
          }
          if ( sum(idx) > 0 )
          {
            What = What[,idx]
            WWhat = t(What)%*%What
            XAWhat = t(XA)%*%What
            CCV[d1,d2,s,2,r] = sum(XA^2) - sum(XAWhat*(XAWhat%*%chol2inv(chol(WWhat))))
          }
          else
            CCV[d1,d2,s,2,r] = sum(XA^2)
        }
    }
    
    opt = which.min(apply(CCV[,,,,r],c(1,2),sum))
    opt_thrW[r] = thrW[(opt-1)%%D1+1]
    opt_thrZ[r] = thrZ[(opt-1)%/%D1+1]
    
    fit = fabia(data$X,L,random=0,center=2)
    fabia_bic = extractBic(fit,opt_thrW[r],opt_thrZ[r])
    
    Shat = list()
    for ( l in 1:L )
      Shat[[l]] =  list(r=fabia_bic$numn[l,1]$numng,c=fabia_bic$numn[l,2]$numnp )
    
    eval = gbcmetric(Shat,S[[r]],p,n)

    bclus[[r]] = Shat
    CE[r] = eval$CE
    FP[r] = eval$FP_CE
    FN[r] = eval$FN_CE
    SEN[r] = eval$SEN_CE
    SPE[r] = eval$SPE_CE
    MCC[r] = eval$MCC_CE
  }
  
  list(overlap=overlap,sigma2=sigma2,L=L,thrW=thrW,thrZ=thrZ,S=S,Shat=bclus,CE=CE,FP=FP,FN=FN,SEN=SEN,SPE=SPE,MCC=MCC,CCV=CCV,opt_thrW=opt_thrW,opt_thrZ=opt_thrZ)
}



SimFABIA_BIC <- function(R,seed,p,n,type,param,overlap,L,thrW,thrZ,batch=0)
{
  D1 = length(thrW)
  D2 = length(thrZ)
  
  BIC = array(0,c(D1,D2,R))
  opt_BIC = rep(Inf,R)
  opt_thrW = rep(0,R)
  opt_thrZ = rep(0,R)
  
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
    data = SimData(seed+batch+r,overlap,sigma2,p)
    
    LL = nrow(data$Z)
    S[[r]] = list()
    for ( l in 1:LL )
      S[[r]][[l]] = list(r=which(data$W[,l]!=0),c=which(data$Z[l,]!=0))
    
    fit = fabia(data$X,L,random=0,center=2)
    fit.m = apply(data$X,1,median)
    
    for ( d1 in 1:D1 )
      for ( d2 in 1:D2 )
      {
        fabia_bic = extractBic(fit,thrW[d1],thrZ[d2])
        
        nParam = p
        What = matrix(0,p,L)
        Zhat = matrix(0,L,n)
        
        for ( l in 1:L )
        {
          Widx = fabia_bic$numn[l,1]$numng
          Zidx = fabia_bic$numn[l,2]$numnp
          if ( length(Widx) != 0 & length(Zidx) != 0 )
          {
            nParam = nParam + length(Widx) + length(Zidx)
            What[Widx,l] = fit@L[Widx,l]
            Zhat[l,Zidx] = fit@Z[l,Zidx]
          }
        }
        
        muhat = What %*% Zhat + fit.m
        
        BIC[d1,d2,r] = -2*llk(data$X,muhat,data$type,data$param) + nParam*log(n*p)
        
        if ( opt_BIC[r] > BIC[d1,d2,r] )
        {
          opt_BIC[r] = BIC[d1,d2,r]
          opt_thrW[r] = thrW[d1]
          opt_thrZ[r] = thrZ[d2]
        }
        
      }
    
    
    fabia_bic = extractBic(fit,opt_thrW[r],opt_thrZ[r])
    
    Shat = list()
    for ( l in 1:L )
      Shat[[l]] =  list(r=fabia_bic$numn[l,1]$numng,c=fabia_bic$numn[l,2]$numnp )
    
    eval = gbcmetric(Shat,S[[r]],p,n)
    
    bclus[[r]] = Shat
    CE[r] = eval$CE
    FP[r] = eval$FP_CE
    FN[r] = eval$FN_CE
    SEN[r] = eval$SEN_CE
    SPE[r] = eval$SPE_CE
    MCC[r] = eval$MCC_CE
  }
  
  list(p=p,n=n,type=type,param=param,overlap=overlap,L=L,thrW=thrW,thrZ=thrZ,S=S,CE=CE,FP=FP,FN=FN,SEN=SEN,SPE=SPE,MCC=MCC,BIC=BIC,opt_thrW=opt_thrW,opt_thrZ=opt_thrZ)
}



SimFABIA_Plain <- function(R,seed,p,n,type,param,overlap,L,thrW,thrZ,batch=0)
{
  D1 = length(thrW)
  D2 = length(thrZ)
  
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
    data = SimData(seed+batch+r,p,n,type,param,overlap)
    
    LL = nrow(data$Z)
    S[[r]] = list()
    for ( l in 1:LL )
      S[[r]][[l]] = list(r=which(data$W[,l]!=0),c=which(data$Z[l,]!=0))
    
    fit = fabia(data$X,L,random=0,center=2)

    for ( d1 in 1:D1 )
      for ( d2 in 1:D2 )
      {
        fabia_bic = extractBic(fit,thrW[d1],thrZ[d2])
        
        Shat = list()
        for ( l in 1:L )
          Shat[[l]] =  list(r=fabia_bic$numn[l,1]$numng,c=fabia_bic$numn[l,2]$numnp )
        
        eval = gbcmetric(Shat,S[[r]],p,n)
    
        CE[d1,d2,r] = eval$CE
        FP[d1,d2,r] = eval$FP_CE
        FN[d1,d2,r] = eval$FN_CE
        SEN[d1,d2,r] = eval$SEN_CE
        SPE[d1,d2,r] = eval$SPE_CE
        MCC[d1,d2,r] = eval$MCC_CE
      }
    
    fits[[r]] = fit
  }
  
  list(p=p,n=n,type=type,param=param,overlap=overlap,L=L,thrW=thrW,thrZ=thrZ,S=S,CE=CE,FP=FP,FN=FN,SEN=SEN,SPE=SPE,MCC=MCC)
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
