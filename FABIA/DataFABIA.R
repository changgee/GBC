library(fabia)

DataFABIA_BCV <- function(X,L,thW,thZ,fold=3,run=NULL)
{
  p = nrow(X)
  n = ncol(X)
  D1 = length(thW)
  D2 = length(thZ)
  
  BCV = array(0,c(D1,D2,fold,fold))
  
  Wfold = fold(p,fold)
  Zfold = fold(n,fold)
  
  for ( s in 1:fold )
    if ( is.null(run) | run[1]==s )
      for ( t in 1:fold )
        if ( is.null(run) | run[2]==t )
        {
          Widx = which(Wfold==s)
          Zidx = which(Zfold==t)
          XA = X[Widx,Zidx]
          XB = X[Widx,-Zidx]
          XC = X[-Widx,Zidx]
          XD = X[-Widx,-Zidx]
        
          fit = fabia(XD,L,random=0,center=0)

          for ( d1 in 1:D1 )
            for ( d2 in 1:D2 )
            {
              fabia_bic = extractBic(fit,thW[d1],thZ[d2])
              
              Shat = list()
              for ( l in 1:L )
                Shat[[l]] =  list(r=fabia_bic$numn[l,1][[1]],c=fabia_bic$numn[l,2][[1]] )
              eval = gbcmetric(Shat,S,p,n)
              
              idx = apply(fit$thetaW>thres,2,sum)>0 & apply(fit$thetaZ>thres,1,sum)>0
              What = fit$W*(fit$thetaW>thres)
              What = matrix(What[,idx],nrow=pD)
              Zhat = fit$Z*(fit$thetaZ>thres)
              Zhat = matrix(Zhat[idx,],ncol=nD)
              WWhat = t(What)%*%What
              ZZhat = Zhat%*%t(Zhat)
              err = XA - (XB%*%t(Zhat))%*%(chol2inv(chol(ZZhat))%*%chol2inv(chol(WWhat)))%*%(t(What)%*%XC)
              BCV[d1,d2,s,t] = sum(err^2)
            }
        }
  
  if ( is.null(run) )
  {
    opt = which.min(apply(BCV,c(1,2),sum))
    opt_v0 = v0[(opt-1)%%D1+1]
    opt_lam = lam[(opt-1)%/%D1+1]
  }
  else
  {
    opt_v0 = v0[run[2]]
    opt_lam = lam[run[3]]
  }
  
  if ( is.null(run) | run[1]==1 )
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
  
  list(L=L,k=k,v0=v0,lam=lam,eta=eta,fold=fold,BCV=BCV,opt_v0=opt_v0,opt_lam=opt_lam,opt_fit=fit,opt_biclus=Shat,time=as.numeric(time[1]))
}  


DataFABIA_BIC <- function(datapath,outpath,name,L,thrW,thrZ,thres=0.5)
{
  D1 = length(L)
  D2 = length(thrW)
  D3 = length(thrZ)

  BIC = array(0,c(D1,D2,D3))
  opt_BIC = Inf
  
  load(datapath)
  
  for ( d1 in 1:D1 )
    for ( d2 in 1:D2 )
      for ( d3 in 1:D3 )
      {
        fname = sprintf("res_%s_FABIA_%02d_%.2f_%.2f",name,L[d1],thrW[d2],thrZ[d3])
        fpath = paste(outpath,fname,sep="/")
        load(fpath)
        
        W = matrix(0,p,L[d1])
        Z = matrix(0,L[d1],n)
        
        nParam = p
        for ( l in 1:L[d1] )
        {
          Widx = bichat$numn[l,1]$numng
          Zidx = bichat$numn[l,2]$numnp
          if ( length(Widx) != 0 & length(Zidx) != 0 )
          {
            nParam = nParam + length(Widx) + length(Zidx)
            W[Widx,l] = What[Widx,l]
            Z[l,Zidx] = Zhat[l,Zidx]
          }
        }
        
        muhat = W %*% Z + mhat
        
        BIC[d1,d2,d3] = n*p*log(mean((X-muhat)^2)) + nParam*log(n*p)
  
        if ( opt_BIC > BIC[d1,d2,d3] )
        {
          opt_BIC = BIC[d1,d2,d3]
          opt_L = L[d1]
          opt_thrW = thrW[d2]
          opt_thrZ = thrZ[d3]
        }
      }
  
  list(name=name,L=L,k=k,v0=v0,lam=lam,eta=eta,BIC=BIC,opt_BIC=opt_BIC,opt_L=opt_L,opt_thrW=opt_thrW,opt_thrZ=opt_thrZ)
}  

DataFABIA_Plain <- function(datapath,outpath,name,L,thrW,thrZ)
{
  D1 = length(thrW)
  D2 = length(thrZ)
  
  load(datapath)
  time = system.time(fit <- fabia(X,L,random=0,center=2))
  
  What = fit@L
  Zhat = fit@Z
  mhat = fit@center

  for ( d1 in 1:D1 )
    for ( d2 in 1:D2 )
    {
      fname = sprintf("res_%s_FABIA_%02d_%.2f_%.2f",name,L,thrW[d1],thrZ[d2])
      fpath = paste(outpath,fname,sep="/")
      if ( !file.exists(fpath) )
      {
        bichat = extractBic(fit,thrW[d1],thrZ[d2])

        save(What,Zhat,mhat,bichat,time,file=fpath)
      }
    }    
  
}  
