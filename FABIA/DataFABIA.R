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

