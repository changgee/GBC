
library(biclust)
source("Eval/gbcmetric.R")
source("Sim/SimData.R")

SimPLAID <- function(R,seed,overlap,sigma2,L,rrel,crel,batch=0)
{
  D1 = length(rrel)
  D2 = length(crel)
  
  CE = array(0,c(D1,D2,R))
  FP = array(0,c(D1,D2,R))
  FN = array(0,c(D1,D2,R))
  SEN = array(0,c(D1,D2,R))
  SPE = array(0,c(D1,D2,R))
  MCC = array(0,c(D1,D2,R))

  S = list()

  for ( r in 1:R )
  {
    data = SimData(seed+batch+r,overlap,sigma2)
    
    p = nrow(data$W)
    n = ncol(data$Z)
    
    LL = nrow(data$Z)
    S[[r]] = list()
    for ( l in 1:LL )
      S[[r]][[l]] = list(r=which(data$W[,l]!=0),c=which(data$Z[l,]!=0))
    
    for ( d1 in 1:D1 )
      for ( d2 in 1:D2 )
      {
        fit = biclust(data$X,method=BCPlaid(),cluster="b",fit.model=y~m+a+b,background=FALSE,row.release=rrel[d1],col.release=crel[d2],max.layers=L)
        
        Shat = list()
        for ( l in 1:ncol(fit@RowxNumber) )
          Shat[[l]] = list(r=which(fit@RowxNumber[,l]),c=which(fit@NumberxCol[l,]))
      
        eval = gbcmetric(Shat,S[[r]],p,n)
        
        CE[d1,d2,r] = eval$CE
        FP[d1,d2,r] = eval$FP_CE
        FN[d1,d2,r] = eval$FN_CE
        SEN[d1,d2,r] = eval$SEN_CE
        SPE[d1,d2,r] = eval$SPE_CE
        MCC[d1,d2,r] = eval$MCC_CE
      }
  }
  
  list(overlap=overlap,sigma2=sigma2,L=L,rrel=rrel,crel=crel,S=S,CE=CE,FP=FP,FN=FN,SEN=SEN,SPE=SPE,MCC=MCC)
}
