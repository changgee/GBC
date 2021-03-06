
if ( !require(isa2) )
{
  install.packages("isa2")
  library(isa2)
}
source("Eval/gbcmetric.R")
source("Sim/SimData.R")

SimISA <- function(R,seed,overlap,sigma2,batch=0)
{
  CE = rep(0,R)
  FP = rep(0,R)
  FN = rep(0,R)
  SEN = rep(0,R)
  SPE = rep(0,R)
  MCC = rep(0,R)
  
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
    
    fit = isa(data$X)
    
    Shat = list()
    for ( l in 1:ncol(fit$rows) )
      Shat[[l]] = list(r=which(fit$rows[,l]!=0),c=which(fit$columns[l,]!=0))
    
    eval = gbcmetric(Shat,S[[r]],p,n)
    
    CE[r] = eval$CE
    FP[r] = eval$FP_CE
    FN[r] = eval$FN_CE
    SEN[r] = eval$SEN_CE
    SPE[r] = eval$SPE_CE
    MCC[r] = eval$MCC_CE
  }
  
  list(overlap=overlap,sigma2=sigma2,S=S,CE=CE,FP=FP,FN=FN,SEN=SEN,SPE=SPE,MCC=MCC)
}
