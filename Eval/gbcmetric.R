if ( !require(clue) )
{
  install.packages("clue")
  library(clue)
}

gbcmetric = function(S1,S2,p,n,clus_dir=2)
{
  # remove empty biclusters
  k1 = length(S1)
  k2 = length(S2)
  
  if ( k1 > 0 )
    for ( i in k1:1 )
      if ( length(S1[[i]]$r)*length(S1[[i]]$c) == 0 )
        S1[[i]] = NULL
  if ( k2 > 0 )
    for ( j in k2:1 )
      if ( length(S2[[j]]$r)*length(S2[[j]]$c) == 0 )
        S2[[j]] = NULL

  k1 = length(S1)
  k2 = length(S2)
  
  if ( clus_dir == 0 )
  {
    n = 1
    if ( k1 > 0 )
      for ( i in 1:k1 )
        S1[[i]]$c = 1
    if ( k2 > 0 )
      for ( j in 1:k2 )
        S2[[j]]$c = 1
  }
  else if ( clus_dir == 1 )
  {
    p = 1
    if ( k1 > 0 )
      for ( i in 1:k1 )
        S1[[i]]$r = 1
    if ( k2 > 0 )
      for ( j in 1:k2 )
        S2[[j]]$r = 1
  }
  
  X1 = matrix(0,p,n)
  X2 = matrix(0,p,n)
  k = max(k1,k2)
  M = matrix(0,k,k) # pairwise true positives
  M1 = matrix(0,k,k) # pairwise false positives
  M2 = matrix(0,k,k) # pairwise false negatives
  
  if ( k1 > 0 )
    for ( i in 1:k1 )
    {
      M1[i,] = length(S1[[i]]$r) * length(S1[[i]]$c)
      X1[S1[[i]]$r,S1[[i]]$c] = X1[S1[[i]]$r,S1[[i]]$c] + 1
    }

  if ( k2 > 0 )
    for ( j in 1:k2 )
    {
      M2[,j] = length(S2[[j]]$r) * length(S2[[j]]$c)
      X2[S2[[j]]$r,S2[[j]]$c] = X2[S2[[j]]$r,S2[[j]]$c] + 1
    }

  if ( k1 > 0 )
    for ( i in 1:k1 )
      if ( k2 > 0 )
        for ( j in 1:k2 )
          M[i,j] = length(intersect(S1[[i]]$r,S2[[j]]$r)) * length(intersect(S1[[i]]$c,S2[[j]]$c))
  M1 = M1 - M
  M2 = M2 - M

  # CE  
  U = sum(pmax(X1,X2))
  asn_ce = solve_LSAP(M,TRUE)
  CE = sum(M[cbind(1:k,asn_ce)])/U

  FP_CE = sum(M1[cbind(1:k,asn_ce)]) * 1.0
  FN_CE = sum(M2[cbind(1:k,asn_ce)]) * 1.0
  TP_CE = sum(M[cbind(1:k,asn_ce)]) * 1.0
  TN_CE = p*n*k - FP_CE - FN_CE - TP_CE
  SEN_CE = TP_CE/sum(X2)
  SPE_CE = TN_CE/(p*n*k-sum(X2))
  if ( (TP_CE+FP_CE == 0) | (TP_CE+FN_CE == 0) | (TN_CE+FP_CE == 0) | (TN_CE+FN_CE == 0) )
    MCC_CE = 0
  else
    MCC_CE = (TP_CE*TN_CE-FP_CE*FN_CE)/sqrt(TP_CE+FP_CE)/sqrt(TP_CE+FN_CE)/sqrt(TN_CE+FP_CE)/sqrt(TN_CE+FN_CE)
  

  # CS  
  S = M/(M1+M2+M) # similarity
  asn_cs = solve_LSAP(S,TRUE)
  CS = mean(S[cbind(1:k,asn_cs)])

  
  # elementwise  
  FP = sum(X1!=0&X2==0) * 1.0
  FN = sum(X1==0&X2!=0) * 1.0
  TP = sum(X1!=0&X2!=0) * 1.0
  TN = p*n-FP-FN-TP
  SEN = TP/(FN+TP)
  SPE = TN/(FP+TN)
  if ( (TP+FP == 0) | (TP+FN == 0) | (TN+FP == 0) | (TN+FN == 0) )
    MCC = 0
  else
    MCC = (TP*TN-FP*FN)/sqrt(TP+FP)/sqrt(TP+FN)/sqrt(TN+FP)/sqrt(TN+FN)
    
  list(CE=CE,FP_CE=FP_CE,FN_CE=FN_CE,TP_CE=TP_CE,TN_CE=TN_CE,SEN_CE=SEN_CE,SPE_CE=SPE_CE,MCC_CE=MCC_CE,CS=CS,FP=FP,FN=FN,TP=TP,TN=TN,SEN=SEN,SPE=SPE,MCC=MCC)
}


# AuC: calculate the area under the ROC curve
# sen: sensitivity must be given in nonincreasing order as the index increases
# spe: specificity must be given in nondecreasing order as the index increases

AuC <- function(sen,spe)
{
  if ( is.matrix(sen) )
  {
    r = nrow(sen)
    c = ncol(sen)
    auc = 0
    for ( i in 1:(r-1) )
      for ( j in 1:(c-1) )
        auc = auc - (sen[i,j]+sen[i,j+1]+sen[i+1,j]+sen[i+1,j+1])*(spe[i,j]-spe[i,j+1]-spe[i+1,j]+spe[i+1,j+1])
    auc = auc/4
  }
  else
  {
    l = length(sen)
    auc = 0
    for ( i in 1:(l-1) )
      auc = auc + (sen[i]+sen[i+1])*(spe[i+1]-spe[i])
    auc = auc/2
  }
  auc
}

