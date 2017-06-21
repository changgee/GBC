

SimData <- function(seed,p,n,type,param,overlap,Edensity=0.05)
{
  size = 50
  L = 4

  meanWZ = 1.5
  sigmaWZ = 0.3

  if ( is.null(type) )
  {
    type = sample(c(0,1,3),p,TRUE)
    param = rep(0,p)
    for ( i in 1:p )
      if ( type[i] == 0 )
        param[i] = 25
      else if ( type[i] < 3 )
        param[i] = sample(1:15,1)
      else
        param[i] = 1000
  }
  else if ( length(type) == 1 )
    type = rep(type,p)
  
  if ( length(param) == 1 )
    param = rep(param,p)
  
  set.seed(seed)

  W = matrix(0,p,L)
  W[1:size,1] = rnorm(size,meanWZ,sigmaWZ)
  W[size+1:size,2] = rnorm(size,meanWZ,sigmaWZ)
  W[2*size+1:size,3] = rnorm(size,meanWZ,sigmaWZ)
  W[3*size+1:size-overlap,4] = rnorm(size,meanWZ,sigmaWZ)

  sizeZ = rpois(L,30)
  Z = matrix(0,L,n)
  for ( l in 1:L )
    Z[l,sample(1:n,sizeZ[l])] = rnorm(sizeZ[l],meanWZ,sigmaWZ)
  
  m = rep(0,p)
  for ( i in 1:p )
    if ( type[i] == 1 )
      m[i] = log(1/49)
    else if ( type[i] == 2 )
      m[i] = log(49)
  
  mu = W%*%Z + m
  X = matrix(0,p,n)
  
  for ( j in 1:p )
    if ( type[j] == 0 )
      X[j,] = mu[j,] + rnorm(n,sd=sqrt(param[j]))
    else
    {
      pr = 1/(1+exp(-mu[j,]))
      if ( type[j] == 1 )
        X[j,] = rbinom(n,param[j],pr)
      if ( type[j] == 2 )
        X[j,] = rnbinom(n,param[j],pr)
      if ( type[j] == 3 )
        X[j,] = rpois(n,exp(mu[j,]))
    }
    
  E = matrix(0,0,2)
  for ( l in 1:(p/50) )
  {
    g = 1:50+(l-1)*50
    for ( i in g )
      for ( j in g )
        if ( i != j & runif(1) < Edensity )
          E = rbind(E,c(i,j))
  }
  
  E = unique(rbind(E,E[,2:1]))
  E = E[order(E[,1],E[,2]),]
  
  list(p=p,n=n,type=type,param=param,W=W,Z=Z,m=m,X=X,E=E)
}
