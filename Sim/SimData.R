

SimData <- function(seed,p,n,type,param,overlap,L,Edensity=0.05)
{
  size = 50

  meanWZ = 1.5
  sigmaWZ = 0.1

  set.seed(seed)

  if ( is.null(type) )
  {
    type = sample(0:3,p,TRUE)
    param = rep(0,p)
    for ( i in 1:p )
      if ( type[i] == 0 )
        param[i] = 4
      else if ( type[i] == 1 )
        param[i] = sample(1:20,1)
      else if ( type[i] == 2 )
        param[i] = sample(1:20,1)
      else
        param[i] = 20000
  }
  else if ( length(type) == 1 )
    type = rep(type,p)
  
  if ( length(param) == 1 )
    param = rep(param,p)
  

  W = matrix(0,p,L)
  if ( L > 1 )
    for ( l in 1:(L-1) )
      W[(l-1)*size+1:(size+overlap),l] = rnorm(size+overlap,meanWZ,sigmaWZ) * sample(c(-1,1),size+overlap,TRUE)
  W[(L-1)*size+1:size,L] = rnorm(size,meanWZ,sigmaWZ) * sample(c(-1,1),size,TRUE)

  sizeZ = rpois(L,30)
  Z = matrix(0,L,n)
  for ( l in 1:L )
    Z[l,sample(1:n,sizeZ[l])] = rnorm(sizeZ[l],meanWZ,sigmaWZ) * sample(c(-1,1),sizeZ[l],TRUE)
  
  m = rep(0,p)
  m[type==3] = 3.5

  mu = W%*%Z + m
  X = matrix(0,p,n)
  
  for ( j in 1:p )
    if ( type[j] == 0 )
      X[j,] = mu[j,] + rnorm(n,sd=sqrt(param[j]))
    else if ( type[j] == 1 )
      X[j,] = rbinom(n,param[j],1/(1+exp(-mu[j,])))
    else if ( type[j] == 2 )
      X[j,] = rnbinom(n,param[j],1/(1+exp(mu[j,])))
    else if ( type[j] == 3 )
      X[j,] = rpois(n,exp(mu[j,]))

  E = matrix(0,0,2)
  for ( l in 1:(p/size) )
  {
    g = 1:size+(l-1)*size
    for ( i in g )
      for ( j in g )
        if ( i != j & runif(1) < Edensity )
          E = rbind(E,c(i,j))
  }
  
  E = unique(rbind(E,E[,2:1]))
  E = E[order(E[,1],E[,2]),]
  
  list(p=p,n=n,type=type,param=param,overlap=overlap,L=L,W=W,Z=Z,m=m,mu=mu,X=X,E=E)
}

