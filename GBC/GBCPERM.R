# function gbcperm finds the clustering on genes and samples
# based on the patterns of W and Z matrices, respecctively.
# Also finds the permutation that groups each cluster together.

gbcperm <- function(W,Z)
{
  L = ncol(W)
  p = nrow(W)
  n = ncol(Z)
  
  arg = as.data.frame(sign(abs(W)))
  clusW = apply(2^(rep(1:L,each=p)-1)*arg,1,sum)
  
  arg$decreasing = T
  permW = do.call(order,arg)

  arg = as.data.frame(sign(abs(t(Z))))
  clusZ = apply(2^(rep(1:L,each=n)-1)*arg,1,sum)

  arg$decreasing = T
  permZ = do.call(order,arg)
  
  list(permW=permW,permZ=permZ,clusW=clusW,clusZ=clusZ)
}

