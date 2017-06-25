source("eval/gbcmetric.R")

nr = 20
nc = 15
S1 = list()
S1[[1]] = list(r=c(1,2,3),c=c(1,2))
S1[[2]] = list(r=c(3,4,5),c=c(2,3))
S1[[3]] = list(r=c(),c=c(3,4))

S2 = list()
S2[[1]] = list(r=c(2,3),c=c(1,2))
S2[[2]] = list(r=c(3,4,5,6),c=c(2,3))

S2 = list()
S2[[1]] = list(r=c(),c=c())
S2[[2]] = list(r=c(),c=c())

gbcmetric(S2,S1,nr,nc)



spe = c(0,0.5,1,1,1)
sen = c(1,1,1,0.5,0)
AuC(sen,spe)



spe = matrix(c(0,0.5,1,0,1,1,0,1,1),3)
sen = matrix(c(1,1,1,1,1,0,1,0,0),3)
AuC(sen,spe)


