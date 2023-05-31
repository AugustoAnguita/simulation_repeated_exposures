library(MASS)
library(Matrix)
library(matrixcalc)

## Correlation matrix taken from Ibon

load(file="total_correlation_cen.RData")

names(totcorr_cen)

# keep only postnatal period (more exposures)

corr = totcorr_cen[totcorr_cen$Period1=="Postnatal" & totcorr_cen$Period2=="Postnatal",]

corr = corr[order(corr$Family1,corr$Exposome1,corr$Family2,corr$Exposome2),]

# put in correlation matrix format

coraux = matrix(corr$Corr,nrow=122,ncol=121,byrow=T)

cormat = matrix(NA,ncol=122,nrow=122)
cormat[1,] = c(1,coraux[1,])
for (i in 2:121) {
  cormat[i,] = c( coraux[i,1:(i-1)], 1, coraux[i,i:121])
}
cormat[122,] = c(coraux[122,],1)

names.exp = corr$Exposome1[!duplicated(corr$Exposome1)]

colnames(cormat) = names.exp
rownames(cormat) = names.exp

names.family = corr$Family1[!duplicated(corr$Exposome1)]

table(names.family)
length(table(names.family))
# 18 families

# Average within-family correlation

avg.corr.family = vector()
for (i in 1:length(unique(names.family))) {
  avg.corr.family[i] = mean(corr$Abs[corr$Family1==unique(names.family)[i] & corr$Family2==unique(names.family)[i]],na.rm=T)
}

within.corr = data.frame(Family1=unique(names.family),avg.corr=avg.corr.family)
within.corr

# sort by absolute within-family mean correlation
# It may be useful if we want to simulate exposures with high or low within-subject corr

aux.df = merge(data.frame(Family1=unique(names.family),avg.corr=avg.corr.family),
               data.frame(order1=1:122,names.exp=names.exp,Family1=names.family), by="Family1")

cormat.ordered = cormat[order(aux.df$avg.corr,aux.df$order1),order(aux.df$avg.corr,aux.df$order1)]

# e.g. low correlation:
cormat.ordered[1:5,1:5]

# e.g. high correlation:
cormat.ordered[117:121,117:121]

Sigma = cormat.ordered


####################################
## Functions to simulate dataset  ##
####################################

check.ICC1 = function(nt,index.exposures,ICC,rho) {
  # Function that checks if the combination of Sigma and ICC will produce
  #  a valid covariance matrix
  # This function checks that correlation between errors does not have to be greater than 1
  
  Sigma = Sigma[index.exposures,index.exposures]
  
  l = length(index.exposures)
  
  rho.min = matrix(NA,nrow=l,ncol=l)
  rows = vector()
  cols = vector()
  ss <- icc1 <- icc2 <- vector()
  k = 1
  for (i in 1:(l-1)){
    for (j in (i+1):l) {
      rho.min[i,j] = max( abs(Sigma[i,j])-ICC[i] , abs(Sigma[i,j])-ICC[j]) / (sqrt(1-ICC[i])*sqrt(1-ICC[j]))
      if (rho.min[i,j]>rho) {
        rows[k] = i
        cols[k] = j
        ss[k] = Sigma[i,j]
        icc1[k] = ICC[i]
        icc2[k] = ICC[j]
        k = k +1
      }
    }
  }
  problem.in = (length(rows)>0)
  
  if (problem.in) {
    res = data.frame(rows,cols,ss,icc1,icc2)
    # increase lowest ICC in all problematic pairs
    newICC = ICC
    for (i in 1:nrow(res)) {
      if (res$icc1[i]<res$icc2[i]) {
       if (res$icc1[i]==.1) newICC[rows] = 0.5 else if (res$icc1[i]==.5) newICC[rows] = 0.9  else if (res$icc1[i]==.9) newICC[rows] = 0.99 
      } else {
       if (res$icc2[i]==.1) newICC[cols] = 0.5 else if (res$icc2[i]==.5) newICC[cols] = 0.9  else if (res$icc2[i]==.9) newICC[cols] = 0.99
      }
    }
    
    # Now check if with the corrections, it is ok
    
    rho.min = matrix(NA,nrow=l,ncol=l)
    rows = vector()
    cols = vector()
    ss <- icc1 <- icc2 <- vector()
    k = 1
    for (i in 1:(l-1)){
      for (j in (i+1):l) {
        rho.min[i,j] = max( abs(Sigma[i,j])-newICC[i] , abs(Sigma[i,j])-newICC[j]) / (sqrt(1-newICC[i])*sqrt(1-newICC[j]))
        if (rho.min[i,j]>rho) {
          rows[k] = i
          cols[k] = j
          ss[k] = Sigma[i,j]
          icc1[k] = newICC[i]
          icc2[k] = newICC[j]
          k = k +1
        }
      }
    }
    problem.out = (length(rows)>0)
    
  } else {
    newICC = ICC
    problem.in = F
    problem.out = F
  }
  
  return(list(newICC=newICC, problem.in=problem.in, problem.out=problem.out))
}


rebuild.ICC = function(nt,index.exposures,ICC,rho) {
  problem.out1 = T
  n.attempts = 1
  newICC = ICC
  while (n.attempts<10 & problem.out1==T) {
    check = check.ICC1(nt=nt,index.exposures=index.exposures,ICC=ICC,rho)
    ICC = check$newICC
    problem.out1 = check$problem.out
    n.attempts = n.attempts + 1
  }
  return(list(newICC=ICC, problem.out= problem.out1))
}


simul = function(N,nt,index.exposures,ICC) {
  # N: sample size
  # nt: number of time points
  # index.exposures: the index for the exposures in cormat.ordered
  #        e.g. if we want to generate 122 exposures, index.exposures = 1:122
  #        e.g. if we want to generate 50 exposures highly correlated, index.exposures = 1:50
  #        e.g. if we want to generate 50 exposures with low correlation, index.exposures = 73:122
  #        e.g. if we want to generate 100 exposures at random, index.exposures = sample(1:122,100,replace=F)
  # ICC: vector of same size than index.exposures that indicates the intraclass correlation of the exposures.
  #       This is ICC[var1] in my notation
  #       High ICC1 means high repatebility of the variable over time.

  # Generate from a multivariate normal distribution with a correlation 
  #    matrix calculated from Sigma and ICC

  # Without loss of generality, we assume all variables have mean 0
  # We assume the same variance = 1 for all exposures - this could be changed but I'm
  #    not sure it's worth it

  Sigma = Sigma[index.exposures,index.exposures]
  
  l = length(index.exposures)

  rho.min = matrix(NA,nrow=l,ncol=l)
  for (i in 1:(l-1)){
    for (j in (i+1):l) {
      rho.min[i,j] <- rho.min[j,i] <- (abs(Sigma[i,j]) - min(ICC[i],ICC[j])) / (sqrt(1-ICC[i])*sqrt(1-ICC[j]))
    }
  }
  diag(rho.min) = 1
  
  rho.max = matrix(NA,nrow=l,ncol=l)
  for (i in 1:(l-1)){
    for (j in (i+1):l) {
      rho.max[i,j] <- rho.max[j,i] <- (abs(Sigma[i,j]) ) / (sqrt(1-ICC[i])*sqrt(1-ICC[j]))
    }
  }
  diag(rho.max) = 1
  
  rho.e = matrix(NA,nrow=l,ncol=l)
  for (i in 1:(l-1)){
    for (j in (i+1):l) {
      rho.e[i,j] <- rho.e[j,i] <- runif(1,max(0,rho.min[i,j]),min(1,rho.max[i,j]))
    }
  }
  diag(rho.e) = 1
  rho.e = rho.e*sign(Sigma)

  Psi = Sigma - rho.e * (sqrt(1-ICC) %*% t(sqrt(1-ICC))) 
  #diag(Psi) = ICC
  
  blocks = list()
  blocks[[1]] = Sigma
  for (i in 2:nt) {
    blocks[[i]] = Psi
  }
  
  m.str <- toeplitz(1:nt)
  
  res <- lapply(1:nt,function(k) {
    res <- matrix(0,ncol=ncol(m.str),nrow=nrow(m.str))
    res[m.str == k] <- 1
    res %x% blocks[[k]]
  })
  
  S = Reduce("+",res)
  
  if (!is.positive.definite(S)) S = as.matrix(nearPD(S,corr=T)$mat)
  
  X = data.frame(mvrnorm(n=N,mu=rep(0,nt*l),Sigma=S))
  names(X) = paste0(paste0("X",1:l),rep(paste0(".",1:nt),each=l))

  return(list(X=X,S=S))
}

N = 1200
set.seed(3409)
# Here I'm taking a random sample of 100 exposures
index.exposures = sample(1:122,100,replace=F)

nvars = length(index.exposures)

# Correlation between repeated measures can be high (0.9), medium (0.5), or low (0.1)
#    Each variable will have one of these 3 correlations
# Here I choose at random for each variable one of the 3 correlations 
ICC = sample(x=c(0.9,0.9,0.9),size=nvars,prob=c(1/3,1/3,1/3),replace=T)

# number of time points
nt = 5

check = rebuild.ICC(nt=nt,index.exposures=index.exposures,ICC=ICC,rho=1) 
# The following needs to be FALSE, we have to check it every time
check$problem.out
# number of changes done in ICC
sum(ICC!=check$newICC)
# we will have fewer cases of low ICCs, because they are not compatible with some correlations
table(ICC)
table(check$newICC)

ICC = check$newICC

simdata = simul(N,nt,index.exposures,ICC) 

X = simdata[[1]]
S = simdata[[2]]

# in case they are needed, names of original variables used to simulate
names.vec = rownames(cormat.ordered)[index.exposures]


## Reshape to long format

lnames = list()
for (i in 1:nvars) {
  lnames[[i]] = matrix(names(X),ncol=l,byrow=T)[,i] 
}

Xlong = reshape(X,direction="long",varying=lnames)
Xlong = Xlong[,c(dim(Xlong)[2],1:(dim(Xlong)[2]-1))]
names(Xlong) = c("id","time",paste0("X",1:nvars)) 
Xlong = Xlong[order(Xlong$id,Xlong$time),]


# What follows has not changed in the new version, except that I removed 
#   scenarios 3 and 4 (the ones involving trajectories of exposures)

#################################
## Generate response variable  ##
#################################

#############################################
#### Here for Y measured at a single point
#############################################

# Some potential scenarios

#### 1) Y depends on all the values (all times) of a set of exposures

# e.g. it depends on 3 exposures

n.true.exposures = 3

# selected at random

true.exposures = sample(1:nvars,n.true.exposures,replace=F)

# suppose the coefficients of the variables decrease with time
true.coefs = seq(2,1,by=-(1/(r-1)))

Xtrue = as.matrix(cbind(X[,paste0(paste0("X",rep(true.exposures,each=r),"."),1:r)]))
Xtrue.names = paste0(paste0("X",rep(true.exposures,each=r),"."),1:r)
  
linpred = Xtrue %*% matrix(rep(true.coefs,n.true.exposures),ncol=1)
# R2 of the model = 10%
resid.var = var(linpred)*9

Y = linpred + rnorm(N,0,sd=sqrt(resid.var))

formul = paste("Y ~", paste(Xtrue.names,sep=" ", collapse="+"))

mod = lm(formul,data=X)
summary(mod)


#### 2) Y depends on a single period of a set of exposures (to simplify I assume the same)

# I assume time=2 is the causal one

n.true.exposures = 3

# selected at random

true.exposures = sample(1:nvars,n.true.exposures,replace=F)

true.coefs = rep(1,n.true.exposures) 

Xtrue = as.matrix(cbind(X[,paste0("X",true.exposures,".2")]))
Xtrue.names = paste0("X",true.exposures,".2")

linpred = Xtrue %*% matrix(true.coefs,ncol=1)
# R2 of the model = 5%
resid.var = var(linpred)*0.95/.05

Y = linpred + rnorm(N,0,sd=sqrt(resid.var))

formul = paste("Y ~", paste(Xtrue.names,sep=" ", collapse="+"))

mod = lm(formul,data=X)
str(summary(mod))













# length(resu.sim.dataX.i)
# [1] 100
# resu.sim.dataX.i[[1]]
# List of 2
# X      :'data.frame':	1200 obs. of  500 variables      <- X
# rho.vec: num [1:100] <- 


dataY1andX
exp3
exp5
exp6
dataY2andX
exp3
exp5
exp6


class(resu.sim.data1.exp3)
str(resu.sim.data1.exp3[[1]])
# List of 5
# $ resp     : num [1:1200, 1] -42.1 29.5 30.6 -14.6 18.9 ...          <- Y
# $ true.pred: chr [1:15] "X52.1" "X52.2" "X52.3" "X52.4" ...          <- Xtrue.names
# $ sigmaE   : num [1, 1] 900                                          <- Sigma o summary(mod)$sigma o resid.var(este)
# $ beta     : num [1:15] 2 1.75 1.5 1.25 1 2 1.75 1.5 1.25 1 ...      <- true.coefs o summary(mod)$coef[2:4,1]
# $ R2       : num 0.1                                                <- summary(mod)$r.squared  o 0.1


