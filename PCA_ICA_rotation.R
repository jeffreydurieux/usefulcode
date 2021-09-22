# Wed Sep 22 09:41:02 2021 ------------------------------



# Author: Jeffrey Durieux, MSc

library(ica)
library(multiway)
library(plotly)

addError<-function(datablock,error, type = "Gaussian", additiontype = 1)
{
  if(type == "Gaussian"){
    errorM<-replicate(ncol(datablock),rnorm(nrow(datablock)))
  }else if(type == "AR"){
    nscan <- ncol(datablock)
    vdim <- nrow(datablock)^(1/3)
    dim <- rep(vdim, 3)
    errorM <- neuRosim::temporalnoise(dim = dim, nscan = nscan, sigma = 1,rho = 0.2)
    errorM <- matrix(errorM, ncol = nscan)
  }
  
  errorM<-SSequal(errorM,datablock)
  errorlevel<-error/(1-error)
  
  if(additiontype == 1){
    res<-datablock + (errorM * sqrt(errorlevel))
  }else{
    res<-datablock + (errorM * errorlevel)
  }
  return(res)
}

SSequal<-function(m1,m2)
{
  #res<-(m1/sqrt(SSsuga(m1)) * sqrt(SSsuga(m2))) #c++
  res<-(m1/sqrt(sum(m1^2)) * sqrt(sum(m2^2))) #R
  return(res)
}



# PCA using least squares 
library(ica)
nobs <- 1000
Amat <- cbind(icasamp("c","rnd",nobs),icasamp("c","rnd",nobs),icasamp("c","rnd",nobs), icasamp("c","rnd",nobs) )
Bmat <- matrix(2*runif(8),4,4)
Xmat <- tcrossprod(Amat,Bmat) 
X <- Xmat
rm(Amat,Bmat,Xmat)
X <- addError(X, error = .1)

X <- scale(X)

Q <- 3

ica <- icafast(X, nc = Q)

cor <- cor(X)
cor
R <- 1/(nobs-1)* t(X) %*% X
R
evd <- eigen(R)

B1 <- evd$vectors[,1:Q] %*% diag(1/sqrt(evd$values[1:Q])) 

B2 <- evd$vectors[,1:Q] %*% diag(1/sqrt(evd$values[1:Q])) %*% ica$R

vm <- varimax(B1)
B3 <- evd$vectors[,1:Q] %*% diag(1/sqrt(evd$values[1:Q])) %*% vm$rotmat

Ff1 <- X %*% B1
Ff2 <- X %*% B2
Ff3 <- X %*% B3


At1 <- solve(t(Ff1)%*%Ff1) %*% t(Ff1) %*% X
At2 <- solve(t(Ff2)%*%Ff2) %*% t(Ff2) %*% X
At3 <- solve(t(Ff3)%*%Ff3) %*% t(Ff3) %*% X

sum( (X - Ff1%*%At1)^2) %>% round(digits = 3)
sum( (X - Ff2%*%At2)^2) %>% round(digits = 3)
sum( (X - Ff3%*%At3)^2) %>% round(digits = 3)

sum( (X - ica$S %*% t(ica$M))^2 )

congru(Ff1, ica$S) %>% round(digits = 3)
congru(Ff2, ica$S) %>% round(digits = 3)  # ICA rotation
congru(Ff3, ica$S) %>% round(digits = 3)

t(B1) %*% R %*% B1 %>% round(digits = 3)
t(B2) %*% R %*% B2 %>% round(digits = 3)
t(B3) %*% R %*% B3 %>% round(digits = 3)

1/(nobs-1) * t(Ff1) %*% Ff1 %>% round(digits = 3)
1/(nobs-1) * t(Ff2) %*% Ff2 %>% round(digits = 3)
1/(nobs-1) * t(Ff3) %*% Ff3 %>% round(digits = 3)

# whitening check
t(ica$Y) %*% ica$Y %>% round(digits = 3)
xw <- X %*% B1
t(xw) %*% xw %>% round(digits = 3)

congru(ica$Y, xw) %>% round(digits = 3)

