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

