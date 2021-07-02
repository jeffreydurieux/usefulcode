# Fri Jul  2 09:13:11 2021
# Author: Jeffrey Durieux, MSc

# What: some code to generate resting state fMRI data

library(neuRosim)
library(e1071)
library(CICA)

##### 1. Generate time courses using neuRosim #####

Gen_TimeCourses <- function(N ,nTS , nscan, SNR = 1, noise="none", type = 'gaussian'){
  # N = number of subjects,
  # nTS = number of signals (number of components)
  # nscan = number of time points
  # SNR = signal to noise ratio
  # noise = type of noise: c("none", "white","temporal", "low-frequency", "physiological", "mixture")
  # type = 'gaussian' of 'rician'
  A <- list()
  for(nobs in 1:N){
    TS <- matrix(data = NA, nrow = nscan, ncol = nTS)
    for(nts in 1:nTS){# add a piece of code such that NaN's do not occur
      options(warn = -1) #temp warning supression
      repeat{
        TS[ ,nts] <- neuRosim::simTSrestingstate(nscan = nscan, TR = 2, SNR = SNR, noise = noise,
                                                 weights = c(.25,.25,.25,.25), type = type)
        if(any(is.finite(TS[,nts])) == TRUE){
          break
        }#end break
      }#end repeat
    }# end for nts
    options(warn = 0) #turn warnings on
    A[[nobs]] <- TS
  }# end for nobs
  return(A)
}

test1 <- Gen_TimeCourses(N = 2, nTS = 2, nscan = 150)

# time courses subject 1
plot(test1[[1]][,1], type = 'l')
plot(test1[[1]][,2], type = 'l')
# time courses subject 2
plot(test1[[2]][,1], type = 'l')
plot(test1[[2]][,2], type = 'l')


# change noise type into for example physiological (heartbeat)
test2 <- Gen_TimeCourses(N = 2, nTS = 2, nscan = 150, noise = 'low-frequency')

# time courses subject 1
plot(test2[[1]][,1], type = 'l')
plot(test2[[1]][,2], type = 'l')
# time courses subject 2
plot(test2[[2]][,1], type = 'l')
plot(test2[[2]][,2], type = 'l')



####### 2 generate spatial components/ non-gaussian signals ######
Sgen <- function(NS, NV){
  # NS = number of signals/components
  # NV = number of voxels
  return(replicate(NS,
                     ica::icasamp(dname = "b",query = "rnd", nsamp = NV) ) )
}

stest <- Sgen(5, 1000)
hist(stest[,1])
e1071::kurtosis(stest[,1])

# generate components for a number of clusters:
NClus <- 3
SR <- replicate(NClus, Sgen(NS=5, NV=1000), simplify = F)


######## 3 mix cluster specific spatial components with subject specific time courses ######
#### X_i = S_r %*% t(A_i)
# in this example 20 subjects will be generated and a 2 cluster structure is used
NClus <-  2

# this line will produce a list of NClus elements, each with 10 subjects containing a 100x5 matrix
Air <- replicate(NClus, Gen_TimeCourses(N = 10, nTS = 5,nscan = 100),simplify = F)

# this line will produce a list of NClus elements, each with a 1000 x 5 matrix
Srs <- replicate(2, Sgen(NS=5, NV=1000), simplify = F)

# mix the data (ICA model)
# Xis is a list of nclus elements each with the mixed data
Xis <- list()
for(clus in 1:length(Srs)){
  Xis[[clus]] <- lapply( seq_along(Air[[clus]]), function(anom) tcrossprod(Srs[[clus]], Air[[clus]][[anom]]))
}

# combined Xis lists into one list
X <- unlist(Xis, recursive = F)


##### add gaussian error (for example 20% noise) to X_is object to retrieve input data #####

addError<-function(datablock,error)
{
  errorM<-replicate(ncol(datablock),rnorm(nrow(datablock)))
  errorM<-SSequal(errorM,datablock)
  errorlevel<-error/(1-error)

  res<-datablock + (errorM * sqrt(errorlevel))
  return(res)
}

SSequal<-function(m1,m2)
{
  res<-(m1/sqrt(sum(m1^2)) * sqrt(sum(m2^2))) #R
  return(res)
}

Xe <- lapply(seq_along(X), function(anom) addError(X[[anom]], error = .2))

#### CICA #####
cicaout <- CICA(DataList = Xe, nStarts = 3, nComp = 5, nClus = 2, verbose = T)
summary(cicaout)

#### two step clustering #####
# step one: single subject ICAs (i.e., CICA with nClus equal to length(Xe))
icas <- CICA(DataList = Xe, nStarts = 1, nComp = 5, nClus = length(Xe), verbose = T)

# step two compute modified RV mat
cluslabs <- c(rep(1,10), rep(2,10))
simmat <- computeRVmat(DataList = icas$Sr)
mds <- cmdscale(simmat)
plot(mds, col = cluslabs)

#based on estimated time courses
simmat <- computeRVmat(DataList = icas$Ais)
mds <- cmdscale(simmat)
plot(mds, col = cluslabs)


#### GROUP ICA #####
# Is a CICA with 1 cluster
gica <- CICA(DataList = Xe, nStarts = 1, nComp = 5, nClus = 1, verbose = T)

# step two compute modified RV mat based on Ais
simmat <- computeRVmat(DataList = icas$Ais)
mds <- cmdscale(simmat)
plot(mds, col = cluslabs)


