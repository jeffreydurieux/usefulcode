# Wed Sep 30 11:01:27 2020
# Author: Jeffrey Durieux, MSc

# Example script for i/o nifti's for trying out stuff
# and do ICA on beckmann templates

library(oro.nifti)
library(ica)
library(multiway)

# example with fMRI data inside oro.nifti package
# also check demo(package = 'oro.nifti')
# e.g., demo(avg152T1)
demo(avg152T1) # plotting can take some time

nif <- readNIfTI(file.path(system.file('nifti',
                                       package = 'oro.nifti'),
                           "filtered_func_data"))

str(nif)

data <- nif@.Data # array (often 4d) check dim of nifti file

#check info of nifti
print(nif)
orthographic(nif)
image(nif)



### example low resolution spatial maps / FC patterns
### see https://www.fmrib.ox.ac.uk/datasets/royalsoc8/
### made low resolution version of the templates to save time/memory
### contact Jeffrey Durieux to receive the low res nifti file
# change path 
sall <- readNIfTI("/Users/jeffreydurieux/Documents/templates/23x_all_beckmann_8_and_csf_and_wm.nii.gz")

print(sall)

# take a slice from template 
s1 <- sall[,,12,5] # default mode network
s2 <- sall[,,12,1] # visual network
s3 <- sall[,,12,7] # right pariatal
s4 <- sall[,,12,8] # left pariatal

s1[s1<0] <- 0 # some thresholding
s2[s2<0] <- 0
s3[s3<0] <- 0
s4[s4<0] <- 0

# check image
image(s1)
image(s2)
image(s3)
image(s4)

# vectorize maps and concatenate in matrix (Voxel by Maps/component)
S <- cbind(as.vector(s1),as.vector(s2)
           ,as.vector(s3),as.vector(s4))

# mixing matrix
set.seed(42)
A <- matrix(runif(100*4, -2, 2), nrow = 100 )

# mix data
X <- tcrossprod(S, A)
dim(X)

result <- icafast(X, nc = 4)

# check congruence (note that maps are permuted, this is due the permutation ambiguity of ICA)
congru(S, result$S)
congru(A, result$M)

estS1 <- matrix(result$S[,1],nrow = 23, ncol = 28)
estS2 <- matrix(result$S[,2],nrow = 23, ncol = 28)
estS3 <- matrix(result$S[,3],nrow = 23, ncol = 28)
estS4 <- matrix(result$S[,4],nrow = 23, ncol = 28)

image(estS1)
image(estS2)
image(estS3)
image(estS4)