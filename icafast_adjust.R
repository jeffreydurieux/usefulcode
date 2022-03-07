# Tue May 18 16:42:13 2021
# Author: Jeffrey Durieux, MSc

# Note that the function below is the icafast function from the ica package written by Nathaniel Helwig
# The original function gave an error when you request nc=1
# I made a small adjustment that deals with this error (see lines: 50:54)



icafast_adjust <- function (X, nc, center = TRUE, maxit = 100, tol = 1e-06, Rmat = diag(nc), 
                            alg = c("par", "def"), fun = c("logcosh", "exp", "kur"), 
                            alpha = 1) 
{
  X <- as.matrix(X)
  nobs <- nrow(X)
  nvar <- ncol(X)
  nc <- as.integer(nc[1])
  if (nc < 1) {
    stop("Must set nc>=1 component.")
  }
  maxit <- as.integer(maxit[1])
  if (maxit < 1) {
    stop("Must set maxit>=1 iteration.")
  }
  tol <- tol[1]
  if (tol <= 0) {
    stop("Must set ctol>0.")
  }
  if (nc > min(nobs, nvar)) {
    stop("Too many components. Set nc<=min(dim(X)).")
  }
  alpha <- alpha[1]
  if (alpha < 1 | alpha > 2) {
    stop("Must set 'alpha' between 1 and 2.")
  }
  if (nrow(Rmat) != nc | ncol(Rmat) != nc) {
    stop("Input 'Rmat' must be nc-by-nc rotation matrix.")
  }
  if (center) 
    X <- scale(X, scale = FALSE)
  xeig <- eigen(crossprod(X)/nobs, symmetric = TRUE)
  nze <- sum(xeig$val > xeig$val[1] * .Machine$double.eps)
  if (nze < nc) {
    warning("Numerical rank of X is less than requested number of components (nc).\n  Number of components has been redefined as the numerical rank of X.")
    nc <- nze
    Rmat <- diag(nc)
  }
  Dmat <- sdiag(sqrt(xeig$val[1:nc]))
  
  if(nc == 1){
    Mprt <- Dmat %*% xeig$vectors[,1]  
  }else{
    Mprt <- tcrossprod(Dmat, xeig$vec[, 1:nc])        
  }
  diag(Dmat) <- 1/diag(Dmat)
  Pmat <- xeig$vec[, 1:nc] %*% Dmat
  Xw <- X %*% Pmat
  if (nc == 1L) {
    return(list(S = Xw, M = t(Mprt), W = t(Pmat), Y = Xw, Q = t(Pmat), 
                R = matrix(1), vafs = (sum(Mprt^2) * nobs)/sum(X^2), 
                iter = NA, alg = alg, fun = fun, alpha = alpha))
  }
  if (fun[1] == "kur") {
    fun1d <- function(x) {
      x^3
    }
    fun2d <- function(x) {
      3 * (x^2)
    }
  }
  else if (fun[1] == "exp") {
    fun1d <- function(x) {
      x * exp(-(x^2)/2)
    }
    fun2d <- function(x) {
      exp(-(x^2)/2) * (1 - x^2)
    }
  }
  else {
    fun1d <- function(x) {
      tanh(alpha * x)
    }
    fun2d <- function(x) {
      alpha * (1 - tanh(alpha * x)^2)
    }
  }
  if (alg[1] == "def") {
    myiters <- rep(NA, nc)
    for (j in 1:nc) {
      if (j < 2) {
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        iter <- 0
        vtol <- 1
        while (vtol > tol && iter < maxit) {
          svec <- Xw %*% Rmat[, j]
          rnew <- colMeans(Xw * matrix(fun1d(svec), nobs, 
                                       nc))
          rnew <- rnew - mean(fun2d(svec)) * Rmat[, j]
          rnew <- rnew/sqrt(sum(rnew^2))
          vtol <- 1 - abs(sum(Rmat[, j] * rnew))
          iter <- iter + 1
          Rmat[, j] <- rnew
        }
        myiters[j] <- iter
      }
      else {
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        svec <- matrix(0, nc, 1)
        for (k in 1:(j - 1)) {
          svec <- svec + sum(Rmat[, k] * Rmat[, j]) * 
            Rmat[, k]
        }
        Rmat[, j] <- Rmat[, j] - svec
        Rmat[, j] <- Rmat[, j]/sqrt(sum(Rmat[, j]^2))
        iter <- 0
        vtol <- 1
        while (vtol > tol && iter < maxit) {
          svec <- Xw %*% Rmat[, j]
          rnew <- colMeans(Xw * matrix(fun1d(svec), nobs, 
                                       nc))
          rnew <- rnew - mean(fun2d(svec)) * Rmat[, j]
          rnew <- rnew/sqrt(sum(rnew^2))
          svec <- matrix(0, nc, 1)
          for (k in 1:(j - 1)) {
            svec <- svec + sum(Rmat[, k] * rnew) * Rmat[, 
                                                        k]
          }
          rnew <- rnew - svec
          rnew <- rnew/sqrt(sum(rnew^2))
          vtol <- 1 - abs(sum(Rmat[, j] * rnew))
          iter <- iter + 1
          Rmat[, j] <- rnew
        }
        myiters[j] <- iter
      }
    }
  }
  else {
    rsvd <- svd(Rmat)
    Rmat <- tcrossprod(rsvd$u, rsvd$v)
    iter <- 0
    vtol <- 1
    while (vtol > tol && iter < maxit) {
      smat <- Xw %*% Rmat
      rnew <- crossprod(Xw, fun1d(smat))/nobs
      rnew <- rnew - Rmat %*% sdiag(colMeans(fun2d(smat)))
      rsvd <- svd(rnew)
      rnew <- tcrossprod(rsvd$u, rsvd$v)
      vtol <- 1 - min(abs(colSums(Rmat * rnew)))
      iter <- iter + 1
      Rmat <- rnew
    }
    myiters <- iter
  }
  M <- crossprod(Rmat, Mprt)
  vafs <- rowSums(M^2)
  ix <- sort(vafs, decreasing = TRUE, index.return = TRUE)$ix
  M <- M[ix, ]
  Rmat <- Rmat[, ix]
  vafs <- (vafs[ix] * nobs)/sum(X^2)
  return(list(S = Xw %*% Rmat, M = t(M), W = t(Pmat %*% Rmat), 
              Y = Xw, Q = t(Pmat), R = Rmat, vafs = vafs, iter = myiters, 
              alg = alg[1], fun = fun, alpha = alpha))
}

# set.seed(123)
# nobs <- 1000
# Amat <- cbind(icasamp("a","rnd",nobs),icasamp("b","rnd",nobs))
# Bmat <- matrix(2*runif(4),2,2)
# Xmat <- tcrossprod(Amat,Bmat)
# 
# # ICA via FastICA with 2 components
# test1 <- icafast(Xmat, 1) # error
# test2 <- icafast_adjust(Xmat, 1)
# 
# multiway::congru(Amat, test2$S)
