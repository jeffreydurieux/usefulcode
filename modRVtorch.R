# Fri Dec 11 14:49:42 2020
# Author: Jeffrey Durieux, MSc

# What: torch version of modifiedRV

library(torch)
modRVtorch <- function(X, Y){

  if(nrow(X) != nrow(Y)){
    stop('Number of rows of input matrices are not equal')
  }

  X <- torch_tensor(data = X, dtype = torch_float())
  Y <- torch_tensor(data = Y, dtype = torch_float())

  #XXtilde <- ( X %*% t(X) ) - diag (diag( X %*% t(X) ) )
  #YYtilde <- ( Y %*% t(Y) ) - diag (diag( Y %*% t(Y) ) )

  XXtilde <- torch_matmul(X, torch_transpose(X, 1,2) ) -
    torch_diag (torch_diag( torch_matmul(X, torch_transpose(X, 1,2) ) ) )

  YYtilde <- torch_matmul(Y, torch_transpose(Y, 1,2) ) -
    torch_diag (torch_diag( torch_matmul(Y, torch_transpose(Y, 1,2) ) ) )



  A <- torch_dot(XXtilde$reshape(-1), YYtilde$reshape(-1))
  B <- torch_dot(XXtilde$reshape(-1), XXtilde$reshape(-1))
  C <- torch_dot(YYtilde$reshape(-1), YYtilde$reshape(-1))

  D <- torch_mul(B,C)
  E <- torch_sqrt(D)
  Ff <- torch_div(A,E)

  res <- as.array(Ff)

  return(res)
}
