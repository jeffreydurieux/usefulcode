ComputePhiMatrix <- function ( Ytarget , Yhat ) # computes Tucker's phi among columns of input matrix
{
  Out = solve( diag( diag( t(Ytarget) %*% Ytarget ) ) ^ .5 ) %*% t(Ytarget) %*% Yhat %*% solve( diag( diag( t(Yhat) %*% Yhat ) ) ^ .5 )
  return( Out )
}
