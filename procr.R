procr <- function ( X , Ytarget )
{
  # find orthonormal rotation of X that fits Ytarget best in least squares sense

  # INPUT
  #   X (n x p1): matrix that will be rotated (orthogonally)
  #   Ytarget (n x p2): target matrix

  # if p1 >= p2: Tmat is orthonormal
  # if p1 < p2: t(Tmat) is orthonormal (Tmat is not orthonormal)

  library(MASS)

  if ( dim(X)[1] == dim(Ytarget)[1]  )
  {
    #source("ComputePhiMatrix.R")
    temp = svd( t(X) %*% Ytarget )
    #LowestDim = min( dim( t(X) %*% Ytarget ) )
    #p = temp$u[ , 1:LowestDim ]
    #q = temp$v[ , 1:LowestDim ]
    #Tmat = p %*% t(q) # rotation matrix (p1 x p2)
    Tmat = temp$u %*% t(temp$v)
    Yhat = X %*% Tmat # (n x p2)
    RotationMatrix = Tmat
    Fit = sum( ( Yhat - Ytarget ) ^ 2 )
    PhiMatrix = ComputePhiMatrix( Ytarget , Yhat )

    RotationMatrixNonorthogonal = ginv( t(X) %*% X ) %*% t(X) %*% Ytarget
    ###RotationMatrixNonorthogonal = solve( t(X) %*% X ) %*% t(X) %*% Ytarget
    YhatNonorthogonal = X %*% RotationMatrixNonorthogonal
    FitNonorthogonal = sum( ( YhatNonorthogonal - Ytarget ) ^ 2 )
    PhiMatrixNonorthogonal = ComputePhiMatrix( Ytarget , YhatNonorthogonal )

    Out = list()
    Out$Yhat = Yhat #orthogonally rotated X matrix
    Out$RotationMatrix = RotationMatrix
    Out$Fit = Fit
    Out$PhiMatrix = PhiMatrix
    Out$TuckerCongruence = mean( abs( diag( PhiMatrix ) ) )

    Out$NonOrthogonal$Yhat = YhatNonorthogonal
    Out$NonOrthogonal$RotationMatrix = RotationMatrixNonorthogonal
    Out$NonOrthogonal$Fit = FitNonorthogonal
    Out$NonOrthogonal$PhiMatrix = PhiMatrixNonorthogonal
    Out$NonOrthogonal$TuckerCongruence = mean( abs( diag( PhiMatrixNonorthogonal ) ) )
  }
  else
  {
    cat(" ",fill=TRUE)
    cat("X and Ytarget should have the same number of rows", fill=TRUE )
    cat(" ",fill=TRUE)
    Out = list()
  }
  return(Out)
}
