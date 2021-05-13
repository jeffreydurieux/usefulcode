FindOptimalClusPermut <- function(Pest, Ptrue){
  # find optimal cluster permutation of estimated clustering
  # compared to simulated clustering
  clus <- length(unique(Pest))

  AllPerms = permutations( n = clus , r = clus)
  nPerms = dim(AllPerms)[1]

  BestRecov = -9999
  BestPerm = -9999
  for( permtel in 1:nPerms )
  {
    if( (permtel%%50) == 0)
    {
      print( paste( "perm: " , permtel , "/" , nPerms ) )
    }

    tp = AllPerms[permtel,]
    tempRecovBlock = matrix( -9999 , 1 , 1 )

    tab <- table(Ptrue, Pest)

    tempRecovBlock[1] = sum( diag( tab[,tp] ) )

    tempRecov = mean(tempRecovBlock)

    if( permtel==1 )
    {
      BestRecov = tempRecov
      BestRecovBlock = tempRecovBlock
      BestPerm = tp
    }
    else
    {
      if( (tempRecov-BestRecov)>.0000000001 )
      {
        BestRecov = tempRecov
        BestRecovBlock = tempRecovBlock
        BestPerm = tp
      }
    }
    rm(tp,tempRecov,tempRecovBlock)
  }

  Out = list()
  Out$BestRecov = BestRecov
  Out$BestRecovBlock = BestRecovBlock
  Out$BestPerm = BestPerm
  return(Out)

}

