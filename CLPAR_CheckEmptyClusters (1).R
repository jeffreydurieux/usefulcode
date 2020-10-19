CLPAR_CheckEmptyClusters <- function( Data3D , tempRes , nObjTotal , nClust , nCompVector , OrthVector , StartSeed )
{
  # source("FindOneEmptyCluster.R")
  # source("CLPAR_ComputeModel.R")

  set.seed( StartSeed )
  AnalSeeds = GenerateSeed( 100 , StartSeed )
  SeedNr = 1

  tempClustering = vector( length = nObjTotal )
  for( tel in 1:nObjTotal )
  {
    tempClustering[tel] = tempRes$Clust[tel]
  }
  tempersol = tempRes

  while( length(unique(tempClustering)) != nClust ) # while empty clusters
  {
    #BadFit = matrix( 0 , 1 , nObjTotal) # initialize vector with fitvalues
    #for( rowtel in 1:nObjTotal ) # for each row-element
    #{
    #  tempElem = which( which( tempClustering==tempClustering[rowtel] ) == rowtel )
    #  tempModel3D = ComputePARAFACmodel( tempersol.A{tempClustering(rowtel)}(tempElem,:) , tempersol.B{tempClustering(rowtel)} , tempersol.C{tempClustering(rowtel)} , MakeCoreMatrixPARAFAC( InputInf.RankArray(tempClustering(rowtel)) ) )
    #  BadFit(1,rowtel) = ssq( DataAr(rowtel,:,:) - tempModel3D ) #compute fit value for each row-element
    #}

    tempCl = sort( tempRes$TotalSol$fObject , decreasing=T , index.return=T ) # list row-element number in descending order for their fit (first row-element is worst fit)
    BlockToSeparateCluster = tempCl$ix
    rm(tempCl)
    position = 1 #position=1 implies the worst fitting row-element (row-elementnr is in 'BlockToSeparateCluster'-vector)

    while( sum( tempClustering == tempClustering[BlockToSeparateCluster[position]] ) == 1 ) # zolang her-allocatie van een row-element tot een lege cluster leidt
    {
      position = position + 1
    }

    tempClustering[BlockToSeparateCluster[position]] = FindOneEmptyCluster( tempClustering , nClust) #search for an empty cluster and assign rater to it
    tempersol = CLPAR_ComputeModel_3D( Data3D , tempClustering , nClust ,  nCompVector , OrthVector , AnalSeeds[SeedNr] , nRuns = 50 , conv = .0000001 , maxit = 100 ) #re-compute model with the 'updated' clustering
    SeedNr = SeedNr + 1
  }

  return(tempersol)
}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

FindOneEmptyCluster <- function(Clust,nClust)
{
  PresentClusters = unique(Clust) #list of all clusters than are non-empty
  AllClusters = 1:nClust #all possible clusters (empty and non-empty)
  m1 = matrix( rep( AllClusters , each = length(PresentClusters) ) , length(PresentClusters) , length(AllClusters) )
  m2 = matrix( rep( PresentClusters , times = length(AllClusters) ) , length(PresentClusters) , length(AllClusters) )
  if(length(unique(Clust)) == 1){
    tempres = which( (m1 == m2) == 0 ) #denotes the empty clusters
  }else{
    tempres = which( sum(m1 == m2) == 0 )
  }

  res = tempres[ ceiling(runif(1)*length(tempres)) ] #pick one empty cluster at random

  # BlockToSeparateCluster = find(BadFit == max(BadFit));
  # if ( length(BlockToSeparateCluster) > 1 )    % neem ??n random block als er meerdere het slechtste zijn
  #     BlockToSeparateCluster = BlockToSeparateCluster( ceil(rand*length(BlockToSeparateCluster)) );
  return(res)
}
