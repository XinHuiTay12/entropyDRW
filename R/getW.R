igraph2Adjmatrix<-function(DirectGraph)
{
  # Convert igraph to row-normalized adjacency matrix
  # input:
  # DirectGraph: the global pathway network
  # output:
  # the row-normalized adjacency matrix of global pathway network

  AdjM <- as_adjacency_matrix(DirectGraph)

  # row-normalization
  for (i in 1:dim(AdjM)[1]){
    sumr <- sum(AdjM[i,])
    if(sumr == 0)
    {
      AdjM[i,] <- numeric(length=length(AdjM[i,]))
    }
    if(sumr > 0)
    {
      AdjM[i,] <- AdjM[i,]/sumr
    }
  }
  return(AdjM)
}
W <- igraph2Adjmatrix(DirectGraph)
W

