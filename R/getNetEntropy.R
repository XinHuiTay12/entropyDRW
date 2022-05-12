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

EntAdjmatrix<-function(W)
{
  # Convert igraph to entropy of row-normalized adjacency matrix
  # input:
  # W: the row-normalized adjacency matrix
  # output:
  # the entropy row-normalized adjacency matrix of global pathway network
  # the entropy of edge-weighted graph

  # Convert W sparse matrix to matrix
  NetEntropy <- as.matrix(W)

  # get adjacency information entropy
  for (i in 1:dim(W)[1]){
    NetEntropy[i,] <- -(W[i,]*(log2(W[i,])))
  }
  # Replace Nan to 0
  NetEntropy[is.nan(NetEntropy)] <- 0
  # Convert matrix to sparse matrix
  NetEntropy <- as(NetEntropy, "sparseMatrix")

  return(NetEntropy)
}
NetEntropy <- EntAdjmatrix(W)
NetEntropy
