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
