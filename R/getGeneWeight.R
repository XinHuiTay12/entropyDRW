getHv0 <- function(GeneWeight, DirectGraph)
{
  # Calculate the initial weight of genes (H(v0)) in directed pathway network
  # input:
  # GeneWeight: the PCTscore of each gene in the directed pathway network
  # DirectGraph: the directed pathway network
  # output:
  # the initial weight of genes in directed pathway network (H(v0))

  Vertexs <- V(DirectGraph)
  W0 <- rep(0, length(Vertexs))
  names(W0) <- Vertexs$name

  for(p in 1:length(GeneWeight))
  {
    for(i in 1:length(W0))
    {
      idx <- which(names(GeneWeight[[p]]) == names(W0[i]))
      if(length(idx) > 0)
      {
        W0[i] <- GeneWeight[[p]][idx]
      }
    }
  }
  return(W0)
}
p0 <- getHv0(GeneWeight, DirectGraph)
p0
