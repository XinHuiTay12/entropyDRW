getEntropyProb <- function(p0)
{
  # Calculate entropy probability vector for each gene in directed pathway network
  # Entropy formula: -(pij*(log2(pij))), pij: probability of each gene (edge weight) across the graph
  # Edge weight is based on the average p0 of two adjacent nodes, (N1 + N2)/2
  # input:
  # p0: initial weight of genes in directed pathway network (H(v0))
  # output:
  # the entropy weight of genes (probability) in global pathway network

  list <- c()
  for (i in 1:length(p0))
  {
    list[[i]] <- (p0[i] + p0[i+1])/2
  }
  EdgeWeight <- unlist(list)
  last <- tail(p0, 1)
  EdgeWeight[length(EdgeWeight)] <- last # last node remained as the initial weight of genes
  prob <- EdgeWeight/sum(EdgeWeight)
  entropy <- -(prob*(log2(prob)))
  entropy[is.nan(entropy)] = 0

  return(entropy)
}
entropyProb <- getEntropyProb(p0)
entropyProb
