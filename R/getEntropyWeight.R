getEntropyWeight <- function(eDRW)
{
  # Calculate entropy weight of genes based on random walk in directed pathway network
  # input:
  # eDRW: the weight of gene after random walk in directed pathway network
  # output:
  # the entropy weight of genes

  for (i in 1:length(eDRW))
  {
    weight <- 1-(eDRW)
  }
  sumt <- sum(weight)
  for (i in 1:length(weight))
  {
    EntWeight <- weight / sumt
  }
  return(EntWeight);
}
EntropyWeight <- getEntropyWeight(eDRW)
EntropyWeight
