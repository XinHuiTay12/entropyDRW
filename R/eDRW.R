rw_direct <- function(NetEntropy, entropyProb, gamma=0.7)
{
  # Calculate the weight of genes based on random walk
  # input:
  # NetEntropy: the entropy edge-weighted adjacency matrix
  # entropyProb: the initial entropy weight assigned to the the global pathway network
  # gamma: the restart probability (0.7)
  # output:
  # the weight of nodes after random walk

  # Add Ground Node, construct new adjacent matrix
  newrow <- matrix(1,1,dim(NetEntropy)[2])
  rownames(newrow) <- c("GN")
  W1 <- rbind(NetEntropy, newrow)
  newcol <- matrix(1,dim(W1)[1],1)
  colnames(newcol) <- c("GN")
  WGN <- cbind(W1,newcol) # adjacency matrix after adding ground node

  entropyProb <- t(as.matrix(entropyProb))
  entropyProb <- cbind(entropyProb,0) # The initial probability of the ground node is 0
  colnames(entropyProb)[dim(entropyProb)[2]] = "GN"

  PT <- entropyProb

  k <- 0
  delta <- 1

  # iteration
  while(is.na(delta > 1e-10))
  {
    PT1 <- (1-gamma)*WGN
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma * entropyProb)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT))
    PT <- PT4
    k <- k + 1
    delta <- FALSE
  }
  cat('converged\n')

  PT <- t(PT)
  rownames(PT) <- NULL

  # distribute the probability of the ground node back to all other nodes
  PT[1:(dim(PT)[1]-1)] <- PT[1:(dim(PT)[1]-1)] + PT[dim(PT)[1]]/(dim(PT)[1]-1)
  res <- drop(PT[1:(dim(PT)[1]-1)])

  return(res)
}
eDRW <- rw_direct(NetEntropy, entropyProb, gamma=0.7)
eDRW
