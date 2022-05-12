calculatePCTScore <- function(Tscore, PCScore)
{
  # Point Biserial Correlation Coefficient and T-test statistics (PCT) of each gene in the expression profile
  # input:
  # T-score: T-test statistic of each gene in the expression profile
  # PCscore: Point Biserial Correlation Coefficient (PBC) of each gene in the expression profile
  # output:
  # the PCT score of each gene in the expression profile

  PCTscore <- matrix(nrow=dim(Tscore)[1], ncol=2, data=Tscore)
  rownames(PCTscore) <- rownames(Tscore)
  for (i in 1:nrow(PCTscore))
  {
    PCTscore[i, 1] <- Tscore[i, 1] * PCScore[i, 1]
  }

  return(PCTscore)
}
PCTScore <- calculatePCTScore(Tscore, PCScore)
PCTScore
