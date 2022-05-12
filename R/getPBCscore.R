calculatePBCScore <- function(normalizedMatrix, LungDiseaseLabel)
{
  # Point Biserial Correlation Coefficient (PBC) of each gene in the expression profile
  # input:
  # the normalized expression profile
  # Disease: the index of disease samples
  # output:
  # the PBC score of each gene in the expression profile

  PCscore <- matrix(NA, nrow=nrow(normalizedMatrix), ncol=2)
  rownames(PCscore) <- rownames(normalizedMatrix)
  for (i in 1:nrow(normalizedMatrix))
  {
    PCscore_tmp <- cor.test(normalizedMatrix[i,], LungDiseaseLabel)
    PCscore[i, 1] <- PCscore_tmp$estimate
    PCscore[i, 2] <- PCscore_tmp$p.value
  }
  return(abs(PCscore))
}
PCScore <- calculatePBCScore(normalizedMatrix, LungDiseaseLabel)
PCScore
