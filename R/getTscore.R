calculatetScore <- function(normalizedMatrix, LungNormSample, LungDiseaseSample)
{
  # T-test statistic of each gene in the expression profile
  # input:
  # the normalized expression profile
  # Normal: the index of normal samples
  # Disease: the index of disease samples
  # output:
  # the t-test statistic of each gene in the expression profile

  tscore <- matrix(NA, nrow=nrow(normalizedMatrix), ncol=2)
  rownames(tscore) <- rownames(normalizedMatrix)
  for (i in 1:nrow(normalizedMatrix))
  {
    tscore_tmp <- t.test(normalizedMatrix[i, LungDiseaseSample], normalizedMatrix[i, LungNormSample], var.equal=TRUE)
    tscore[i, 1] <- tscore_tmp$statistic
    tscore[i, 2] <- tscore_tmp$p.value
  }
  return(tscore)
}
Tscore <- calculatetScore(normalizedMatrix, LungNormSample, LungDiseaseSample)
Tscore
