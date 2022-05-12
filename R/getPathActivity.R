getPathwayRWTrain <- function(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScore)
{
  # Infer pathway expression profile for training set
  # input:
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # mRNA_matrix_training: the expression profile of training set
  # vertexWeight: the weight of genes from random walk
  # vertexTScore: the PCT score of each gene in the directed pathway network
  # output:
  # the pathway expression profile for training set
  # the t-test statistics of each pathways for training set

  # infer non metabolic pathway expression profile
  pathwayRW_training_gNonMetabolic <- c()
  TValueRW_training_gNonMetabolic <- matrix(NA, nrow=length(gNonMetabolic), ncol=1)
  rownames(TValueRW_training_gNonMetabolic) <- names(gNonMetabolic)
  for (i in 1 : length(gNonMetabolic))
  {
    Vpathwayi <- get.vertex.attribute(gNonMetabolic[[i]], "name", index=V(gNonMetabolic[[i]]))
    if (length(Vpathwayi) > 0)
    {
      n <- 0    # the number of differential genes in ith pathway
      pathway_training_tmp <- matrix(nrow=1, ncol=dim(mRNA_matrix_training)[2], data=0)
      TValueRW_training_tmp <- 0

      Idx_pathwayi <- c()
      for (j in 1 : length(Vpathwayi))
      {
        Idx <- which(rownames(mRNA_matrix_training)==Vpathwayi[j])
        if (length(Idx) > 0)
        {
          if (rownames(mRNA_matrix_training)[Idx] %in% names(vertexWeight))
          {
            if(vertexTScore[rownames(mRNA_matrix_training)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_training_tmp <- pathway_training_tmp + sign(vertexTScore[rownames(mRNA_matrix_training)[Idx],1]) * sum(vertexWeight[rownames(mRNA_matrix_training)[Idx]]) * mRNA_matrix_training[Idx,]
              TValueRW_training_tmp <- TValueRW_training_tmp + vertexWeight[rownames(mRNA_matrix_training)[Idx]] * abs(vertexTScore[rownames(mRNA_matrix_training)[Idx],1])
              n <- n + 1
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx)
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_training_tmp <- pathway_training_tmp / sqrt(sum(vertexWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]]^2))
        TValueRW_training_gNonMetabolic[i, 1] <- TValueRW_training_tmp / sum(vertexWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]])
        rownames(pathway_training_tmp) <- names(gNonMetabolic)[i]
        pathwayRW_training_gNonMetabolic <- rbind(pathwayRW_training_gNonMetabolic, pathway_training_tmp)
      }#else
    }#else
  }


  # infer metabolic pathway expression profile
  pathwayRW_training_gMetabolic <- c()
  TValueRW_training_gMetabolic <- matrix(NA, nrow=length(gMetabolic), ncol=1)
  rownames(TValueRW_training_gMetabolic) <- names(gMetabolic)
  for (i in 1 : length(gMetabolic))
  {
    Vpathwayi <- get.vertex.attribute(gMetabolic[[i]], "name", index=V(gMetabolic[[i]]))
    if (length(Vpathwayi) > 0)
    {
      n <- 0     # the number of differential genes in ith pathway
      pathway_training_tmp <- matrix(nrow=1, ncol=dim(mRNA_matrix_training)[2], data=0)
      TValueRW_training_tmp <- 0
      Idx_pathwayi <- c()   # record the index of nodes
      for (j in 1 : length(Vpathwayi))
      {
        Idx <- which(rownames(mRNA_matrix_training)==Vpathwayi[j])
        if (length(Idx) > 0)
        {
          if (rownames(mRNA_matrix_training)[Idx] %in% names(vertexWeight))
          {
            if(vertexTScore[rownames(mRNA_matrix_training)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_training_tmp <- pathway_training_tmp + sign(vertexTScore[rownames(mRNA_matrix_training)[Idx],1]) * sum(vertexWeight[rownames(mRNA_matrix_training)[Idx]]) * mRNA_matrix_training[Idx,]
              TValueRW_training_tmp <- TValueRW_training_tmp + vertexWeight[rownames(mRNA_matrix_training)[Idx]] * abs(vertexTScore[rownames(mRNA_matrix_training)[Idx],1])
              n <- n + 1
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx)
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_training_tmp <- pathway_training_tmp / sqrt(sum(vertexWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]]^2))
        TValueRW_training_gMetabolic[i, 1] <- TValueRW_training_tmp / sum(vertexWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]])
        rownames(pathway_training_tmp) <- names(gMetabolic)[i]
        pathwayRW_training_gMetabolic <- rbind(pathwayRW_training_gMetabolic, pathway_training_tmp)
      }#else
    }#else
  }

  # combine the nonMetabolic pathway and metabolic pathway expression profiles
  pathwayRW_training <- rbind(pathwayRW_training_gNonMetabolic, pathwayRW_training_gMetabolic)
  TValueRW_training <- rbind(TValueRW_training_gNonMetabolic, TValueRW_training_gMetabolic)
  TValueRW_training <- TValueRW_training[!is.na(TValueRW_training), ]

  Idx <- sort(TValueRW_training, decreasing = TRUE, index.return=TRUE)$ix
  TValueRW_training <- TValueRW_training[Idx]

  return(list(pathwayRW_training, TValueRW_training))
}
PathwayRWTrain <- getPathwayRWTrain(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScore)
PathwayRWTrain
