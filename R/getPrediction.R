predictExprCross <- function(mRNA_matrix_training, mRNA_matrix_test,  mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
{
  # classification
  # mRNA_matrix_training: the expression profile of training set
  # mRNA_matrix_test: the expression profile of feature evaluation set
  # mRNA_matrix_crosstest: the expression profile of test set
  # Normal: the index of normal samples in training set
  # Disease: the index of disease samples in training set
  # Normal_test: the index of normal samples in validation set
  # Disease_test: the index of disease samples in validation set
  # Normal_crosstest: the index of normal samples in test set
  # Disease_crosstest: the index of disease samples in test set
  # output:
  # the prediction results

  tScoreTrain <- calculatetScore(mRNA_matrix_training, LungNormSample, LungDiseaseSample)
  tScoreTest <- calculatetScore(mRNA_matrix_test, LungNormSample_test, LungDiseaseSample_test)
  tScoreCrossTest <- calculatetScore(mRNA_matrix_crosstest, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  vertexTScoreTrain <- getVertexTScore(tScoreTrain, DirectGraph, rownames(mRNA_matrix_training))
  vertexTScoreTest <- getVertexTScore(tScoreTest, DirectGraph, rownames(mRNA_matrix_test))
  vertexTScoreCrossTest <- getVertexTScore(tScoreCrossTest, DirectGraph, rownames(mRNA_matrix_crosstest))

  # extract pathway expression profiles
  PathwayRWTrain <- getPathwayRWTrain(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScoreTrain)
  PathwayRWTest <- getPathwayRWTest(gNonMetabolic, gMetabolic, mRNA_matrix_test, vertexWeight, vertexTScoreTest)
  PathwayRWCross <- getPathwayRWCross(gNonMetabolic, gMetabolic, mRNA_matrix_crosstest, vertexWeight, vertexTScoreCrossTest)

  # calculate the t-test statistic of each pathway in the pathway expression profiles for training set
  tScoreP_pathway_train <- calculatetScore(PathwayRWTrain[[1]], LungNormSample, LungDiseaseSample)
  tScore_pathway_train <- tScoreP_pathway_train[ ,1]
  Idx <- sort(tScore_pathway_train, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_train <- tScore_pathway_train[Idx]

  # calculate the t-test statistic of each pathway in the pathway expression profiles for feature set
  tScoreP_pathway_test <- calculatetScore(PathwayRWTest[[1]], LungNormSample_test, LungDiseaseSample_test)
  tScore_pathway_test <- tScoreP_pathway_test[ ,1]
  Idx <- sort(tScore_pathway_test, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_test <- tScore_pathway_test[Idx]

  # calculate the t-test statistic of each pathway in the pathway expression profiles for test set
  tScoreP_pathway_cross <- calculatetScore(PathwayRWCross[[1]], LungNormSample_crosstest, LungDiseaseSample_crosstest)
  tScore_pathway_cross <- tScoreP_pathway_cross[ ,1]
  Idx <- sort(tScore_pathway_cross, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_cross <- tScore_pathway_cross[Idx]

  tScore_pathway <- c(tScore_pathway_train, tScore_pathway_test, tScore_pathway_cross)

  classType_training <- rep(NA, (length(LungNormSample) + length(LungDiseaseSample)))
  classType_training[LungNormSample] <- "normal"
  classType_training[LungDiseaseSample] <- "disease"
  arffRW_training <- data.frame(t(PathwayRWTrain[[1]]), "class"=classType_training, check.names=F)

  classType_test <- rep(NA, (length(LungNormSample_test) + length(LungDiseaseSample_test)))
  classType_test[LungNormSample_test] <- "normal"
  classType_test[LungDiseaseSample_test] <- "disease"
  arffRW_test <- data.frame(t(PathwayRWTest[[1]]), "class"=classType_test, check.names=F)

  classType_crosstest <- rep(NA, (length(LungNormSample_crosstest) + length(LungDiseaseSample_crosstest)))
  classType_crosstest[LungNormSample_crosstest] <- "normal"
  classType_crosstest[LungDiseaseSample_crosstest] <- "disease"
  arffRW_crosstest <- data.frame(t(PathwayRWCross[[1]]), "class" = classType_crosstest, check.names=F)

  resPredict <- list(arffRW_training, arffRW_test, arffRW_crosstest, tScore_pathway_train, tScore_pathway_test, tScore_pathway_cross)
  return(resPredict)
}
predictExpr <- predictExprCross(mRNA_matrix_training, mRNA_matrix_test,  mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
predictExpr
