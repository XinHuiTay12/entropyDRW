# entropyDRW
An Entropy-based Directed Random Walk for Pathway Activity Inference Using Topological Importance and Gene Interactions

entropyDRW is an entropy-based pathway activity inference method using directed random walks on graph. It implements entropy as a parameter variable for random walking in a biological network and Entropy Weight Method (EWM) is applied in for pathway activity inference.  

Dependent packages: Depends: R (>= 3.5), igraph, Matrix, lattice, caret, e1071 

Getting started

1. Load entropyDRW package

     library(entropyDRW)

2. Read gene expression dataset and load dataset sample index and directed pathway network

    a. Read gene exp data (Gene ID) excel file as matrix
     
     There are three gene expression profiles:
     - GSE10072.txt : Lung expression profiles
     - GSE13911.txt : Stomach expression profiles
     - GSE17856.txt : Liver expression profiles
     
     Try our sample data (Lung cancer dataset) as below
     
     LungExpData_ID <- as.matrix(read.table(file="GSE10072.txt", sep="", header=T))
     LungExpData_ID

    b. Get lung cancer dataset sample index and KEGG pathways
     
     load(entropyDRW.RData)

3. Infer pathway activities for classification

    Run the function to obtain pathway expression profiles for further classification

     PathwayRWTrain <- getPathwayRWTrain(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScore)
     PathwayRWTrain
