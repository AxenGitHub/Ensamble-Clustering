library(mclust) #Expectation Maximization package
library(cluster.datasets)

rockM <- as.matrix(rock)
modelName = "EEE"
data = iris[,1:4]
modelName <- Mclust(data)$modelName
z = unmap(iris[,5])
msEst <- mstep(modelName, data, z)
names(msEst)
modelName = msEst$modelName
parameters = msEst$parameters
em(modelName, data, parameters)
parameters$variance$G <- 4



#idea, use the expextation maximization to generate the clusters
clus<-em(modelName, data, parameters)
clusNum <- clus$G
numObs <- clus$n
clusMatrix <- clus$z

for (iter in 1:numObs){
  maxVal <- max(clus$z[iter,])
  clus$z[iter,] <- (clus$z[iter,] == maxVal)+0
}

res <- mclustBIC(data, G=1:14, modelNames=modelName)
resMclust <-Mclust(data, G=4)
vector <- clus$z[,1]
length(vector)
clusMatrix
dim(clusMatrix)
