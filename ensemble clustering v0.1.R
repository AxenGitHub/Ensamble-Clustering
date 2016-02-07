library(cluster.datasets)
library(skmeans)

#approach:
#user can define the number of iterations for the clustering algorithm
#and if to change the number of clusters at each iteration to simualte different clustering solutions
#in order to do so, change the chCl flag to TRUE and set the minimum number of clusters to start from
#the max number of clusters will be the one specified in numCl

#steps:
#generate consensus matrix
#introduce drop tolerance (less than 10% of the times items are clustered together, set the correlation to 0)
#markov chain? divide rows by their sum?
#run principal component analysis -> extract eigenvalues
#plot eighenvalues to determine the number of clusters  -> the ones closer to 1 before a gap should represent the "right" number of clusters
#run ensemble clustering on the consensus matrix multiple times to refine and ultimately determine the "right" amount of clusters in the data

#webpage describing clustering algorithm implemented in the cluster package
#http://cran.r-project.org/web/views/Cluster.html



#using sunspots as a simple dataset
#sunspots - Monthly Sunspot Numbers, 1749-1983

#set up variables
#number of clusters
numCl = 14
#number of runs of kmeans algorithm with different seeds
numIter = 20
#random number boundaries
minRan = 1
maxRan = 10000

#Tolerance: minimum percentage of times items needs to be clustered together to be considered,
#if the number is below tolerance, their correlation will be set to 0
tol = 0.10

#control if the algorithm runs with the same number of clusters and different seeds or if the number of clusters changes at each iteration
#TRUE - FALSE
chgCl = TRUE
minCl = 5

#use spherical kmeans, if FALSE, eucledian kmeans is used
sphere = FALSE

#standardize dataset, if TRUE, standardization is applied (not the ranged one, the variance one)
std = FALSE

#parameter: the dataset to analyze
DatasetAnalysis <- function(){
  #set the minimum number of clusters to the max number of clusters in order to execute only once the internal do-loop
  if(chgCl == FALSE){
    minCl = numCl
  }
  
  #un-comment the dataset to analyze, comment the others
  if(std == FALSE){
    #dataset <- sunspots
    dataset <- iris[,1:4]
    #dataset <- rock
  }
  else{
    #dataset <- scale(sunspots,center=FALSE,scale=TRUE)
    dataset <- scale(iris[,1:4],center=FALSE,scale=TRUE)
    #dataset <- scale(rock,center=FALSE,scale=TRUE)
  }
  

#un-comment the dataset to analyze, comment the others
#   if(std == TRUE){
#     dataset <- scale(dataset,center=FALSE,scale=TRUE)
#   }

  #run the k-means cluster numIter times and each time initialize the seed with a random number
  #if instead, the flag to change the number of clusters is active, minCl will be set to the user value and the number of iterations will change.
  for (iter in 1:numIter){
    
    #set seed
    set.seed(sample(minRan:maxRan,1))
    #run kmeans algorithm 
    
    
    for (currNumCl in minCl:numCl){
      if (sphere == TRUE){
        cl <- skmeans(dataset, currNumCl)
      }else{
        cl <- kmeans(dataset, currNumCl)
      }
      #data(european.foods)    #run this only one time to import the data
      #cl <- kmeans(european.foods[,3:18], numCl)
      
      
      #generate consensus matrix
      for (i in 1:currNumCl){
        #(cl$cluster == i)  returns TRUE for each element that is part of cluster i, FALSE otherwise
        #for the consensus matrix, we need ones and zeros, the +0 trick transforms TRUE in 1s and FALSE in 0
        clVector <- (cl$cluster == i) + 0
        #generate consensus matrix for current cluster
        consMatrixCl <- clVector %*% t(clVector)
        
        #if it's the first cluster of the first iteration, there is no general consensus matrix yet, the cluster consensus matrix is also the general one
        if(iter==1 && i==1){
          consMatrixAlg <- consMatrixCl
        }else{
          #from the second cluster onward, sum cluster consensus matrices to the general one
          consMatrixAlg <- consMatrixAlg + consMatrixCl
        }
      }
    }
  }
  
  #free unnecessary memory
  remove(clVector,consMatrixCl)
  
  #determine the minimum number of times items need to be clustered together to be considered valid clustering solutions
  #use the ceiling function to determine the first smallest integer above the tolerance limit
  if(chgCl==FALSE){
    minNumCl <- ceiling(numIter * tol)
  }else{
    minNumCl <- ceiling(numIter * minCl * tol)
  }
  #apply tolerance filtering only if it is above 1, otherwise it is an uneccessary waste of computational time
  #invisible is used to avoid printing the entire matrix after the values are replaced
  if (minNumCl > 1){
    invisible(replace(consMatrixAlg, consMatrixAlg < minNumCl, 0))
  }
  
  #divide each row by its row sum (Markov Chain appraoch)
  consMatrixAlg <- consMatrixAlg/rowSums(consMatrixAlg)
  
  #calculate eigenvalues for the consensus matrix
  #only.values -> we care only about the eigenvalues here.
  evalues <- eigen(consMatrixAlg, symmetric=TRUE,only.values=TRUE)
  
  if(chgCl==TRUE){
    title <- paste("Analysis on Dataset \n chgCl=",chgCl," minCl=",minCl,"numCl=",numCl," numIter=",numIter)
  }else{
    title <- paste("Analysis on Dataset \n chgCl=",chgCl,"numCl=",numCl," numIter=",numIter)
  }
  plot(evalues$values[1:numCl], main=title)
  return (consMatrixAlg)
}


# #analyze data with principal
# #is it necessary to center the data?
# #svddata<-svd(scale(consMatrixAlg,center=TRUE,scale=FALSE))
# 
# #should we even do a SVD or we just need the eigenvalues for this?
# svddata<-svd(consMatrixAlg)
# plot(svddata$d)
# 
# #PCA does not work for this matrix, it is sparse
# pcadata <-prcomp(consMatrixAlg)




#run algorithm on consensus matix
#run the k-means cluster numIter times and each time initialize the seed with a random number
#if instead, the flag to change the number of clusters is active, minCl will be set to the user value and the number of iterations will change.
ConsensusMatrixAnalysis <- function(){
  for (iter in 1:numIter){
    
    #set seed
    set.seed(sample(minRan:maxRan,1))
    #run kmeans algorithm 
    
    
    for (currNumCl in minCl:numCl){
      if (sphere == TRUE){
        cl <- skmeans(consMatrixAlg, currNumCl)
      }else{
        cl <- kmeans(consMatrixAlg, currNumCl)
      }
      #data(european.foods)    #run this only one time to import the data
      #cl <- kmeans(european.foods[,3:18], numCl)
      
      
      #generate consensus matrix
      for (i in 1:currNumCl){
        #(cl$cluster == i)  returns TRUE for each element that is part of cluster i, FALSE otherwise
        #for the consensus matrix, we need ones and zeros, the +0 trick transforms TRUE in 1s and FALSE in 0
        clVector <- (cl$cluster == i) + 0
        #generate consensus matrix for current cluster
        consMatrixCl <- clVector %*% t(clVector)
        
        #if it's the first cluster of the first iteration, there is no general consensus matrix yet, the cluster consensus matrix is also the general one
        if(iter==1 && i==1){
          consMatrixAlg2 <- consMatrixCl
        }else{
          #from the second cluster onward, sum cluster consensus matrices to the general one
          consMatrixAlg2 <- consMatrixAlg2 + consMatrixCl
        }
        
        #free unnecessary memory
        #remove(clVector,consMatrixCl)
      }
    }
  }
  
  consMatrixAlg <- consMatrixAlg2
  #free unnecessary memory
  remove(clVector,consMatrixCl, consMatrixAlg2)
  
  #determine the minimum number of times items need to be clustered together to be considered valid clustering solutions
  #use the ceiling function to determine the first smallest integer above the tolerance limit
  if(chgCl==FALSE){
    minNumCl <- ceiling(numIter * tol)
  }else{
    minNumCl <- ceiling(numIter * minCl * tol)
  }
  #a
  #apply tolerance filtering only if it is above 1, otherwise it is an uneccessary waste of computational time
  #invisible is used to avoid printing the entire matrix after the values are replaced
  if (minNumCl > 1){
    invisible(replace(consMatrixAlg, consMatrixAlg < minNumCl, 0))
  }
  
  #divide each row by its row sum (Markov Chain appraoch) is this how it is done?
  #is it necessary to do it? it seems that the eigenvalues are the same even if the MC approach is not used...
  consMatrixAlg <- consMatrixAlg/rowSums(consMatrixAlg)
  
  #calculate eigenvalues for the consensus matrix
  #only.values -> we care only about the eigenvalues here.
  evalues <- eigen(consMatrixAlg, symmetric=TRUE,only.values=TRUE)
  if(chgCl==TRUE){
    title <- paste("Analysis on Consensus Matrix \n chgCl=",chgCl," minCl=",minCl,"numCl=",numCl," numIter=",numIter)
  }else{
    title <- paste("Analysis on Consensus Matrix \n chgCl=",chgCl,"numCl=",numCl," numIter=",numIter)
  }
  plot(evalues$values[1:numCl],main=title)
  
  return (consMatrixAlg)
}

#The three dataset i've tested this implementation with are:
# sunspots
# iris[,1:4]
# rock

#execute functions
#execute first analysis on dataset
system.time(consMatrixAlg <-DatasetAnalysis())

#execute analysis on Consensus Matrix
system.time(consMatrixAlg<-ConsensusMatrixAnalysis())
