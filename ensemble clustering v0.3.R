library(cluster.datasets)
library(skmeans) #spherical kmeans package
library(mclust) #Expectation Maximization package

#approach:
#user can define the number of iterations for the clustering algorithm
#and if to change the number of clusters at each iteration to simualte different clustering solutions
#in order to do so, change the chCl flag to TRUE and set the minimum number of clusters to start from
#the max number of clusters will be the one specified in numCl
#choose between spherical or euclidean clusterting algorithm
#choose if to standardize or not the dataset

#steps:
#generate consensus matrix
#introduce drop tolerance (less than 10% of the times items are clustered together, set the correlation to 0)
#markov chain? divide rows by their sum?
#run principal component analysis -> extract eigenvalues
#plot eighenvalues to determine the number of clusters  -> the ones closer to 1 before a gap should represent the "right" number of clusters
#run ensemble clustering on the consensus matrix multiple times to refine and ultimately determine the "right" amount of clusters in the data

#webpage describing clustering algorithm implemented in the cluster package
#http://cran.r-project.org/web/views/Cluster.html


########################################
# Set up global variables
#######################################

#max number of clusters
numCl_global = 14
#number of runs of kmeans algorithm with different seeds
numIter_global = 20
#random number generator boundaries
minRan = 1
maxRan = 10000

#Tolerance: minimum percentage of times items needs to be clustered together to be considered,
#if the number is below tolerance, their correlation will be set to 0
tol_global = 0.20

#control if the algorithm runs with the same number of clusters if the number of clusters changes at each iteration
#TRUE - FALSE
chgCl_global = TRUE
#min number of clusters to start from if chgCl is set to TRUE
minCl_global = 4

#use spherical kmeans, if FALSE, eucledian kmeans is used
sphere_global = TRUE

#standardize dataset, if TRUE, standardization is applied (not the ranged one, the variance one)
std_global = FALSE

#######################################
# Define functions 
#######################################

#Generate first consensus matrix from original dataset
#parameter: 
#dataset: the dataset to analyze in matrix format
#datasetName: the name of the dataset to use in the plot, by default is null
#all the other parameters are a replica of the global parameter
DatasetAnalysis <- function(dataset, datasetName=NULL, chgCl=chgCl_global, minCl=minCl_global, numCl=numCl_global, numIter=numIter_global, std=std_global, sphere=sphere_global, tol=tol_global){
  #set the minimum number of clusters to the max number of clusters in order to execute only once the internal do-loop
  if(chgCl == FALSE){
    minCl = numCl
  }
  
  if(std == TRUE){
    dataset <- scale(dataset,center=FALSE,scale=TRUE)
  }
  
  #if skemans is requested, remove the rows and columns that sum up to zero, otherwise the skemans function won't work
  if(sphere == TRUE){
    #dataset <- as.matrix(dataset)
    dataset <- dataset[rowSums(dataset) != 0, colSums(dataset) != 0]
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
        cl <- skmeans(dataset, currNumCl, method='pclust', control=list(verbose=FALSE))
      }else{
        cl <- kmeans(dataset, currNumCl)
      }
            
      
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
    title <- paste("Analysis on Dataset ",datasetName,"\n chgCl=",chgCl," minCl=",minCl,"numCl=",numCl," numIter=",numIter,"\n Sphere=",sphere," Std=",std)
  }else{
    title <- paste("Analysis on Dataset ",datasetName,"\n chgCl=",chgCl,"numCl=",numCl," numIter=",numIter,"\n Sphere=",sphere," Std=",std)
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




#Run algorithm on consensus matix
#parameter: 
#datasetName: the name of the dataset to use in the plot, by default is null
#all the other parameters are a replica of the global parameter

ConsensusMatrixAnalysis <- function(datasetName=NULL, chgCl=chgCl_global, minCl=minCl_global, numCl=numCl_global, numIter=numIter_global, sphere=sphere_global, tol=tol_global){
  #set the minimum number of clusters to the max number of clusters in order to execute only once the internal do-loop
  if(chgCl == FALSE){
    minCl = numCl
  }
  
  for (iter in 1:numIter){ 
    #set seed
    set.seed(sample(minRan:maxRan,1))
    #run kmeans algorithm 
    
    for (currNumCl in minCl:numCl){
      if (sphere == TRUE){
        cl <- skmeans(consMatrixAlg, currNumCl, method='pclust', control=list(verbose=FALSE))
      }else{
        cl <- kmeans(consMatrixAlg, currNumCl)
      }

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
      }
    }
  }
  
  consMatrixAlg <- consMatrixAlg2
  #free unecessary memory
  remove(clVector,consMatrixCl, consMatrixAlg2)
  
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
  if (is.null(datasetName) == TRUE){
    if(chgCl==TRUE){
      title <- paste("Analysis on Consensus Matrix \n chgCl=",chgCl," minCl=",minCl,"numCl=",numCl," numIter=",numIter,"\n Sphere=",sphere," Std=",std)
    }else{
      title <- paste("Analysis on Consensus Matrix \n chgCl=",chgCl,"numCl=",numCl," numIter=",numIter,"\n Sphere=",sphere," Std=",std)
    }
  }else{
    if(chgCl==TRUE){
      title <- paste("Dataset: ",datasetName,"\n Analysis on Consensus Matrix \n chgCl=",chgCl," minCl=",minCl,"numCl=",numCl," numIter=",numIter,"\n Sphere=",sphere," Std=",std)
    }else{
      title <- paste("Dataset: ",datasetName,"\n Analysis on Consensus Matrix \n chgCl=",chgCl,"numCl=",numCl," numIter=",numIter,"\n Sphere=",sphere," Std=",std)
    }    
  }
  plot(evalues$values[1:numCl],main=title)
  
  return (consMatrixAlg)
}


#######################################
# Call functions
#######################################

#The three dataset i've tested this implementation with are:
# sunspots -> matrix(as.matrix(sunspots),ncol=12,byrow=TRUE))
# iris[,1:4] -> as.matrix(iris[.1:4])
# rock -> as.matrix(rock)

#from http://archive.ics.uci.edu/ml/datasets/Turkiye+Student+Evaluation#
student_eval <- read.csv(file="turkiye-student-evaluation_R_Specific.csv", header=TRUE, sep=",")

#execute functions
#execute first analysis on dataset -> be sure to transform it into a matrix for the skmeans to work
#for some reason sunspots dataset is not transformed properly, it needs to be formatted
system.time(consMatrixAlg <-DatasetAnalysis(as.matrix(iris[.1:4]), datasetName = "Iris"))

#execute analysis on Consensus Matrix
#Indicate if to use the spherical kmeans or the euclidean kmeans (like the glogal sphere flag, TRUE is spherical Kmeans, FALSE euclidean)
system.time(consMatrixAlg<-ConsensusMatrixAnalysis(sphere=TRUE, datasetName = "Iris", chgCl=TRUE))
