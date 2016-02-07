library("NMF")
# generate a synthetic dataset with known classes
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts)
# perform a 3-rank NMF using the default algorithm
res <- nmf(V, 4)
basismap(res)
coefmap(res)


library("mclust")

mc <- Mclust(iris[,1:4], 3)
plot(mc, data=iris[,1:4], what=c('classification'), dimens=c(3,4))
table(iris$Species, mc$classification)
