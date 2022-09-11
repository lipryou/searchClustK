## create share object using the following code in your bash
##  R CMD SHLIB -c test_twomeans.c twomeans.c list.c

dyn.load("test_twomeans.so")

X <- as.matrix(iris[,-5])
n <- nrow(X)
p <- ncol(X)

Kmax <- 2
m1 <- X[10, ]
m2 <- X[50,]

km <- kmeans(X, centers=rbind(m1, m2))

## test create cluster array
ret <- .C("test_create_cluster_array",
          C=as.integer(km$cluster),
          n=as.integer(n),
          Kmax=as.integer(Kmax))
## test two means

ret <- .C("test_two_means",
          X = as.double(t(X)),
          p = as.integer(p),
          c = as.integer(1:n-1),
          n = as.integer(n),
          cnt1 = as.double(m1),
          cnt2 = as.double(m2))

ret$cnt1 - km$centers[1,]
ret$cnt2 - km$centers[2,]
