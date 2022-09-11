## create share object using the following code in your bash
##  R CMD SHLIB -c test_diptest.c diptest.c dip.c utils.c

dyn.load("test_diptest.so")

## test euc_dists
set.seed(100)
X <- as.matrix(iris[sample(1:nrow(iris), 15), -5])
row.names(X) <- NULL
n <- nrow(X)
p <- ncol(X)
D <- as.matrix(dist(X))

ret <- .C("test_euc_dists",
          X = as.double(t(X)),
          n = as.integer(n),
          p = as.integer(p),
          D = double(0.5*n*(n-1))
          )

ret$D - D[upper.tri(D)]

## test diptst

x <- sort(rnorm(20))

diptest::dip(x)

ret <- .C("test_diptst",
          d = as.double(x),
          n = as.integer(length(x)),
          dip = double(1))


## test unimodarity test

seed <- 124
set.seed(seed)
Q <- MixSim::MixSim(0.01, p=2, K=2, sph=T, hom=T)
A <- MixSim::simdataset(1000, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X
n <- nrow(X)
p <- ncol(X)

D <- as.matrix(dist(X))

B <- 1000
dip_b <- rep(0, B)

set.seed(seed)
for (i in 1:B) dip_b[i] <- diptest::dip(runif(n-1))

dips <- rep(0, n)
for (i in 1:n) {
    d <- D[i, -i]
    dips[i] <- diptest::dip(d)
}

pvals <- rep(0, n)
for (i in 1:n)
    pvals[i] <- sum(dips[i] < dip_b) / B

print(sum(pvals < 1e-16) / n)

set.seed(seed)
ret <- .C("test_testCluster_unimodarity",
          X = as.double(t(X)),
          n = as.integer(n),
          p = as.integer(p),
          c = as.integer(0:(n-1)),
          a = as.double(1e-16),
          D = double(0.5*n*(n-1)))
