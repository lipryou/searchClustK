dyn.load("test_xmeans.so")

## test1 esitmate_sigma
X <- as.matrix(iris[,-5])
n <- nrow(X)
p <- ncol(X)
mu <- colMeans(X)

ret <- .C("test_estimate_sigma",
          X = as.double(t(X)),
          p = as.integer(p),
          c = as.integer(0:(n-1)),
          c_len = as.integer(n),
          mu = as.double(mu),
          df = as.integer(n-1),
          sigma = double(1))

sum(colSums((t(X) - mu)**2)) / (n-1) - ret$sigma

## test2 esimate_sigma
km <- kmeans(X, 2)

cluster1 <- which(km$cluster == 1) - 1
cluster2 <- which(km$cluster == 2) - 1
centers <- km$center

nk <- table(km$cluster)

ret1 <- .C("test_estimate_sigma",
          X = as.double(t(X)),
          p = as.integer(p),
          c = as.integer(cluster1),
          c_len = as.integer(nk[1]),
          mu = as.double(centers[1,]),
          df = as.integer(n-2),
          sigma = double(1))

ret2 <- .C("test_estimate_sigma",
           X = as.double(t(X)),
           p = as.integer(p),
           c = as.integer(cluster2),
           c_len = as.integer(nk[2]),
           mu = as.double(centers[2,]),
           df = as.integer(n-2),
           sigma = double(1))

ret1$sigma + ret2$sigma - sum((X - centers[km$cluster, ])**2) / (n-2)

## test loglik
X <- as.matrix(iris[,-5])
n <- nrow(X)
p <- ncol(X)
mu <- colMeans(X)

sigma <- sum(colSums((t(X) - mu)**2)) / (n-1)

ret <- .C("test_loglik",
          n = as.integer(n),
          p = as.integer(p),
          k = as.integer(1),
          sigma = as.double(sigma),
          l = double(1))

l <- - 0.5 * n * log(2*pi) - 0.5 * (n * p) * log(sigma) - 0.5 * (n - 1)
ret$l - l

## test testCluster_bic
loglik <- function(n, p, k, sigma) {
    - 0.5 * n * log(2*pi) - 0.5 * (n * p) * log(sigma) - 0.5 * (n - k)
}

X <- as.matrix(iris[,-5])
n <- nrow(X)
p <- ncol(X)
mu <- colMeans(X)
mu1 <- X[sample(1:n, 1),]
mu2 <- X[sample(1:n, 1),]

ret <- .C("test_testCluster_bic",
          X = as.double(t(X)),
          p = as.integer(p),
          c = as.integer(0:(n-1)),
          n = as.integer(n),
          mu = as.double(mu),
          vec1 = as.double(mu1),
          vec2 = as.double(mu2),
          bic_p = double(1),
          bic_c = double(1)
          )

km <- kmeans(X, rbind(mu1, mu2))
nk <- table(km$cluster)

sigma_p <- sum(colSums((t(X) - mu)**2)) / (n-1)
sigma_c <- sum((X - km$centers[km$cluster, ])**2) / (n-2)

l_p <- loglik(n, p, 1, sigma_p)
bic_p <- l_p - 0.5 * (p + 1) * log(n)

l_c <- loglik(nk[1], p, 2, sigma_c) + nk[1] * log(nk[1]) - nk[1] * log(n)
l_c <- l_c + loglik(nk[2], p, 2, sigma_c) + nk[2] * log(nk[2]) - nk[2] * log(n)
bic_c <- l_c - 0.5 * (2*p + 2) * log(n)

ret$bic_p - bic_p
ret$bic_c - bic_c

## test xmeans
Q <- MixSim::MixSim(0.01, p=2, K=3, sph=T, hom=T)
A <- MixSim::simdataset(3000, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

xm <- xmeans(X, debug=1)

tt.max <- 100
mat <- matrix(0, tt.max, 2)

for (tt in 1:tt.max) {
    set.seed(tt)
    xm <- xmeans(X, debug=0)
    mat[tt, 1] <- xm$M
    mat[tt, 2] <- xm$bic
}
