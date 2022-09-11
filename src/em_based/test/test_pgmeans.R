## create share object using the following code in your bash
##  R CMD SHLIB -c test_pgmeans.c em_gmm.c kstest.c
dyn.load("test_pgmeans.so")

test_ksone <- function(X, M, Q) {
    pi <- Q$Pi
    mu <- as.vector(Q$Mu)
    sigma <- as.vector(Q$S)

    ret <- .C("test_ksone",
              x = as.double(X),
              n = as.integer(length(X)),
              M = as.integer(M),
              pi = as.double(pi),
              mu = as.double(mu),
              sigma = as.double(sigma)
              )
}

pmixnorm <- function(q, pi, mean, sigma) {
    G <- length(pi)

    cdf <- rep(0, length(q))
    for (g in 1:G) cdf <- cdf + pi[g] * pnorm(q, mean[g], sqrt(sigma[g]))

    return(cdf)
}

M <- 3
n <- 6000
Q <- MixSim::MixSim(0.05, K=M, p=1, sph=T, hom=F)
A <- MixSim::simdataset(n, Q$Pi, Q$Mu, Q$S)
X <- as.vector(A$X)

system.time(ksvalue <- ks.test(X, pmixnorm,
                               pi=Q$Pi, mean=as.vector(Q$Mu), sigma=as.vector(Q$S)))

system.time(test_ksone(X, M, Q))


## pgmeans

## test create_projection
p <- 2
set.seed(123)
ret <- .C("test_create_projection", p=as.integer(p))

set.seed(123)
x <- rnorm(p, 0, 1/p)
x / sqrt(sum(x*x))

## test projection
set.seed(123)
p <- 2
M <- 3
P <- rnorm(p, 0, 1/p)
P <- P / sqrt(sum(P*P))

## test prj mean
mus <- matrix(1:6, M, p)
prj_mus <- rep(0, M)

ret <- .C("test_projection_mean",
          P=as.double(P),
          mus = as.double(t(mus)),
          prj_mus = as.double(prj_mus),
          p = as.integer(p),
          M = as.integer(M))

ret$prj_mus - mus %*% P

## test prj fullcov
covs <- array(1:(M*p*p), dim=c(p, p, M))
prj_sigmas <- rep(0, M)

ret <- .C("test_projection_fullcov",
          P=as.double(P),
          covs = as.double(covs),
          prj_sigmas = as.double(prj_sigmas),
          p = as.integer(p),
          M = as.integer(M))

ret$prj_sigmas - apply(covs, 3, function(x) t(P) %*% x %*% P)

## test prj diagcov
sigs <- matrix(1:(M*p), M, p)
prj_sigmas <- rep(0, M)

ret <- .C("test_projection_diagcov",
          P=as.double(P),
          sigs = as.double(t(sigs)),
          prj_sigmas = as.double(prj_sigmas),
          p = as.integer(p),
          M = as.integer(M))

ret$prj_sigmas - sigs %*% (P * P)

## test prj x
set.seed(123)
n <- 10
X <- matrix(rnorm(n*p), n, p)
prj_x <- rep(0, n)

ret <- .C("test_projection_x",
          P=as.double(P),
          data = as.double(t(X)),
          prj_x = as.double(prj_x),
          n = as.integer(n),
          p = as.integer(p))

ret$prj_x - X %*% P

## test add_onecluster
p <- 2
M <- 3

mus <- matrix(1:(M*p), M, p)
covs <- array(1:(M*p*p), dim=c(p, p, M))
sigs <- matrix(1:(M*p), M, p)
props <- rep(1/M, M)

x <- c(7, 8)

new_props <- rep(0, (M+1))
new_mus <- rep(0, (M+1)*p)
new_covs <- rep(0, p*p*(M+1))
new_sigs <- rep(0, (M+1)*p)

## test add props

ret <- .C("test_add_onecluster_props",
          props = as.double(props),
          new_props = as.double(new_props),
          M = as.integer(M))

ret$new_props - c(props, 1/3) / sum(c(props, 1/3))

ret <- .C("test_add_onecluster_mean",
          xi = as.double(x),
          mus = as.double(t(mus)),
          new_mus = as.double(new_mus),
          p = as.integer(p),
          M = as.integer(M))

matrix(ret$new_mus, 4, 2, byrow=T) - rbind(mus, x)

ret <- .C("test_add_onecluster_fullcov",
          covs = as.double(covs),
          new_covs = as.double(new_covs),
          p = as.integer(p),
          M = as.integer(M))

array(ret$new_covs, dim=c(2, 2, 4))[,,1:3] - covs
array(ret$new_covs, dim=c(2, 2, 4))[,,4] - apply(covs, c(1, 2), mean)

ret <- .C("test_add_onecluster_diagcov",
          sigs = as.double(t(sigs)),
          new_sigs = as.double(new_sigs),
          p = as.integer(p),
          M = as.integer(M))

matrix(ret$new_sigs, M+1, p, byrow=T) - rbind(sigs, colMeans(sigs))

## test copy
p <- 2
M <- 3

new_mus <- matrix(1:(M*p), M, p)
new_covs <- array(1:(M*p*p), dim=c(p, p, M))
new_sigs <- matrix(1:(M*p), M, p)
new_props <- rep(1/M, M)

props <- rep(0, M)
mus <- matrix(0, M, p)
covs <- array(0, dim=c(p, p, M))
sigs <- matrix(0, M, p)

ret <- .C("test_copy_prop",
          new_props = as.double(new_props),
          props = as.double(props),
          M = as.integer(M))

ret$props - ret$new_props

ret <- .C("test_copy_mean",
          new_mus = as.double(new_mus),
          mus = as.double(mus),
          p = as.integer(p),
          M = as.integer(M))

ret$mus - ret$new_mus

ret <- .C("test_copy_fullcov",
          new_covs = as.double(new_covs),
          covs = as.double(covs),
          p = as.integer(p),
          M = as.integer(M))

ret$covs - ret$new_covs

ret <- .C("test_copy_diagcov",
          new_sigs = as.double(new_sigs),
          sigs = as.double(sigs),
          p = as.integer(p),
          M = as.integer(M))

ret$sigs - ret$new_sigs

## run pgmeans
n <- 3000
p <- 2
M <- 6
Q <- MixSim::MixSim(0.03, K=M, p=p, sph=T, hom=F)
A <- MixSim::simdataset(n, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

set.seed(123)
system.time(ret <- pgmeans(X, M=1, Mmax=20, cov_type="full", debug=1))
