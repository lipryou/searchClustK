## create share object using the following code in your bash
##  R CMD SHLIB -c test_em.c em_gmm.c
dyn.load("test_em.so")

## Test for EM

test_list_fullcov <- function(x, mus, covs) {
    M <- nrow(mus)
    p <- ncol(covs)

    ret <- .C("test_list_fullcov",
              x = as.double(x),
              mus = as.double(t(mus)),
              covs = as.double(covs),
              M = as.integer(M),
              p = as.integer(p))
}

M <- 3
Q <- MixSim::MixSim(0.05, K=M, p=2, sph=F, hom=F)
x <- matrix(c(0.3, 0.3), 1, 2)

test_list_fullcov(x, Q$Mu, Q$S)
for (m in 1:M) print(mclust::dmvnorm(x, Q$Mu[m,], Q$S[,,m], log=T))

test_list_diagcov <- function(x, mus, sigs) {
    M <- nrow(mus)
    p <- ncol(mus)

    ret <- .C("test_list_diagcov",
              x = as.double(x),
              mus = as.double(t(mus)),
              covs = as.double(t(sigs)),
              M = as.integer(M),
              p = as.integer(p))
}

M <- 3
Q <- MixSim::MixSim(0.05, K=M, p=2, sph=T, hom=F)
x <- matrix(c(0.3, 0.3), 1, 2)

s <- Q$Mu
for (m in 1:M) s[m, ] <- diag(Q$S[,,m])

test_list_diagcov(x, Q$Mu, s)
for (m in 1:M) print(mclust::dmvnorm(x, Q$Mu[m,], Q$S[,,m], log=T))

test_estep_fullcov <- function(X, props, mus, covs) {
    n <- nrow(X)
    p <- ncol(X)

    M <- nrow(mus)

    weights <- matrix(0, n, M)

    ret <- .C("test_estep_fullcov",
              data = as.double(t(X)),
              weights = as.double(weights),
              props = as.double(props),
              mus = as.double(t(mus)),
              covs = as.double(covs),
              n = as.integer(n),
              M = as.integer(M),
              p = as.integer(p))
    return(matrix(ret$weights, n, M))
}

## Estep test
n <- 3000
M <- 3
p <- 2
Q <- MixSim::MixSim(0.05, K=M, p=p, sph=F, hom=F)
A <- MixSim::simdataset(n, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

Z <- matrix(0, n, M)
mus <- X[sample(1:n, M), ]
covs <- array(cov(X), dim=c(p, p, M))
props <- rep(1/M, M)

for (m in 1:M) Z[,m] <- props[m] * mclust::dmvnorm(X, mus[m, ], covs[,,m])
Z <- Z / rowSums(Z)

ret <- test_estep_fullcov(X, props, mus, covs)

norm(Z - ret, "2")

test_estep_diagcov <- function(X, props, mus, sigs) {
    n <- nrow(X)
    p <- ncol(X)

    M <- nrow(mus)

    weights <- matrix(0, n, M)

    ret <- .C("test_estep_diagcov",
              data = as.double(t(X)),
              weights = as.double(weights),
              props = as.double(props),
              mus = as.double(t(mus)),
              covs = as.double(t(sigs)),
              n = as.integer(n),
              M = as.integer(M),
              p = as.integer(p))
    return(matrix(ret$weights, n, M))
}

## Estep test
n <- 3000
M <- 3
p <- 2
Q <- MixSim::MixSim(0.05, K=M, p=p, sph=T, hom=F)
A <- MixSim::simdataset(n, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

Z <- matrix(0, n, M)
mus <- X[sample(1:n, M), ]
covs <- Q$S
props <- rep(1/M, M)

for (m in 1:M) Z[,m] <- props[m] * mclust::dmvnorm(X, mus[m, ], covs[,,m])
Z <- Z / rowSums(Z)

s <- Q$Mu
for (m in 1:M) s[m, ] <- diag(covs[,,m])

ret <- test_estep_diagcov(X, props, mus, s)

norm(Z - ret, "2")

test_mstep_fullcov <- function(X, props, mus, covs) {
    n <- nrow(X)
    p <- ncol(X)
    M <- nrow(mus)

    weights <- matrix(0, n, M)

    ret <- .C("test_mstep_fullcov",
              data = as.double(t(X)),
              weights = as.double(weights),
              props = as.double(props),
              mus = as.double(t(mus)),
              covs = as.double(covs),
              n = as.integer(n),
              M = as.integer(M),
              p = as.integer(p))

    return(list(props=ret$props,
                mus=matrix(ret$mus, M, p, byrow=T),
                covs=array(ret$covs, dim=c(p, p, M))))

}

## Mstep test
n <- 3000
M <- 3
p <- 2
Q <- MixSim::MixSim(0.05, K=M, p=p, sph=F, hom=F)
A <- MixSim::simdataset(n, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

Z <- matrix(0, n, M)
mus <- X[sample(1:n, M), ]
covs <- array(cov(X), dim=c(p, p, M))
props <- rep(1/M, M)

for (m in 1:M) Z[,m] <- props[m] * mclust::dmvnorm(X, mus[m, ], covs[,,m])
Z <- Z / rowSums(Z)

mus_hat <- t(Z) %*% X / colSums(Z)
covs_hat <- covs
for (m in 1:M) {
    tmp <- t(X) - mus_hat[m, ]
    covs_hat[,,m] <- tmp %*% (Z[, m] * t(tmp)) / sum(Z[, m])
}
props_hat <- colSums(Z) / n

ret <- test_mstep_fullcov(X, props, mus, covs)

norm(ret$mus - mus_hat, "2")
norm(ret$covs - covs_hat, "2")
norm(ret$props - props_hat, "2")

test_mstep_diagcov <- function(X, props, mus, sigs) {
    n <- nrow(X)
    p <- ncol(X)
    M <- nrow(mus)

    weights <- matrix(0, n, M)

    ret <- .C("test_mstep_diagcov",
              data = as.double(t(X)),
              weights = as.double(weights),
              props = as.double(props),
              mus = as.double(t(mus)),
              sigs = as.double(t(sigs)),
              n = as.integer(n),
              M = as.integer(M),
              p = as.integer(p))

    return(list(props=ret$props,
                mus=matrix(ret$mus, M, p, byrow=T),
                sigs=matrix(ret$sigs, M, p, byrow=T)
                )
           )
}

## Mstep test
n <- 3000
M <- 3
p <- 2
Q <- MixSim::MixSim(0.05, K=M, p=p, sph=T, hom=F)
A <- MixSim::simdataset(n, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

Z <- matrix(0, n, M)
mus <- X[sample(1:n, M), ]
covs <- array(diag(p), dim=c(p, p, M))
props <- rep(1/M, M)

for (m in 1:M) Z[,m] <- props[m] * mclust::dmvnorm(X, mus[m, ], covs[,,m])
Z <- Z / rowSums(Z)

mus_hat <- t(Z) %*% X / colSums(Z)
covs_hat <- covs
for (m in 1:M) {
    tmp <- t(X) - mus_hat[m, ]
    covs_hat[,,m] <- tmp %*% (Z[, m] * t(tmp)) / sum(Z[, m])
}
props_hat <- colSums(Z) / n

sigs_hat <- mus_hat
for (m in 1:M) sigs_hat[m, ] <- diag(covs_hat[,,m])

s <- mus
for (m in 1:M) s[m, ] <- diag(covs[,,m])

ret <- test_mstep_diagcov(X, props, mus, s)

norm(ret$mus - mus_hat, "2")
norm(ret$sigs - sigs_hat, "2")
norm(ret$props - props_hat, "2")

em_gmm_fullcov <- function(X, M, silent=0, itmax=500, th=1e-4) {
    if (is.vector(X))
        X <- matrix(X, length(X), 1)
    n <- nrow(X)
    p <- ncol(X)

    dm <- p*(p+3)/2 ## the number of parameters

    mus <- X[sample(1:n, M),,drop=F]
    props <- rep(1, M) / M

    covs <- array(diag(p), dim=c(p, p, M))

    weights <- matrix(0, n, M)
    logliks <- rep(0, itmax)
    ret <- .C("em_gmm",
              data = as.double(t(X)),
              weights = as.double(weights),
              mus = as.double(t(mus)),
              covs = as.double(covs),
              props = as.double(props),
              cov_type = as.integer(0),
              n = as.integer(n),
              p = as.integer(p),
              M = as.integer(M),
              th = as.double(th),
              countf = integer(1),
              logliks = as.double(logliks),
              itmax = as.integer(itmax),
              silent = as.integer(silent))

    countf <- ret$countf

    logliks <- ret$logliks[1:countf]
    loglik <- tail(logliks, 1)

    weights <- matrix(ret$weights, n, M)

    mus <- matrix(ret$mus, M, p, byrow=T)
    covs <- array(ret$covs, dim=c(p, p, M))
    props <- ret$props

    list(weights=weights, props=props, mus=mus, covs=covs,
         logliks=logliks, loglik=loglik, countf=countf, cov_type = ret$cov_type)
}

## EM fullcov
n <- 3000
M <- 3
p <- 2
Q <- MixSim::MixSim(0.05, K=M, p=p, sph=T, hom=F)
A <- MixSim::simdataset(n, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

ret <- em_gmm_fullcov(X, M, silent=0, itmax=1000, th=1e-8)

props <- ret$props
mus <- ret$mus
covs <- ret$covs

Z <- matrix(0, n, M)
for (m in 1:M) Z[,m] <- props[m] * mclust::dmvnorm(X, mus[m, ], covs[,,m])

sum(log(rowSums(Z))) - ret$loglik

em_gmm_diagcov <- function(X, M, silent=0, itmax=500, th=1e-4) {
    if (is.vector(X))
        X <- matrix(X, length(X), 1)
    n <- nrow(X)
    p <- ncol(X)

    mus <- X[sample(1:n, M),,drop=F]
    props <- rep(1, M) / M

    sigs <- matrix(1, M, p)

    weights <- matrix(0, n, M)
    logliks <- rep(0, itmax)
    ret <- .C("em_gmm",
              data = as.double(t(X)),
              weights = as.double(weights),
              mus = as.double(t(mus)),
              sigs = as.double(t(sigs)),
              props = as.double(props),
              cov_type = as.integer(1),
              n = as.integer(n),
              p = as.integer(p),
              M = as.integer(M),
              th = as.double(th),
              countf = integer(1),
              logliks = as.double(logliks),
              itmax = as.integer(itmax),
              silent = as.integer(silent))

    countf <- ret$countf

    logliks <- ret$logliks[1:countf]
    loglik <- tail(logliks, 1)

    weights <- matrix(ret$weights, n, M)

    mus <- matrix(ret$mus, M, p, byrow=T)
    sigs <- matrix(ret$sigs, M, p, byrow=T)
    props <- ret$props

    list(weights=weights, props=props, mus=mus, sigs=sigs,
         logliks=logliks, loglik=loglik, countf=countf, cov_type=ret$cov_type)
}

## EM diagcov
n <- 3000
M <- 3
p <- 2
Q <- MixSim::MixSim(0.05, K=M, p=p, sph=T, hom=F)
A <- MixSim::simdataset(n, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

ret <- em_gmm_diagcov(X, M, silent=0, itmax=1000, th=1e-8)

props <- ret$props
mus <- ret$mus
sigs <- ret$sigs

covs <- array(0, dim=c(p, p, M))
for (m in 1:M) covs[,,m] <- diag(sigs[m,])

Z <- matrix(0, n, M)
for (m in 1:M) Z[,m] <- props[m] * mclust::dmvnorm(X, mus[m, ], covs[,,m])

sum(log(rowSums(Z))) - ret$loglik
