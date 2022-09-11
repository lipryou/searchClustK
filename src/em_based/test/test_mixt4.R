##
## R CMD SHLIB -c test_mixt4.c mixt4.c em_gmm.c utils.c

dyn.load("test_mixt4.so")
source("../../mixt4.R")

## test estimate mixprop
M <- 5
n <- 100
p <- 2
dmover2 <- 0.5 * p * (p+3) / 2
n_m <- 2
props <- c(0.5, 0.2, 0.25, 0.049, 0.001)
comp <- 5

ret <- .C("test_estimate_mixprop",
          comp = as.integer(comp-1),
          props = as.double(props),
          n_m = as.double(n_m),
          dmover2 = as.double(dmover2),
          n = as.integer(n),
          M = as.integer(M)
          )

props[comp] <- max(0, (n_m - dmover2) / n)
props <- props / sum(props)
ret$props - props

## test lshift_weights
n <- 10
M <- 2
Mmax <- 2
mat <- matrix(1:(n*Mmax), n, Mmax)

comp <- 1
ret <- .C("test_lshift_weights",
          weights = as.double(mat),
          comp = as.integer(comp-1),
          n = as.integer(n),
          M = as.integer(M),
          Mmax = as.integer(Mmax))

matrix(ret$weights, n, Mmax)

## X
Q <- MixSim::MixSim(0.05, p=2, K=6, sph=T, hom=F)
A <- MixSim::simdataset(3000, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

tt.max <- 50
mat <- matrix(0, tt.max, 6)
for (tt in 1:tt.max) {
    set.seed(tt + 100)
    elp1 <- system.time(ret1 <- mixtures4(X, Mmax=9, cov_type="full", silent=1, th=1e-4))
    set.seed(tt + 100)
    elp2 <- system.time(ret2 <- mixtures4(X, Mmax=9, cov_type="full", silent=1, th=1e-5))
    set.seed(tt + 100)
    elp3 <- system.time(ret3 <- mixtures4(X, Mmax=9, cov_type="full", silent=1, th=1e-6))
    mat[tt, 1] <- aricode::ARI(ret1$cluster, A$id)
    mat[tt, 2] <- elp1[3]
    mat[tt, 3] <- aricode::ARI(ret2$cluster, A$id)
    mat[tt, 4] <- elp2[3]
    mat[tt, 5] <- aricode::ARI(ret3$cluster, A$id)
    mat[tt, 6] <- elp3[3]
}
