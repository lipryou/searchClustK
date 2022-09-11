source("SMLSOM.R")
source("xmeans.R")
source("gmeans.R")
source("dipmeans.R")
source("pgmeans.R")
source("mixt4.R")

Q <- MixSim::MixSim(0.01, p=2, K=4, sph=T, hom=F)
A <- MixSim::simdataset(3000, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

xm <- xmeans(X, debug=1)
gm <- gmeans(X, debug=1)
dipm <- dipmeans(X, simulate=T, debug=1)
pgm <- pgmeans(X, M=1, cov_type="full", debug=1)
mem <- mixtures4(X, Mmax=20, cov_type = "full", silent=0)
sobj <- smlsom.mvn(X, 4, 5, init.type="random", silent=0)

plot(X, col=A$id+1)
points(xm$centers, pch=16)

plot(X, col=A$id+1)
points(gm$centers, pch=16)

plot(X, col=A$id+1)
points(dipm$centers, pch=16)

plot(X, col=A$id+1)
points(pgm$mus, pch=16)

plot(X, col=A$id+1)
points(mem$mus, pch=16)

plot(X, col=A$id+1)
points(sobj$mu1, pch=16)
