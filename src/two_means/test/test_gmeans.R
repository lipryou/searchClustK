## create share object using the following code in your bash
##  R CMD SHLIB -c test_gmeans.c gmeans.c twomeans.c list.c utils.c

dyn.load("test_gmeans.so")

X <- as.matrix(iris[,-5])
mu <- colMeans(X)
n <- nrow(X)
p <- ncol(X)
ck <- 1:n

ret <- .C("test_first_prcomp",
          X = as.double(t(X)),
          p = as.integer(p),
          ck = as.integer(ck-1),
          nk = as.integer(length(ck)),
          mu = as.double(mu),
          vec1 = double(p),
          vec2 = double(p))

X.cov <- (t(X) - mu) %*% t(t(X) - mu) / n
cov.e <- eigen(X.cov)
cov.e$values[1]
vec - ret$vec1

## normal CFD
ret <- .C("test_normalCFD", value=as.double(1.2))
pnorm(1.2)

## test ADstatic
x <- rnorm(100)
n <- length(x)

ret <- .C("test_ADstatic",
          x=as.double(x),
          n=as.integer(n),
          d=double(1))

z <- pnorm(sort(x))
A <- -n - mean((2 * (1:n) - 1) * (log(z) + log(1-rev(z))))
ret$d - A * (1 + 4/n - 25 / (n*n))

## test projection
X <- as.matrix(iris[, -5])
n <- nrow(X)
p <- ncol(X)
km <- kmeans(X, 2)
vec1 <- km$centers[1,]
vec2 <- km$centers[2,]
ck <- 1:n

ret <- .C("test_projection",
          X = as.double(t(X)),
          p = as.integer(p),
          ck = as.integer(ck-1),
          nk = as.integer(n),
          vec1 = as.double(vec1),
          vec2 = as.double(vec2),
          prj_x = double(n)
          )

vec <- vec1 - vec2
prj_x <- scale(X %*% vec / sum(vec**2))
norm(ret$prj_x - prj_x, "2")

## test gmeans
Q <- MixSim::MixSim(0.05, p=2, K=3, sph=T, hom=F)
A <- MixSim::simdataset(1000, Pi=Q$Pi, Mu=Q$Mu, S=Q$S)
X <- A$X

km <- gmeans(X, debug=1)
