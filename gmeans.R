gmeans <- function(X, M=1, Mmax=20, alpha=0.01, debug=0) {
    if (alpha < 0 | alpha > 1)
        stop("alpha should be in [0, 1]")
    if (is.vector(X)) {
        X <- matrix(X, length(X), 1)
        message("Input is one-dimensional. This algorithm is suitable for tasks where p>1.")
    }

    n <- nrow(X)
    p <- ncol(X)

    cv <- 1.8692

    centers <- matrix(0, Mmax, p)

    km <- kmeans(X, M)
    centers[1:M,] <- km$centers
    cluster <- km$cluster

    ret <- .C("gmeans",
              X = as.double(t(X)),
              n = as.integer(n),
              p = as.integer(p),
              M = as.double(t(centers)),
              K = as.integer(M),
              cluster = as.integer(cluster),
              cv = as.double(cv),
              Kmax = as.integer(Mmax),
              debug = as.integer(debug))

    M <- ret$K
    centers <- matrix(ret$M, Mmax, p, byrow=T)
    centers <- centers[1:M,,drop=F]

    clusters <- as.integer(ret$cluster+1)

    list(M = M, centers=centers, clusters = clusters, cv = ret$cv)
}

dyn.load("lib/gmeans.so")
