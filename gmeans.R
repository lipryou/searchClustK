gmeans <- function(X, M=1, Mmax=20, cv=1.8692, debug=0) {
    ## X: input matrix
    ## M: start cluster number
    ## Mmax: maximum number of clusters
    ## cv: critical value
    ## debug: debug flag

    if (is.vector(X)) {
        X <- matrix(X, length(X), 1)
        message("Input is one-dimensional. This algorithm is suitable for tasks where p>1.")
    }

    n <- nrow(X)
    p <- ncol(X)

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
