xmeans <- function(X, M=1, Mmax=20, debug=0) {
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

    ret <- .C("xmeans",
              X = as.double(t(X)),
              n = as.integer(n),
              p = as.integer(p),
              M = as.double(t(centers)),
              K = as.integer(M),
              cluster = as.integer(cluster),
              Kmax = as.integer(Mmax),
              debug = as.integer(debug))

    M <- ret$K
    centers <- matrix(ret$M, Mmax, p, byrow=T)
    centers <- centers[1:M,,drop=F]

    clusters <- as.integer(ret$cluster+1)

    bic <- xmeans.bic(X, clusters, centers)

    list(M=M, centers=centers, clusters = clusters, bic=bic)
}

xmeans.bic <- function(X, clusters, centers) {
    loglik <- function(n, p, k, sigma) {
        - 0.5 * n * log(2*pi) - 0.5 * (n * p) * log(sigma) - 0.5 * (n - k)
    }

    n <- nrow(X)
    p <- ncol(X)
    K <- nrow(centers)
    nk <- table(clusters)

    sigma <- sum((X - centers[clusters, ])**2) / (n-K)

    ll <- rep(0, K)
    for (k in 1:K)
        ll[k] <- loglik(nk[k], p, K, sigma) + nk[k] * log(nk[k] / n)

    return(sum(ll) - 0.5 * ((K - 1) + K * p + 1) * log(n))
}

dyn.load("lib/xmeans.so")
