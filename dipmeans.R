dipmeans <- function(X, M=1, Mmax=20, alpha=1e-16, vthd = 0.01,
                     itermax=20, debug=0, simulate=FALSE) {
    if (alpha < 0 | alpha > 1)
        stop("alpha should be in [0, 1]")
    if (vthd < 0 | vthd > 1)
        stop("vthd should be in [0, 1]")

    n <- nrow(X)
    p <- ncol(X)

    km <- kmeans(X, M)
    centers <- matrix(0, Mmax, p)
    centers[1:M,] <- km$centers
    cluster <- km$cluster

    if (simulate) {
        ret <- .C("dipmeans_sim",
                  X = as.double(t(X)),
                  n = as.integer(n),
                  p = as.integer(p),
                  M = as.double(t(centers)),
                  K = as.integer(M),
                  cluster = as.integer(cluster),
                  alpha = as.double(alpha),
                  vthd = as.double(vthd),
                  Kmax = as.integer(Mmax),
                  debug = as.integer(debug))
    } else {
        data(qDiptab, package="diptest")
        betas <- as.double(colnames(qDiptab))
        nn <- as.integer(rownames(qDiptab))

        selectp <- which.min(abs(betas + alpha - 1))

        if (length(selectp) == 0)
            stop("beta error")

        baseDips <- qDiptab[, selectp]

        ret <- .C("dipmeans_th",
                  X = as.double(t(X)),
                  n = as.integer(n),
                  p = as.integer(p),
                  M = as.double(t(centers)),
                  K = as.integer(M),
                  cluster = as.integer(cluster),
                  nn_len = as.integer(length(nn)),
                  nn = as.integer(nn),
                  qdips = as.double(baseDips),
                  vthd = as.double(vthd),
                  Kmax = as.integer(Mmax),
                  debug = as.integer(debug))

    }

    M <- ret$K
    centers <- matrix(ret$M, Mmax, p, byrow=T)
    centers <- centers[1:M,,drop=F]

    clusters <- as.integer(ret$cluster+1)

    list(M = M, centers=centers, clusters = clusters, alpha = ret$alpha, vthd = ret$vthd)
}

dyn.load("lib/dipmeans.so")
