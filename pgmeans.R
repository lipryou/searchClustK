pgmeans <- function(X, M, Mmax=20, alpha=0.001, nprj=12, cov_type = c("full", "diag"),
                    em.restart=10, itmax=1000, th=1e-4, debug=0) {
    if (alpha < 0 | alpha > 1)
        stop("alpha should be in [0, 1]")
    if (nprj < 1)
        stop("the number of projections should be greater than 1")
    if (is.vector(X)) {
        X <- matrix(X, length(X), 1)
        message("Input is one-dimensional. This algorithm is suitable for tasks where p>1.")
    }

    cov_type <- match.arg(cov_type)

    n <- nrow(X)
    p <- ncol(X)

    props <- rep(1/M, Mmax)
    mus <- X[sample(1:n, Mmax),,drop=F]

    if (cov_type == "full") {
        cov.type <- 0
        covs <- array(diag(p), dim=c(p, p, Mmax))
    } else if (cov_type == "diag") {
        cov.type <- 1
        covs <- matrix(1, Mmax, p)
    } else {
        stop("such `cov_type` does not exist.")
    }

    weights <- matrix(0, n, Mmax)

    alpha <- 0.001

    ret <- .C("pgmeans",
              data = as.double(t(X)),
              weights = as.double(weights),
              props = as.double(props),
              mus = as.double(t(mus)),
              covs = as.double(covs),
              cov_type = as.integer(cov.type),
              n = as.integer(n),
              p = as.integer(p),
              M = as.integer(M),
              Mmax = as.integer(Mmax),
              alpha = as.double(alpha),
              nprj = as.integer(nprj),
              em_restart = as.integer(em.restart),
              itmax = as.integer(itmax),
              th = as.double(th),
              debug = as.integer(debug)
              )

    M <- ret$M

    props <- ret$props[1:M]

    mus <- matrix(ret$mus[], Mmax, p, byrow=T)
    mus <- mus[1:M, ]

    if (cov_type == "full") {
        covs <- array(ret$covs, dim=c(p, p, Mmax))
        covs <- covs[,,1:M]

        Z <- matrix(0, n, M)
        for (m in 1:M) Z[,m] <- props[m] * mclust::dmvnorm(X, mus[m, ], covs[,,m])
        Z <- Z / rowSums(Z)
    } else if (cov_type == "diag") {
        covs <- matrix(ret$covs, Mmax, p, byrow=T)
        covs <- covs[1:M, ]

        Z <- matrix(0, n, M)
        for (m in 1:M) Z[,m] <- props[m] * mclust::dmvnorm(X, mus[m, ], diag(covs[m,]))
        Z <- Z / rowSums(Z)
    }

    clusters <- apply(Z, 1, which.max)

    return(list(M=M, Z=Z, props=props, mus=mus, covs=covs, clusters=clusters, cov_type=cov_type))
}


pgmeans_r <- function(X, G=2, modelNames="VVV", max.G=100, alpha=0.001, nprj=12, em.restart=10) {
    d <- ncol(X)

    z <- matrix(runif(n*G), n, G)
    z <- z / rowSums(z)
    mc <- mclust::me(X, modelNames, z)

    while (G < max.G) {
        mc.pi <- mc$parameters$pro
        mc.mean <- mc$parameters$mean
        mc.cov <- mc$parameters$variance$sigma
        G <- mc$G

        ## projection
        P <- matrix(rnorm(nprj*d, 0, 1/d), nprj, d)
        P <- P / sqrt(rowSums(P*P))

        mc.mu.prj <- P %*% mc.mean

        mc.sigma.prj <- matrix(0, nprj, G)
        for (g in 1:G) mc.sigma.prj[, g] <- rowSums(P %*% mc.cov[,,g] * P)

        X.prjs <- X %*% t(P)

        ## KS test
        test_flag <- FALSE
        for (prj_i in 1:nprj) {
            ksvalue <- ks.test(X.prjs[,prj_i], pmixnorm,
                               pi=mc.pi, mean=mc.mu.prj[prj_i,], sigma=mc.sigma.prj[prj_i,])
            test_flag <- ksvalue$p.value < alpha
            if (test_flag) {
                print(paste("G =", G, ":", "KS rejected at", ksvalue$p.value))
                break
            }
        }
        if (!test_flag)
            return(mc)

        ## add new cluster
        new_pi <- c(mc.pi, 1/G)
        new_pi <- new_pi / sum(new_pi)

        new_mean <- cbind(mc.mean, X[sample(1:n, 1),])

        new_cov <- array(0, dim=c(d, d, G+1))
        new_cov[,,1:G] <- mc.cov
        new_cov[,,G+1] <- apply(mc.cov, c(1,2), mean)

        ## execute EM algorithm with new cluster
        loglik <- -Inf
        for (tt in 1:em.restart) {
            new_mean <- cbind(mc.mean, X[sample(1:n, 1),])

            z <- matrix(0, n, G+1)
            for (g in 1:(G+1)) {
                z[, g] <- new_pi[g] * mclust::dmvnorm(X, new_mean[,g], new_cov[,,g])
            }
            z <- z / rowSums(z)

            new_mc <- mclust::me(X, "VVV", z)

            if (loglik < new_mc$loglik) {
                loglik <- new_mc$loglik
                mc <- new_mc
            }
        }
    }
    return(mc)
}

pmixnorm <- function(q, pi, mean, sigma) {
    G <- length(pi)

    cdf <- rep(0, length(q))
    for (g in 1:G) cdf <- cdf + pi[g] * pnorm(q, mean[g], sqrt(sigma[g]))

    return(cdf)
}

dyn.load("lib/pgmeans.so")
