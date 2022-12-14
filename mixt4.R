mixtures4 <- function(X, Mmax, Mmin=1, cov_type = c("full", "diag"),
                      itmax=2000, th=1e-4, silent=0) {

    n <- nrow(X)
    p <- ncol(X)

    cov_type = match.arg(cov_type)

    ## initialization model parameters

    ## props: uniform mixing probability
    props <- rep(1 / Mmax, Mmax)

    ## means: randomly chosen data points
    mus <- X[sample(1:n, Mmax),]

    if (cov_type == "full") {

        dm <- p * (p + 3) / 2
        cov.type <- 0

        ## covs: initialization is taken from authors code.
        covs <- array(0, dim=c(p, p, Mmax))

        ## diagonal matrices proportional to 1/10 of the maximum global variance
        gcov <- cov(X) ## global covariance
        for (m in 1:Mmax)
            covs[,,m] <- diag(p) * max(diag(gcov/10))

    } else if (cov_type == "diag") {
        dm <- 2*p
        cov.type <- 1

        gcov <- cov(X)
        covs <- matrix(max(diag(gcov/10)), Mmax, p)
    }

    likelihoods <- matrix(0, n, Mmax)
    dl <- rep(0, itmax)
    logliks <- rep(0, itmax)
    kappas <- c(Mmax, rep(0,itmax-1))

    ret <- .C("mixtures4",
              data = as.double(t(X)),
              likelihoods = as.double(likelihoods),
              mus = as.double(t(mus)),
              covs = as.double(covs),
              props = as.double(props),
              cov_type = as.integer(cov.type),
              n = as.integer(n),
              p = as.integer(p),
              dmover2 = as.double(0.5 * dm),
              M = as.integer(Mmax),
              Mmin = as.integer(Mmin),
              th = as.double(th),
              countf = integer(1),
              dl = as.double(dl),
              logliks = as.double(logliks),
              kappas = as.integer(kappas),
              trans1 = integer(itmax),
              trans2 = integer(itmax),
              lives = as.integer(rep(1, Mmax)),
              blives = as.integer(rep(1, Mmax)),
              bmus = double(length(mus)),
              bcovs = double(length(covs)),
              bprops = double(Mmax),
              itmax = as.integer(itmax),
              silent = as.integer(silent))

    countf <- ret$countf
    dl <- ret$dl[1:countf]
    logliks <- ret$logliks[1:countf]
    kappas <- ret$kappas[1:countf]
    trans1 <- which(ret$trans1 == 1)
    trans2 <- which(ret$trans2 == 1)

    ## extract model parameters
    blives <- ret$blives == 1

    M = sum(blives)

    props <- ret$bprops[1:M]

    likelihoods <- matrix(likelihoods, n, Mmax)
    likelihoods <- likelihoods[,1:M, drop=F]

    mus <- matrix(ret$bmus, Mmax, p, byrow=T)
    mus <- mus[blives,, drop=F]

    likelihoods <- matrix(0, n, M)
    if (cov_type == "full") {
        covs <- array(ret$bcovs, dim=c(p, p, Mmax))
        covs <- covs[,,blives, drop=F]

        for (m in 1:M)
            likelihoods[, m] <- props[m] * mclust::dmvnorm(X, mus[m, ], covs[,,m])

    } else if (cov_type == "diag") {
        covs <- matrix(ret$bcovs, Mmax, p)
        covs <- covs[blives,, drop=F]

        for (m in 1:M)
            likelihoods[, m] <- props[m] * mclust::dmvnorm(X, mus[m, ], diag(covs[m, ]))
    }

    ## calculation posterior probabilities
    likelihoods <- pmax(likelihoods, 1e-300)
    z <- likelihoods / rowSums(likelihoods)

    clusters <- apply(z, 1, which.max)

    list(M=M, z=z, props=props, mus=mus, covs=covs, clusters=clusters,
         dl=dl, logliks=logliks, kappas=kappas, trans1=trans1, trans2=trans2)
}

dyn.load("lib/mixt4.so")
