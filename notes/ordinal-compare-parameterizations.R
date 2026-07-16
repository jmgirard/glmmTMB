## Empirical comparison of threshold parameterizations for cumulative link
## models, responding to Bolker's question about a softmax approach
## (Koslik et al 2025, arXiv:2511.17071) vs cumsum-exp vs raw thresholds.
##
## All three share the same likelihood and optimizer (nlminb, numeric
## gradients); they differ only in how unconstrained parameters map to
## ordered thresholds. Starts are identical *on the threshold scale* and
## mapped exactly into each parameterization.

set.seed(2026)

theta_to_par <- function(theta, param) {
    switch(param,
           raw       = theta,
           cumsumexp = c(theta[1], log(diff(theta))),
           softmax   = {                     # theta -> cumulative probs -> simplex -> log-ratios
               p <- diff(c(0, plogis(theta), 1))
               log(p[-length(p)] / p[length(p)])
           })
}
par_to_theta <- function(psi, param, K) {
    switch(param,
           raw       = psi,
           cumsumexp = cumsum(c(psi[1], exp(psi[-1]))),
           softmax   = {
               p <- exp(c(psi, 0)); p <- p / sum(p)
               qlogis(cumsum(p)[seq_len(K - 1)])
           })
}

nll_factory <- function(y, X, param) {
    K <- max(y); nb <- ncol(X); nt <- K - 1
    function(par) {
        theta <- par_to_theta(par[seq_len(nt)], param, K)
        if (param == "raw" && is.unsorted(theta)) return(Inf)
        eta <- drop(X %*% par[nt + seq_len(nb)])
        cdf <- function(j)
            ifelse(j <= 0, 0, ifelse(j >= K, 1,
                   plogis(theta[pmax(pmin(j, K - 1), 1)] - eta)))
        -sum(log(pmax(cdf(y) - cdf(y - 1), 1e-300)))
    }
}

simdata <- function(n, K, beta = 1, sparse = FALSE) {
    x <- rnorm(n)
    ## sparse: pile mass in middle categories, starving the extremes
    theta <- if (sparse) qlogis(c(0.005, seq(0.15, 0.85, length.out = K - 3), 0.995))
             else qlogis(seq_len(K - 1) / K)
    eta <- beta * x
    u <- runif(n)
    y <- 1 + colSums(sweep(plogis(outer(theta, eta, "-")), 2, u, "<"))
    list(y = y, X = cbind(x), theta = theta)
}

run_cell <- function(n, K, sparse, nrep = 50) {
    out <- list()
    for (r in seq_len(nrep)) {
        d <- simdata(n, K, sparse = sparse)
        ## random start on the threshold scale (same for all params):
        ## sorted normal draws, sd 2 -- deliberately rough
        theta0 <- sort(rnorm(K - 1, 0, 2))
        ## enforce minimum gap so cumsumexp/softmax mappings are finite
        theta0 <- theta0 + cumsum(c(0, pmax(0, 0.05 - diff(theta0))))
        best <- Inf
        cell <- list()
        for (pm in c("raw", "cumsumexp", "softmax")) {
            f <- nll_factory(d$y, d$X, pm)
            p0 <- c(theta_to_par(theta0, pm), 0)
            fit <- suppressWarnings(
                try(nlminb(p0, f, control = list(iter.max = 500, eval.max = 1000)),
                    silent = TRUE))
            if (inherits(fit, "try-error")) {
                cell[[pm]] <- data.frame(param = pm, conv = FALSE, nll = NA,
                                         iter = NA, feval = NA)
            } else {
                cell[[pm]] <- data.frame(param = pm, conv = fit$convergence == 0,
                                         nll = fit$objective,
                                         iter = fit$iterations,
                                         feval = fit$evaluations[1])
                best <- min(best, fit$objective, na.rm = TRUE)
            }
        }
        cc <- do.call(rbind, cell)
        cc$hit_opt <- !is.na(cc$nll) & (cc$nll - best) < 1e-4
        cc$rep <- r
        out[[r]] <- cc
    }
    cbind(n = n, K = K, sparse = sparse, do.call(rbind, out))
}

grid <- expand.grid(n = c(200, 2000), K = c(4, 7), sparse = c(FALSE, TRUE))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    run_cell(grid$n[i], grid$K[i], grid$sparse[i])))

agg <- aggregate(cbind(hit_opt = hit_opt, conv = conv, iter = iter, feval = feval)
                 ~ param + n + K + sparse, data = res, FUN = mean)
agg$iter <- round(agg$iter, 1); agg$feval <- round(agg$feval, 1)
agg <- agg[order(agg$sparse, agg$K, agg$n, agg$param), ]
print(agg, row.names = FALSE)

cat("\n--- overall ---\n")
ov <- aggregate(cbind(hit_opt, conv, iter, feval) ~ param, data = res, FUN = mean)
ov$iter <- round(ov$iter, 1); ov$feval <- round(ov$feval, 1)
print(ov, row.names = FALSE)
