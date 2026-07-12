## Validation of glmmTMB ordinal family against MASS::polr and ordinal::clm/clmm
options(warn = 1)
suppressMessages(devtools::load_all("/Users/jmgirard/GitHub/glmmTMB/glmmTMB", quiet = TRUE))
library(MASS)
library(ordinal)

banner <- function(x) cat("\n=====", x, "=====\n")
chk <- function(label, a, b, tol = 1e-4) {
    rel <- max(abs(a - b) / pmax(abs(b), 1e-8))
    cat(sprintf("%-55s max rel diff = %.2e  %s\n", label, rel,
                if (rel < tol) "OK" else "** MISMATCH **"))
    invisible(rel < tol)
}

## ---------------------------------------------------------------
banner("1. Fixed effects, logit link: housing data (weights)")
fit_polr <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
fit_clm  <- ordinal::clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
fit_tmb  <- glmmTMB(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
                    family = ordinal())

chk("logLik vs clm",  logLik(fit_tmb), logLik(fit_clm))
chk("coefs vs clm",   fixef(fit_tmb)$cond[-1], coef(fit_clm)[-(1:2)])
chk("thresholds vs clm", unname(family_params(fit_tmb)), unname(fit_clm$alpha))
chk("coefs vs polr",  fixef(fit_tmb)$cond[-1], coef(fit_polr))
cat("threshold names:", names(family_params(fit_tmb)), "\n")

## ---------------------------------------------------------------
banner("2. Probit link")
fit_clm_p <- ordinal::clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
                          link = "probit")
fit_tmb_p <- glmmTMB(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
                     family = ordinal(link = "probit"))
chk("probit logLik vs clm", logLik(fit_tmb_p), logLik(fit_clm_p))
chk("probit coefs vs clm",  fixef(fit_tmb_p)$cond[-1], coef(fit_clm_p)[-(1:2)])

banner("3. cloglog link")
fit_clm_c <- ordinal::clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
                          link = "cloglog")
fit_tmb_c <- glmmTMB(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
                     family = ordinal(link = "cloglog"))
chk("cloglog logLik vs clm", logLik(fit_tmb_c), logLik(fit_clm_c))
chk("cloglog coefs vs clm",  fixef(fit_tmb_c)$cond[-1], coef(fit_clm_c)[-(1:2)])

## ---------------------------------------------------------------
banner("4. Mixed model: wine data, rating ~ temp + contact + (1|judge)")
data(wine, package = "ordinal")
fit_clmm <- ordinal::clmm(rating ~ temp + contact + (1 | judge), data = wine)
fit_tmbm <- glmmTMB(rating ~ temp + contact + (1 | judge), data = wine,
                    family = ordinal())
chk("clmm logLik",      logLik(fit_tmbm), logLik(fit_clmm), tol = 1e-3)
chk("clmm coefs",       fixef(fit_tmbm)$cond[-1], fit_clmm$beta, tol = 1e-2)
chk("clmm thresholds",  unname(family_params(fit_tmbm)), unname(fit_clmm$alpha), tol = 1e-2)
sd_tmb <- attr(VarCorr(fit_tmbm)$cond$judge, "stddev")
sd_clmm <- sqrt(ordinal::VarCorr(fit_clmm)$judge[1, 1])
chk("clmm RE sd",       unname(sd_tmb), unname(sd_clmm), tol = 1e-2)

## ---------------------------------------------------------------
banner("5. predict: probs sum to 1, match ordinal::clm probabilities")
pr_tmb <- predict(fit_tmb, type = "probs")
cat("prob matrix dim:", dim(pr_tmb), " colnames:", colnames(pr_tmb), "\n")
chk("probs sum to 1", rowSums(pr_tmb), rep(1, nrow(pr_tmb)))
## clm fitted probabilities for each obs category pattern
pr_clm <- predict(fit_clm, newdata = housing[, c("Infl","Type","Cont")], type = "prob")$fit
chk("probs vs clm", as.numeric(pr_tmb), as.numeric(pr_clm), tol = 1e-3)

## E[Y] response prediction consistent with probs
ey <- as.numeric(pr_tmb %*% seq_len(ncol(pr_tmb)))
chk("type='response' = E[category]", predict(fit_tmb, type = "response"), ey)

## se.fit for probs works
pr_se <- predict(fit_tmbm, type = "probs", se.fit = TRUE)
cat("probs se.fit dims:", dim(pr_se$fit), "/", dim(pr_se$se.fit),
    " (any NA se:", anyNA(pr_se$se.fit), ")\n")

## newdata prediction
nd <- wine[1:4, ]
pr_nd <- predict(fit_tmbm, newdata = nd, type = "probs")
cat("newdata probs (first 2 rows):\n"); print(round(pr_nd[1:2, ], 3))

## ---------------------------------------------------------------
banner("6. simulate + refit round trip (parameter recovery)")
set.seed(101)
n_grp <- 100; n_per <- 20
simdat <- data.frame(x = rnorm(n_grp * n_per),
                     g = factor(rep(seq_len(n_grp), each = n_per)))
b_true <- 1.2; sd_true <- 0.8; theta_true <- c(-1, 0.5, 2)
eta <- b_true * simdat$x + rnorm(n_grp)[as.integer(simdat$g)] * sd_true
u <- runif(nrow(simdat))
cum <- plogis(outer(theta_true, eta, "-"))  # 3 x n: P(Y<=j)
simdat$y <- ordered(1 + colSums(sweep(cum, 2, u, "<")), levels = 1:4)
fit_sim <- glmmTMB(y ~ x + (1 | g), data = simdat, family = ordinal())
cat(sprintf("beta: true %.2f est %.3f | RE sd: true %.2f est %.3f\n",
            b_true, fixef(fit_sim)$cond["x"], sd_true,
            attr(VarCorr(fit_sim)$cond$g, "stddev")))
chk("sim thresholds", unname(family_params(fit_sim)), theta_true, tol = 0.15)

## simulate() from fitted model returns ordered factor w/ right levels
ss <- simulate(fit_sim, nsim = 2, seed = 1)
stopifnot(is.ordered(ss[[1]]), identical(levels(ss[[1]]), levels(simdat$y)))
cat("simulate() returns ordered factor: OK; table of sim_1:\n")
print(table(ss[[1]]))

## refit to simulated response
simdat2 <- simdat; simdat2$y <- ss[[1]]
fit_sim2 <- update(fit_sim, data = simdat2)
cat(sprintf("refit beta: %.3f (orig est %.3f)\n",
            fixef(fit_sim2)$cond["x"], fixef(fit_sim)$cond["x"]))

## ---------------------------------------------------------------
banner("7. residuals + DHARMa")
r_ds <- residuals(fit_tmbm, type = "dunn-smyth")
cat("dunn-smyth residuals: mean", round(mean(r_ds), 3),
    "sd", round(sd(r_ds), 3), "(expect ~0/~1)\n")
r_resp <- residuals(fit_tmbm, type = "response")
cat("response residuals range:", round(range(r_resp), 2), "\n")

if (requireNamespace("DHARMa", quietly = TRUE)) {
    sr <- DHARMa::simulateResiduals(fit_sim, n = 200)
    ks <- suppressWarnings(DHARMa::testUniformity(sr, plot = FALSE))
    cat("DHARMa simulateResiduals ran; KS p-value =", round(ks$p.value, 3),
        "(expect > 0.05 for well-specified model)\n")
}

## ---------------------------------------------------------------
banner("8. misc: summary/print, sigma, intercept fixed, errors")
cat("sigma():", sigma(fit_tmbm), "(expect 1)\n")
cat("intercept estimate (should be absent/0):",
    fixef(fit_tmbm)$cond["(Intercept)"], "\n")
s <- capture.output(print(summary(fit_tmbm)))
cat(head(s, 25), sep = "\n")

## error paths
e1 <- tryCatch(glmmTMB(rating ~ temp, ziformula = ~1, data = wine,
                       family = ordinal()), error = conditionMessage)
cat("zi error:", e1, "\n")
e2 <- tryCatch(glmmTMB(response ~ temp, data = transform(wine, response = as.numeric(rating) - 1),
                       family = ordinal()), error = conditionMessage)
cat("bad-response error:", substr(e2, 1, 80), "\n")

cat("\nAll validation checks completed.\n")
