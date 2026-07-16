## Pre-PR correctness gauntlet for the ordinal family.
## Focus: everything NOT covered by prior validation (SEs, edge cases,
## bookkeeping, stress).
options(warn = 1)
suppressMessages(devtools::load_all("/Users/jmgirard/GitHub/glmmTMB/glmmTMB", quiet = TRUE))
suppressMessages(library(ordinal))
data(wine, package = "ordinal"); data(housing, package = "MASS")

chk <- function(label, a, b, tol = 1e-3) {
    rel <- max(abs(a - b) / pmax(abs(b), 1e-8))
    cat(sprintf("%-58s max rel diff = %.2e  %s\n", label, rel,
                if (rel < tol) "OK" else "** MISMATCH **"))
}
hdr <- function(x) cat("\n=====", x, "=====\n")

## ---------------------------------------------------------------
hdr("A1. Coefficient SEs vs clm (fixed effects, housing)")
f_clm <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
f_tmb <- glmmTMB(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
                 family = ordinal())
se_clm <- summary(f_clm)$coefficients[, "Std. Error"]
se_tmb <- summary(f_tmb)$coefficients$cond[, "Std. Error"]
chk("beta SEs vs clm", se_tmb, se_clm[-(1:2)][names(se_tmb)])

hdr("A2. Threshold SEs vs clm (delta method through softmax)")
Vpsi <- vcov(f_tmb, full = TRUE)
## psi rows are labeled with family_params() names (threshold labels)
psi_i <- match(names(family_params(f_tmb)), rownames(Vpsi))
Vpsi <- Vpsi[psi_i, psi_i]
psi <- f_tmb$fit$parfull[names(f_tmb$fit$parfull) == "psi"]
th_fun <- function(p) { w <- exp(c(p, 0)); qlogis(cumsum(w / sum(w))[seq_along(p)]) }
J <- numDeriv::jacobian(th_fun, psi)
se_th_tmb <- sqrt(diag(J %*% Vpsi %*% t(J)))
se_th_clm <- summary(f_clm)$coefficients[1:2, "Std. Error"]
chk("threshold SEs vs clm", se_th_tmb, unname(se_th_clm))

hdr("A3. Coefficient SEs vs clmm (mixed, wine)")
f_clmm <- clmm(rating ~ temp + contact + (1 | judge), data = wine)
f_tmbm <- glmmTMB(rating ~ temp + contact + (1 | judge), data = wine,
                  family = ordinal())
se_clmm <- summary(f_clmm)$coefficients[, "Std. Error"]
se_tmbm <- summary(f_tmbm)$coefficients$cond[, "Std. Error"]
chk("mixed beta SEs vs clmm", se_tmbm, se_clmm[names(se_tmbm)], tol = 5e-3)

hdr("A4. Predicted-probability SEs vs clm")
pr_tmb <- predict(f_tmb, type = "probs", se.fit = TRUE)
pr_clm <- predict(f_clm, newdata = housing[, c("Infl","Type","Cont")],
                  type = "prob", se.fit = TRUE)
chk("prob SEs vs clm", as.numeric(pr_tmb$se.fit), as.numeric(pr_clm$se.fit),
    tol = 5e-3)

## ---------------------------------------------------------------
hdr("B. Random slopes vs clmm (simulated)")
set.seed(42)
ng <- 60; np <- 25
d <- data.frame(x = rnorm(ng * np), g = factor(rep(1:ng, each = np)))
b0 <- rnorm(ng, 0, 0.8); b1 <- rnorm(ng, 0, 0.5)
eta <- (1 + b1[as.integer(d$g)]) * d$x + b0[as.integer(d$g)]
u <- runif(nrow(d))
d$y <- ordered(1 + colSums(sweep(plogis(outer(c(-1, 0.6, 2), eta, "-")), 2, u, "<")))
f1 <- glmmTMB(y ~ x + (1 + x | g), data = d, family = ordinal())
f2 <- clmm(y ~ x + (1 + x | g), data = d)
chk("slope model logLik", c(logLik(f1)), c(logLik(f2)), tol = 1e-4)
chk("slope model beta", unname(fixef(f1)$cond["x"]), unname(f2$beta), tol = 5e-3)
vc1 <- VarCorr(f1)$cond$g; vc2 <- ordinal::VarCorr(f2)$g
chk("slope model RE SDs", attr(vc1, "stddev"), sqrt(diag(vc2)), tol = 1e-2)
chk("slope model RE corr", attr(vc1, "correlation")[1,2],
    cov2cor(vc2)[1,2], tol = 2e-2)

## ---------------------------------------------------------------
hdr("C. K=2 reduces exactly to binomial logit")
d2 <- data.frame(x = rnorm(500))
d2$y01 <- rbinom(500, 1, plogis(0.5 + d2$x))
d2$yord <- ordered(d2$y01 + 1)
fb <- glmmTMB(y01 ~ x, family = binomial, data = d2)
fo <- glmmTMB(yord ~ x, family = ordinal(), data = d2)
chk("K=2 logLik == binomial", c(logLik(fo)), c(logLik(fb)), tol = 1e-6)
chk("K=2 beta == binomial", unname(fixef(fo)$cond["x"]),
    unname(fixef(fb)$cond["x"]), tol = 1e-4)
chk("K=2 threshold == -intercept", unname(family_params(fo)),
    -unname(fixef(fb)$cond["(Intercept)"]), tol = 1e-4)

## ---------------------------------------------------------------
hdr("D. AIC / df bookkeeping vs clm & clmm")
chk("AIC vs clm", AIC(f_tmb), AIC(f_clm), tol = 1e-6)
chk("AIC vs clmm", AIC(f_tmbm), AIC(f_clmm), tol = 1e-5)
chk("logLik df vs clm", attr(logLik(f_tmb), "df"), f_clm$edf, tol = 1e-9)

## ---------------------------------------------------------------
hdr("E. Offset support vs clm")
housing$off <- rep(c(0.1, -0.2, 0), length.out = nrow(housing))
fo_clm <- clm(Sat ~ Infl + offset(off), weights = Freq, data = housing)
fo_tmb <- glmmTMB(Sat ~ Infl + offset(off), weights = Freq, data = housing,
                  family = ordinal())
chk("offset logLik vs clm", c(logLik(fo_tmb)), c(logLik(fo_clm)), tol = 1e-6)
chk("offset coefs vs clm", unname(fixef(fo_tmb)$cond[-1]),
    unname(coef(fo_clm)[-(1:2)]), tol = 1e-4)

## ---------------------------------------------------------------
hdr("F. NA handling & prediction paths")
wine_na <- wine; wine_na$rating[c(3, 10)] <- NA; wine_na$temp[5] <- NA
f_na <- glmmTMB(rating ~ temp + contact + (1 | judge), data = wine_na,
                family = ordinal())
cat("NA fit ok:", f_na$fit$convergence == 0, " nobs =", nobs(f_na), "(expect 69)\n")
p_na <- predict(f_na, newdata = wine_na[1:6, ], type = "probs")
cat("predict newdata w/ NA predictor: rows =", nrow(p_na),
    "; NA row propagated:", anyNA(p_na[5, ]), "\n")
p_pop <- predict(f_tmbm, newdata = data.frame(temp = "warm", contact = "no",
                                              judge = NA),
                 type = "probs", re.form = NA, allow.new.levels = TRUE)
cat("population-level newdata predict: ", round(p_pop, 3), "\n")
r_na <- residuals(f_na, type = "response")
cat("residuals length with na.action:", length(r_na), "(expect 72 w/ NA)\n")

## ---------------------------------------------------------------
hdr("G. save / reload round trip")
tf <- tempfile(fileext = ".rds")
saveRDS(f_tmbm, tf)
f_re <- up2date(readRDS(tf))
chk("reloaded predict == original",
    predict(f_re, type = "response"), predict(f_tmbm, type = "response"),
    tol = 1e-10)
chk("reloaded probs == original",
    as.numeric(predict(f_re, type = "probs")),
    as.numeric(predict(f_tmbm, type = "probs")), tol = 1e-10)

## ---------------------------------------------------------------
hdr("H. REML behavior")
fr <- try(glmmTMB(rating ~ temp + contact + (1 | judge), data = wine,
                  family = ordinal(), REML = TRUE), silent = TRUE)
if (inherits(fr, "try-error")) {
    cat("REML errors:", substr(attr(fr, "condition")$message, 1, 80), "\n")
} else {
    cat("REML runs; logLik =", c(logLik(fr)),
        "; beta =", round(unname(fixef(fr)$cond[-1]), 3),
        "; RE sd =", round(attr(VarCorr(fr)$cond$judge, "stddev"), 3), "\n")
}

## ---------------------------------------------------------------
hdr("I. Stress: K=9 with sparse extremes; separation")
set.seed(7)
n <- 3000
xs <- rnorm(n)
ths <- qlogis(c(0.002, 0.01, seq(0.1, 0.9, length.out = 4), 0.99, 0.998))
us <- runif(n)
ys <- ordered(1 + colSums(sweep(plogis(outer(ths, 1.5 * xs, "-")), 2, us, "<")),
              levels = 1:9)
cat("K=9 category counts:", table(ys), "\n")
fs <- glmmTMB(ys ~ xs, family = ordinal())
cat("K=9 sparse fit: convergence =", fs$fit$convergence,
    "; max |grad| =", format(max(abs(fs$sdr$gradient.fixed)), digits = 3), "\n")
fs_clm <- clm(ys ~ xs)
chk("K=9 logLik vs clm", c(logLik(fs)), c(logLik(fs_clm)), tol = 1e-6)
chk("K=9 thresholds vs clm", unname(family_params(fs)),
    unname(fs_clm$alpha), tol = 1e-3)

## complete separation: estimates should blow up gracefully (like clm)
dsep <- data.frame(x = c(rep(0, 20), rep(1, 20)),
                   y = ordered(c(rep(1, 20), rep(3, 20)), levels = 1:3))
fsep <- suppressWarnings(try(glmmTMB(y ~ x, data = dsep, family = ordinal()),
                             silent = TRUE))
if (inherits(fsep, "try-error")) {
    cat("separation: hard error (clm-comparable?)\n")
} else {
    cat("separation: fits with |beta| =",
        round(abs(unname(fixef(fsep)$cond["x"])), 1),
        "(expect large), SE =",
        format(summary(fsep)$coefficients$cond["x","Std. Error"], digits=2), "\n")
}

## ---------------------------------------------------------------
hdr("J. Priors on psi (thresholds)")
fp <- try(glmmTMB(rating ~ temp + contact + (1 | judge), data = wine,
                  family = ordinal(),
                  priors = data.frame(prior = "normal(0, 3)", class = "psi",
                                      coef = "")), silent = TRUE)
if (inherits(fp, "try-error")) {
    cat("psi prior errors:", substr(attr(fp, "condition")$message, 1, 100), "\n")
} else {
    cat("psi prior runs; thresholds shrink toward equiprobable:\n")
    print(round(rbind(no_prior = family_params(f_tmbm),
                      prior = family_params(fp)), 3))
}

cat("\nGauntlet complete.\n")
