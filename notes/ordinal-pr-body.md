Implements the `ordinal` family (cumulative link / proportional odds) discussed in #514, following the design sketched there (thresholds in `psi`) and incorporating the feedback on my implementation comment. Closes #514.

## What's included

- `ordinal(link = "logit"|"probit"|"cloglog")` family for ordered-factor responses (1-based integer codes also accepted), supporting the full range of glmmTMB fixed/random-effect structures
- `family_params()` returns the K−1 thresholds named by adjacent level pairs (`"Low|Medium"`, ...)
- `predict(type = "probs")`: n × K per-category probability matrix, with `se.fit` support (ADREPORT)
- `type = "response"` = expected category index E[Y]; `type = "link"` = latent-scale linear predictor
- `simulate()` returns ordered factors on the original response levels
- Dunn-Smyth residuals via the discrete PIT of the cumulative-link CDF
- Tests: validation against `MASS::polr` and `ordinal::clm`/`clmm` (coefficients, thresholds, logLik to ~1e-5, all three links), prediction identities, simulate/refit round trip, error paths, and downstream (emmeans/car) regression tests

## Responses to the review points raised on #514

**1. `length(psi) > 2` and downstream packages.** I ran a downstream audit on an ordinal fit (wine data, `rating ~ temp + contact + (1|judge)`) and a tweedie control (length-1 psi): emmeans (EMMs, contrasts, `type="response"`), `car::Anova`, `broom.mixed::tidy/glance`, `insight::get_parameters/find_parameters/get_predicted/model_info`, `parameters::model_parameters`, `performance::r2/check_convergence`, `confint`, `drop1`, `diagnose` — all pass; the tweedie control is unaffected. (Two of these needed the fix described under point 3.) The easystats folks have signalled on #514 that any remaining special-casing on their side is easy.

**2. Threshold parameterization.** I took the softmax suggestion (Koslik et al 2025) seriously and ran a controlled comparison — same likelihood, same optimizer (nlminb), starts identical on the threshold scale, 50 random starts × 8 conditions (n ∈ {200, 2000} × K ∈ {4, 7} × balanced/sparse-extreme-categories):

| parameterization | found optimum | mean iterations |
|---|---|---|
| softmax (Koslik-style) | 100% | 24.5 |
| raw + Inf-on-violation (`ordinal`-style) | 100% | 27.0 |
| cumsum-exp (`polr`-style) | 100% | 36.1 |

All three are equally *robust*, but softmax is the most efficient, decisively so when extreme categories are sparse (~19–29 vs ~64–72 iterations for cumsum-exp at n=200). **This PR therefore uses the softmax parameterization**: `psi` holds log-weights of the K baseline category probabilities (last fixed to 0) and `theta = qlogis(cumsum(softmax(c(psi, 0))))`. Two side benefits: the equiprobable default start is exactly `psi = 0`, and `psi` has a direct interpretation (baseline log probability ratios vs. the last category). The comparison script is available on request.

**3. Internal `map` for the intercept.** Verified that `modelInfo$map` records only the *user's* map (`NULL` for a plain ordinal fit); the internal map lives only in `obj$env$map`, and user-supplied psi maps compose correctly with it (tested). The audit surfaced that emmeans and `car::Anova` mishandle map-fixed coefficients — for *any* mapped model, not just ordinal (`vcov` drops/NA-fills mapped entries while `fixef` keeps them, so emmeans hits a dimension mismatch and Anova's hypothesis algebra propagates `0*NA = NA` into blank chisq values). Added a small `pad_mapped_vcov()` helper that treats map-fixed coefficients as known constants (zero variance), used by `emm_basis.glmmTMB` and `Anova.glmmTMB`. Regression tests assert EMM contrasts equal the fixed-effect coefficients and Anova chisq equals z². The cosmetic `(Intercept) 0 NA` summary row is now suppressed (only when the map was internal).

**4. E[Y] as the reported mean.** Kept, with the reasoning documented: some scalar per-observation summary is needed for `fitted()`/response residuals/`mu_predict`, and E[Y] is the standard choice (it's what DHARMa-style diagnostics want as `fittedPredictedResponse`, and the alternatives — modal or median category — discard more information). The "real" predictions are `type = "probs"`; docs point users there.

## Not included (possible follow-ups)

- nominal (category-specific) and scale effects (scale plausibly maps onto `dispformula` later)
- automatic DHARMa support — blocked on a small DHARMa-side fix (its `getObservedResponse()` applies a binomial factor convention, `as.numeric(f) − 1`, while simulations arrive 1-based); the manual `createDHARMa()` route gives uniform residuals (KS p = 0.86) and I'm happy to send DHARMa a patch
- OSA residuals

## Notes for review

- Numerical robustness: cumulative log-probabilities go through `logit_inverse_linkfun()` (accurate tails for probit/cloglog via `logit_pnorm`/`logit_invcloglog`) and the middle-category likelihood uses the branch-free identity `log(plogis(s1)−plogis(s2)) = logspace_sub(s1,s2) − log1pexp(s1) − log1pexp(s2)`, so probit/cloglog fits keep finite NLL/gradients at extreme η (tested); thresholds are computed by exact prefix/suffix logsumexp, so no overflow when one category weight dominates (tested)
- Beyond estimate validation, the test/validation battery also covers: coefficient/threshold/predicted-probability **standard errors** vs `clm`/`clmm` (~1e-5, thresholds via delta method), random slopes + RE correlation vs `clmm`, exact K=2 ≡ binomial-logit reduction, exact AIC/df bookkeeping vs `clm`/`clmm`, offsets, NA handling, save/reload, K=9 with sparse extreme categories, and complete-separation behavior (diverges gracefully like `clm`)
- `confint()` reports delta-method Wald CIs for thresholds (analytic Jacobian of the softmax transform); `car::Anova` type III gives an NA row for the fixed intercept rather than an error; integer-coded responses warn that K is inferred as `max(code)`
- REML runs and behaves sanely (RE SD increases as expected) but note the thresholds (`psi`) are *not* integrated out alongside `beta` — happy to take direction on whether that's the definition you want, or whether REML should warn/error for this family
- No new TMB DATA entries and no changes to existing families' parameter structures, so no `up2date()` burden
- `enum.R` regenerated via `make enum-update`; docs via roxygen (kept at RoxygenNote 7.3.3 formatting where possible)
- Full test suite + validation results in the PR checks; happy to walk through any part of the diff

*Disclosure, as on #514: developed with substantial help from an AI assistant (Claude), with design directed by the #514 discussion and every numeric claim validated against `MASS::polr` / `ordinal::clm` / `clmm`.*
