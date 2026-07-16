Implements the `ordinal` family (cumulative link / proportional odds) discussed in #514, following the design sketched there (thresholds in `psi`) and incorporating the feedback on my implementation comment. Closes #514.

## What's included

- `ordinal(link = "logit"|"probit"|"cloglog")` family for ordered-factor responses (1-based integer codes also accepted), supporting the full range of glmmTMB fixed/random-effect structures
- `family_params()` returns the K−1 thresholds named by adjacent level pairs (`"Low|Medium"`, ...)
- `predict(type = "probs")`: n × K per-category probability matrix, with `se.fit` support (ADREPORT)
- `type = "response"` = expected category index E[Y]; `type = "link"` = latent-scale linear predictor
- `simulate()` returns ordered factors on the original response levels
- Dunn-Smyth residuals via the discrete PIT of the cumulative-link CDF
- Tests/validation (against `MASS::polr` and `ordinal::clm`/`clmm`, all three links): coefficients, thresholds, logLik, **and standard errors** (coefficient, delta-method threshold, and predicted-probability SEs) agree to ~1e-5; random slopes + RE correlations match `clmm`; K=2 reduces exactly to binomial-logit; AIC/df bookkeeping matches `clm`/`clmm` exactly; offsets, NA handling, save/reload, K=9 with sparse extreme categories, and complete-separation behavior (diverges gracefully, like `clm`) all covered, plus prediction identities, simulate/refit round trips, error paths, and downstream (emmeans/car) regression tests

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

## Not included — proposed sequencing for follow-ups

- **Scale effects** (heteroscedastic latent scale, à la `clm(scale = ~...)`; requested by @mmrabe above): planned as the *next* PR. Implementation is small — divide `(theta − eta)` by `exp(eta_disp)` in the likelihood and let the existing `dispformula` machinery drive it — but it changes the semantics of `sigma()`/`predict(type="disp")` for this family and, more importantly, the headline combination (scale effects **with** random effects) has no frequentist reference implementation to validate against (`clm` has scale but no REs; `clmm` has REs but no scale), so it needs simulation-based validation rather than riding along here.
- **Nominal (category-specific) effects** (requested by @qdread above): deferred pending an API design decision — these need a new user-facing formula argument (`nominal = ~ ...` or similar), and per the maintainers' criteria on #982, API additions deserve their own future-proofing discussion. I'd propose opening a follow-up issue with a design sketch once this PR settles, rather than baking in an interface unilaterally.
- **Automatic DHARMa support**: blocked on a small DHARMa-side fix, not anything here — DHARMa's `getObservedResponse()` applies its binomial factor convention (`as.numeric(f) − 1`, 0-based) while simulations arrive as 1-based category codes, an off-by-one that makes `simulateResiduals()` degenerate; it also doesn't know to treat the ordinal family as an integer response. I'll submit a patch to DHARMa once this merges. Meanwhile the manual route works correctly (KS uniformity p ≈ 0.86 on a well-specified model):
  <details><summary>working <code>createDHARMa()</code> recipe</summary>

  ```r
  sr <- DHARMa::createDHARMa(
      simulatedResponse = sapply(simulate(fit, 250), as.numeric),
      observedResponse  = as.numeric(model.response(model.frame(fit))),
      fittedPredictedResponse = predict(fit, type = "response"),
      integerResponse = TRUE)
  ```
  </details>
- **OSA residuals**: feasible in principle via TMB's discrete one-step-prediction machinery (the `keep` indicator is already in place for this family), but Dunn-Smyth residuals and the DHARMa route above already cover the residual-diagnostic need, so I'd treat this as demand-driven rather than planned.

## Notes for review

- Numerical robustness: cumulative log-probabilities go through `logit_inverse_linkfun()` (accurate tails for probit/cloglog via `logit_pnorm`/`logit_invcloglog`) and the middle-category likelihood uses the branch-free identity `log(plogis(s1)−plogis(s2)) = logspace_sub(s1,s2) − log1pexp(s1) − log1pexp(s2)`, so probit/cloglog fits keep finite NLL/gradients at extreme η; thresholds are computed by exact prefix/suffix logsumexp, so no overflow when one category weight dominates (both tested)
- `confint()` reports delta-method Wald CIs for thresholds (analytic Jacobian of the softmax transform); `car::Anova` type III gives an NA row for the fixed intercept rather than an error; integer-coded responses warn that K is inferred as `max(code)`
- REML runs and behaves sanely (RE SD increases as expected) but note the thresholds (`psi`) are *not* integrated out alongside `beta` — happy to take direction on whether that's the definition you want, or whether REML should warn/error for this family
- No new TMB DATA entries and no changes to existing families' parameter structures, so no `up2date()` burden
- `enum.R` regenerated via `make enum-update`; docs via roxygen (kept at RoxygenNote 7.3.3 formatting where possible)
- Full test suite + validation results in the PR checks; happy to walk through any part of the diff

*Disclosure, as on #514: developed with substantial help from an AI assistant (Claude), with design directed by the #514 discussion and every numeric claim validated against `MASS::polr` / `ordinal::clm` / `clmm`.*
