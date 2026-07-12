Following up on this thread: I've implemented a working version of an `ordinal` family (proportional-odds cumulative link) for glmmTMB, following the design @bbolker sketched in the [Oct 2024 comment](https://github.com/glmmTMB/glmmTMB/issues/514#issuecomment-2430571238) above (thresholds in `psi`) and the RTMB prototype in bbolker/bbmisc.

Branch: https://github.com/jmgirard/glmmTMB/tree/ordinal-family ([diff vs. master](https://github.com/glmmTMB/glmmTMB/compare/master...jmgirard:glmmTMB:ordinal-family), single commit, +374/−22)

**Usage:**

```r
fit <- glmmTMB(rating ~ temp + contact + (1 | judge),
               data = wine, family = ordinal())  # logit/probit/cloglog
family_params(fit)             # thresholds, named "1|2", "2|3", ...
predict(fit, type = "probs")   # n x K probability matrix (se.fit supported)
```

**Validation** (all three links):

- `ordinal::clm` / `MASS::polr` fixed-effects (housing data): logLik, coefficients, and thresholds agree to ~1e-5 or better
- `ordinal::clmm` (wine data, `(1|judge)`): logLik, coefs, thresholds, RE sd agree to ~1e-5
- predicted probabilities match `clm` to ~1e-5; simulate-and-recover on n=2000 with random intercepts matches `clmm` to 4 decimals
- 27 new testthat assertions pass; no regressions in the existing suite (basics/families/predict/simulate/methods/disp)

**Design choices, per the Oct 2024 comment — happy to change any of these:**

1. Thresholds live in `psi` (option 3 from that comment), length K−1 (data-dependent, so `find_psi()` is bypassed for this family).
2. Monotonicity is enforced by parameterization: `psi = c(theta[1], log(diff(theta)))`, MASS::polr-style. I know the `ordinal` package deliberately avoids this for convergence reasons; starting values are equiprobable-category thresholds (as in `clm`) and I haven't seen convergence problems in testing, but this is easy to swap for the NaN-on-violation approach if you prefer.
3. The fixed-effect intercept is fixed to 0 via an internal `map` entry (absorbed into the thresholds) rather than dropping the X column — this keeps the model matrices and predict machinery untouched. Cosmetic cost: summary shows an `(Intercept) 0 NA` row.
4. `mu` (and hence `type="response"`, `fitted()`, response residuals) is the expected category index E[Y]; per-category probabilities are REPORTed and exposed via a new `predict(type = "probs")` (ADREPORTed when `se.fit = TRUE`).
5. The response must be an ordered factor (or 1-based integer codes); level labels are stored in `modelInfo$ord_levels`; `simulate()` returns ordered factors on the original levels.
6. Dunn-Smyth residuals are implemented via the discrete PIT of the CLM CDF, computed R-side from the thresholds and latent predictor.

No changes to the TMB data structures (so no `up2date()` burden) and no new family-specific arguments in the API.

**Known limitations / follow-ups:** no nominal (category-specific) or scale effects yet (scale would plausibly map onto `dispformula` later, per @DrJerryTAO's analysis above); DHARMa's automatic path needs a small DHARMa-side fix (its `getObservedResponse()` does `as.numeric(factor) - 1`, a binomial assumption, while simulations come through 1-based — the manual `createDHARMa()` route gives uniform residuals, KS p = 0.86); no emmeans support yet.

Would you be open to a PR? If so, which of the design points above would you like handled differently first?
