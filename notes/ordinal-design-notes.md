# Ordinal family: design notes & review crib sheet

*Branch `ordinal-family` (commit e7cc3c36), July 2026. Companion to
[GH #514](https://github.com/glmmTMB/glmmTMB/issues/514). This file is
deliberately untracked — keep it out of the PR.*

## What it is

Proportional-odds cumulative link model: `P(Y <= j) = linkinv(theta_j - eta)`,
j = 1..K-1, with ordered thresholds `theta` and the usual glmmTMB linear
predictor `eta` (fixed + random effects, any covstruct). Links: logit,
probit, cloglog. The sign convention (`theta - eta`, not `eta - theta`)
matches `ordinal::clm` and `MASS::polr`, so positive coefficients push the
response toward higher categories.

## Decisions and why

| # | Decision | Rationale | Where |
|---|----------|-----------|-------|
| 1 | Thresholds stored in `psi` (extra-family-parameter vector) | Bolker's stated preference (option 3 of his Oct 2024 comment on #514); reuses existing R<->C++ transport, priors machinery (`psi_vprior`), `map`/`start` support | `glmmTMB.cpp` PARAMETER_VECTOR(psi); `R/glmmTMB.R` psi_init block |
| 2 | Monotone parameterization `psi = c(theta[1], log(diff(theta)))` | Guarantees ordered thresholds with unconstrained optimization; same trick as `MASS::polr`. NOTE: `ordinal` pkg deliberately uses raw thresholds + Inf-on-violation for a more quadratic surface; we saw no convergence issues, but this is the most likely review discussion point. Swap = ~10 lines in C++ + start values. | `glmmTMB.cpp` theta_ord construction (search "theta_ord") |
| 3 | K inferred as `psi.size() + 1` in C++ | Avoids adding a DATA_INTEGER — no TMB data-structure change, no `up2date()` burden for stored fits | `glmmTMB.cpp` n_ord_levels |
| 4 | Fixed-effect intercept fixed to 0 via internal `map` entry | Redundant with thresholds (any intercept shift is absorbed). Mapping (vs dropping the X column) keeps model matrices, predict-on-newdata, and rank checks untouched. Cost: summary shows "(Intercept) 0 NA" row. Only applied when the user hasn't supplied their own beta map. | `R/glmmTMB.R` "intercept is redundant" block after `parameters <-` |
| 5 | `mu` redefined as expected category index E[Y] = K - sum_j P(Y<=j) | Keeps every scalar-mean downstream path working (fitted, response residuals, mu_predict). Per-category probabilities handled separately (#6). | `glmmTMB.cpp` mu override before obs loop |
| 6 | `predict(type="probs")` returns n x K matrix; REPORTed always, ADREPORTed when doPredict==1 (se.fit) | mu_predict machinery is one-scalar-per-obs; a matrix report sidesteps it. sdreport flattens column-major; predict.R reshapes with `matrix(pred, ncol=K)`. Guarded: ordinal-only, no aggregate, no cov.fit/bias.correct. | `glmmTMB.cpp` ordinal_probs block; `R/predict.R` type=="probs" branches |
| 7 | Response: ordered factor (or 1-based integer codes) -> numeric 1..K; labels kept in `modelInfo$ord_levels` | Factor conversion happens in mkTMBStruc next to the binomial factor handling. Unordered factor = warning (levels used in current order), like `clm`'s permissiveness but louder. Integer codes accepted so simulate->refit round-trips work. | `R/glmmTMB.R` ord_levels block; `finalizeTMB` modelInfo |
| 8 | psi start values: equiprobable-category thresholds `qlogis(j/K)` via `family$linkfun`, transformed to the log-diff scale | Same default as `ordinal::clm`; robust in testing | `R/glmmTMB.R` psi_init block |
| 9 | `simulate()` returns ordered factors on original levels | Friendlier for posterior predictive checks; precedent: stats::simulate for binomial factor responses. Codes are recoverable via `as.numeric()`. | `R/methods.R` simulate.glmmTMB ordinal branch |
| 10 | Dunn-Smyth residuals computed R-side from thresholds + `predict(type="link")` | The generic `dunnsmyth_resids()` is (y, mu, phi)-based; the CLM CDF needs thresholds, cleanest as a special case in `residuals.glmmTMB`. Discrete PIT: u ~ U(F(y-1), F(y)). | `R/methods.R` residuals "dunn-smyth" ordinal branch |
| 11 | No dispersion (`.noDispersionFamilies`), no zero-inflation (hard error) | sigma() returns 1; dispformula ignored like binomial/poisson. Scale effects (heteroscedastic latent SD) would later map onto `dispformula` — see DrJerryTAO's #514 comment: sd multiplier = exp(z*zeta) dividing (theta-eta). | `R/glmmTMB.R` .noDispersionFamilies; zi check in glmmTMB() |
| 12 | Family name `ordinal()` | Plan default; `cumulative()` (brms/VGAM naming) is the obvious alternative and leaves room for acat/cratio families later. Cheap to rename pre-merge; expensive after. | `R/family.R` |

## File map (all changes)

- `glmmTMB/src/glmmTMB.cpp` — enum `ordinal_family = 1100`; theta_ord + mu
  override block (before obs loop); likelihood case + SIMULATE (after bell);
  ordinal_probs REPORT/ADREPORT (after mu_predict)
- `glmmTMB/R/enum.R` — auto-generated (`make enum-update`), do not hand-edit
- `glmmTMB/R/family.R` — `ordinal()` constructor + roxygen family docs
- `glmmTMB/R/glmmTMB.R` — zi guard; ord_levels response conversion; psi
  init; intercept map; `.noDispersionFamilies`; ord_levels in return list +
  modelInfo
- `glmmTMB/R/predict.R` — type="probs" (choices, guards, return_par,
  reshape); roxygen for type
- `glmmTMB/R/methods.R` — family_params ordinal branch (threshold
  back-transform + "lo|hi" names); simulate factor mapping; residuals factor
  handling + dunn-smyth branch
- `glmmTMB/tests/testthat/test-ordinal.R` — 27 assertions
- `glmmTMB/inst/NEWS.Rd`, `NAMESPACE`, `man/nbinom2.Rd`,
  `man/predict.glmmTMB.Rd` — docs/exports

## Validation evidence

- housing data (weights): logLik/coefs/thresholds vs `clm` and `polr`
  agree to ~1e-5 (logit, probit, cloglog)
- wine data `(1|judge)`: logLik/coefs/thresholds/RE-sd vs `clmm` ~1e-5
- simulated n=2000 crossed data: glmmTMB == clmm to 4 decimals (deviation
  from *true* values is sampling noise; verified by fitting clmm to the
  same data)
- predicted probs match `clm` ~1e-5; rowSums == 1; E[Y] == probs %*% 1:K
- Dunn-Smyth residuals ~N(0,1); manual DHARMa KS p = 0.86
- No regressions: basics/families/predict/simulate/methods/disp/
  distributions/VarCorr all pass

Validation script lives in the session scratchpad (`validate_ordinal.R`);
re-create from this list if lost.

## Known issues / follow-ups

1. **DHARMa auto path broken upstream**: `DHARMa:::getObservedResponse`
   does `as.numeric(factor) - 1` (binomial assumption) while simulations
   arrive 1-based via `data.matrix()` -> off-by-one -> KS p ~ 0. Fix is a
   small DHARMa PR (teach it the ordinal family: observed =
   `as.numeric(factor)`, integerResponse = TRUE). Workaround to document:
   `createDHARMa(sapply(simulate(fit, 200), as.numeric), as.numeric(y),
   predict(fit, type="response"), integerResponse = TRUE)`.
2. Cosmetic `(Intercept) 0 NA` row in summary (could suppress in
   print.summary for ordinal).
3. Not implemented: nominal (category-specific) effects, scale effects,
   emmeans support, REML sanity check with the mapped intercept.
4. Prediction with newdata relies on the combined old+new model frame for
   factor-level consistency; document that newdata responses should use the
   same ordered levels.

## Likely review questions -> answers

- *"Why cumsum-exp when ordinal pkg avoids it?"* — polr precedent; no
  observed convergence issues with clm-style starts; happy to switch to
  Inf-on-violation (decision #2).
- *"Why is find_psi bypassed?"* — psi length is data-dependent (K-1);
  `.extraParamFamilies` is fixed-length-per-name. Generalizing find_psi to
  take data was more invasive than a one-branch special case.
- *"Does this break stored fits?"* — No: no new DATA entries, no changes to
  existing families' parameters; `up2date()` untouched.
- *"OSA residuals?"* — Not implemented; likelihood uses
  `CppAD::Integer(yobs(i))` for indexing, which is fine for taping but OSA
  oneStepGeneric would need the CDF method; Dunn-Smyth covers the
  diagnostic need meanwhile.
- *"Weights?"* — Standard per-observation likelihood weights work (used in
  the housing validation, matching clm's `weights=Freq` exactly).
