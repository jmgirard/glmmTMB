## Downstream package compatibility audit for the ordinal family
suppressMessages(devtools::load_all("/Users/jmgirard/GitHub/glmmTMB/glmmTMB", quiet = TRUE))
data(wine, package = "ordinal")

fit  <- glmmTMB(rating ~ temp + contact + (1 | judge), data = wine,
                family = ordinal())
## psi length-1 control model (tweedie) to confirm no regression
data(Salamanders, package = "glmmTMB")
fit_tw <- glmmTMB(count ~ mined + (1 | site), family = tweedie,
                  data = Salamanders)

try_quiet <- function(label, expr) {
    res <- tryCatch(
        withCallingHandlers({
            v <- expr
            capture.output(print(v))  ## force lazy printing failures here
            v
        }, warning = function(w) invokeRestart("muffleWarning")),
        error = function(e) structure(conditionMessage(e), class = "audit_err"))
    status <- if (inherits(res, "audit_err")) paste("ERROR:", substr(res, 1, 120))
              else "OK"
    cat(sprintf("%-52s %s\n", label, status))
    invisible(res)
}

cat("\n############ ORDINAL MODEL ############\n")
cat("\n-- emmeans --\n")
library(emmeans)
em <- try_quiet("emmeans(fit, ~temp)", emmeans(fit, ~ temp))
if (!inherits(em, "audit_err")) print(em)
try_quiet("emmeans contrast (pairs)", pairs(emmeans(fit, ~ temp)))
try_quiet("emmeans type='response'", summary(emmeans(fit, ~ temp), type = "response"))

cat("\n-- car --\n")
library(car)
aa <- try_quiet("car::Anova(fit)", car::Anova(fit))
if (!inherits(aa, "audit_err")) print(aa)

cat("\n-- broom.mixed --\n")
tt <- try_quiet("broom.mixed::tidy(fit)", broom.mixed::tidy(fit))
if (!inherits(tt, "audit_err")) print(as.data.frame(tt))
try_quiet("broom.mixed::glance(fit)", broom.mixed::glance(fit))

cat("\n-- insight --\n")
try_quiet("insight::get_parameters(fit)", insight::get_parameters(fit))
try_quiet("insight::find_parameters(fit)", insight::find_parameters(fit))
try_quiet("insight::get_predicted(fit)", insight::get_predicted(fit))
try_quiet("insight::model_info(fit)$is_ordinal", insight::model_info(fit))

cat("\n-- parameters --\n")
pp <- try_quiet("parameters::model_parameters(fit)", parameters::model_parameters(fit))
if (!inherits(pp, "audit_err")) print(pp)

cat("\n-- performance --\n")
try_quiet("performance::r2(fit)", performance::r2(fit))
try_quiet("performance::check_convergence(fit)", performance::check_convergence(fit))

cat("\n-- core methods --\n")
try_quiet("confint(fit)", confint(fit))
try_quiet("profile-free vcov(fit, full=TRUE)", vcov(fit, full = TRUE))
try_quiet("drop1(fit, test='Chisq')", drop1(fit, test = "Chisq"))
try_quiet("diagnose(fit)", capture.output(diagnose(fit)))

cat("\n############ TWEEDIE CONTROL (psi length 1) ############\n")
try_quiet("emmeans(fit_tw, ~mined)", emmeans(fit_tw, ~ mined))
try_quiet("car::Anova(fit_tw)", car::Anova(fit_tw))
try_quiet("broom.mixed::tidy(fit_tw)", broom.mixed::tidy(fit_tw))
try_quiet("parameters::model_parameters(fit_tw)", parameters::model_parameters(fit_tw))
try_quiet("family_params(fit_tw)", family_params(fit_tw))

cat("\nAudit complete.\n")
