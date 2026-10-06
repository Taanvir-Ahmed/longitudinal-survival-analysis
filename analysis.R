## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 8, fig.height = 4.8, dpi = 150, fig.align = "center"
)


## ----packages-----------------------------------------------------------------
library(tidyverse)   # data handling and plotting
library(afex)        # repeated-measures ANOVA
library(emmeans)     # estimated marginal means
library(lme4)        # linear mixed models
library(lmerTest)    # p-values for mixed-model fixed effects
library(geepack)     # generalized estimating equations
library(survival)    # Kaplan-Meier and Cox models
library(broom)
library(broom.mixed)
library(knitr)

theme_set(theme_minimal(base_size = 12))
trt_cols   <- c("Standard" = "#2C7FB8", "New" = "#D95F0E")
dir.create("figures", showWarnings = FALSE)
save_fig <- function(p, name, w = 8, h = 4.8) {
  ggsave(file.path("figures", name), p, width = w, height = h, dpi = 150, bg = "white")
  p
}


## ----read---------------------------------------------------------------------
data_wide <- read.table("data/extra_data.txt", header = TRUE) |>
  mutate(trt = factor(treatment, levels = 0:1, labels = c("Standard", "New")))

dim(data_wide)
anyNA(data_wide)


## ----long---------------------------------------------------------------------
data_long <- data_wide |>
  pivot_longer(starts_with("bio_"), names_to = "visit", values_to = "biomarker") |>
  mutate(
    time_numeric = as.numeric(sub("bio_", "", visit)),
    time_factor  = factor(time_numeric, levels = c(0, 3, 6, 12)),
    id           = factor(id),
    age_c        = age - mean(data_wide$age)   # centred age, used in the mixed models
  )

table(data_long$time_factor)
stopifnot(nrow(data_long) == 240 * 4)


## ----baseline-----------------------------------------------------------------
data_wide |>
  group_by(Treatment = trt) |>
  summarise(
    n            = n(),
    `Age, mean (SD)`        = sprintf("%.1f (%.1f)", mean(age), sd(age)),
    `Male, %`               = sprintf("%.0f", 100 * mean(sex)),
    `Baseline biomarker, mean (SD)` = sprintf("%.1f (%.1f)", mean(bio_0), sd(bio_0)),
    Deaths       = sum(death)
  ) |>
  kable()

t.test(bio_0 ~ trt, data = data_wide)


## ----fig-trajectories, fig.height = 4.5---------------------------------------
mean_profiles <- data_long |>
  group_by(trt, time_numeric) |>
  summarise(mean = mean(biomarker), se = sd(biomarker) / sqrt(n()), .groups = "drop")

p_traj <- ggplot(data_long, aes(time_numeric, biomarker)) +
  geom_line(aes(group = id, colour = trt), alpha = 0.18, linewidth = 0.35) +
  geom_line(data = mean_profiles, aes(y = mean, colour = trt), linewidth = 1.4) +
  geom_point(data = mean_profiles, aes(y = mean, colour = trt), size = 2.6) +
  scale_colour_manual(values = trt_cols) +
  scale_x_continuous(breaks = c(0, 3, 6, 12)) +
  labs(x = "Months since baseline", y = "Biomarker", colour = "Treatment",
       title = "Individual biomarker trajectories with group means") +
  theme(legend.position = "top")

save_fig(p_traj, "01_trajectories.png")


## ----mean-table---------------------------------------------------------------
mean_profiles |>
  select(-se) |>
  pivot_wider(names_from = time_numeric, values_from = mean, names_prefix = "Month ") |>
  mutate(`Change 0→12` = `Month 12` - `Month 0`) |>
  kable(digits = 1)


## ----rmanova------------------------------------------------------------------
anova_rm <- aov_car(
  biomarker ~ trt * time_factor + Error(id / time_factor),
  data = data_long,
  anova_table = list(correction = "none", es = "pes")
)
anova_rm


## ----rmanova-emm--------------------------------------------------------------
emmeans(anova_rm, ~ trt | time_factor) |> as_tibble() |>
  select(time_factor, trt, emmean, SE) |> kable(digits = 2)


## ----sphericity---------------------------------------------------------------
summary(anova_rm)$sphericity.tests
summary(anova_rm)$pval.adjustments


## ----gg-corrected-------------------------------------------------------------
nice(anova_rm, correction = "GG")


## ----cov-structure------------------------------------------------------------
wide_bio <- data_wide |> select(bio_0, bio_3, bio_6, bio_12)
round(apply(wide_bio, 2, var), 1)       # variances by visit
round(cor(wide_bio), 2)                 # correlations between visits


## ----fig-rm-resid, fig.height = 3.6-------------------------------------------
rm_resid <- tibble(resid = residuals(anova_rm))

p_qq <- ggplot(rm_resid, aes(sample = resid)) +
  stat_qq(alpha = 0.4, size = 0.9) + stat_qq_line(colour = "#D95F0E") +
  labs(title = "Normal Q-Q plot of RM-ANOVA residuals",
       x = "Theoretical quantiles", y = "Residual")
save_fig(p_qq, "02_rmanova_qq.png", h = 3.6)


## ----lmm-ri-------------------------------------------------------------------
model_ri <- lmer(
  biomarker ~ time_numeric * treatment + age_c + sex + (1 | id),
  data = data_long, REML = TRUE
)

fixef_tab <- tidy(model_ri, effects = "fixed", conf.int = TRUE) |>
  select(term, estimate, std.error, conf.low, conf.high, p.value)
kable(fixef_tab, digits = c(0, 3, 3, 3, 3, 4),
      col.names = c("Term", "Estimate", "SE", "95% CI low", "95% CI high", "p"))


## ----lmm-slopes---------------------------------------------------------------
b <- fixef(model_ri)
slopes <- c(Standard = b[["time_numeric"]],
            New      = b[["time_numeric"]] + b[["time_numeric:treatment"]])
round(slopes, 3)


## ----icc----------------------------------------------------------------------
vc <- as.data.frame(VarCorr(model_ri))
sigma2_b <- vc$vcov[vc$grp == "id"]
sigma2_e <- vc$vcov[vc$grp == "Residual"]
icc <- sigma2_b / (sigma2_b + sigma2_e)

tibble(Component = c("Between-patient (random intercept)", "Within-patient (residual)"),
       Variance  = c(sigma2_b, sigma2_e),
       SD        = sqrt(c(sigma2_b, sigma2_e))) |>
  kable(digits = 2)
round(icc, 3)


## ----fig-ranef, fig.height = 3.4----------------------------------------------
re <- ranef(model_ri)$id |> rownames_to_column("id") |> rename(b0 = `(Intercept)`)
p_re <- ggplot(re, aes(sample = b0)) +
  stat_qq(alpha = 0.5, size = 1) + stat_qq_line(colour = "#D95F0E") +
  labs(title = "Predicted random intercepts (Q-Q)", x = "Theoretical quantiles", y = "b0")
save_fig(p_re, "03_random_intercepts_qq.png", h = 3.4)


## ----lmm-rs-------------------------------------------------------------------
model_rs <- lmer(
  biomarker ~ time_numeric * treatment + age_c + sex + (time_numeric | id),
  data = data_long, REML = TRUE
)
VarCorr(model_rs)
isSingular(model_rs)


## ----lrt----------------------------------------------------------------------
lrt <- anova(model_ri, model_rs, refit = FALSE)
lrt

# The null hypothesis puts the slope variance on the boundary of its space
# (sigma^2_slope = 0), so the naive chi-square(2) p-value is conservative.
# The correct reference is a 50:50 mixture of chi-square(1) and chi-square(2).
chisq <- lrt$Chisq[2]
p_mixture <- 0.5 * pchisq(chisq, 1, lower.tail = FALSE) +
             0.5 * pchisq(chisq, 2, lower.tail = FALSE)
round(c(LRT = chisq, p_naive = lrt$`Pr(>Chisq)`[2], p_mixture = p_mixture), 3)


## ----fig-lmm-resid, fig.height = 3.6------------------------------------------
p_fit <- tibble(fitted = fitted(model_ri), resid = resid(model_ri)) |>
  ggplot(aes(fitted, resid)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(se = FALSE, colour = "#D95F0E", method = "loess") +
  labs(title = "Random-intercept model: residuals vs fitted", x = "Fitted", y = "Residual")
save_fig(p_fit, "04_lmm_residuals.png", h = 3.6)


## ----gee----------------------------------------------------------------------
data_gee <- data_long |> arrange(id, time_numeric)

gee_ind  <- geeglm(biomarker ~ time_numeric * treatment, id = id, data = data_gee,
                   family = gaussian, corstr = "independence")
gee_exch <- geeglm(biomarker ~ time_numeric * treatment, id = id, data = data_gee,
                   family = gaussian, corstr = "exchangeable")
gee_ar1  <- geeglm(biomarker ~ time_numeric * treatment, id = id, data = data_gee,
                   family = gaussian, corstr = "ar1")

gee_tab <- map_dfr(
  list(Independence = gee_ind, Exchangeable = gee_exch, `AR(1)` = gee_ar1),
  ~ tidy(.x, conf.int = TRUE), .id = "Working correlation"
) |>
  select(`Working correlation`, term, estimate, std.error, conf.low, conf.high)
kable(gee_tab, digits = 3)

c(exchangeable_alpha = gee_exch$geese$alpha, ar1_alpha = gee_ar1$geese$alpha) |> round(3)
sapply(list(Independence = gee_ind, Exchangeable = gee_exch, `AR(1)` = gee_ar1),
       function(m) QIC(m)[c("QIC", "CIC")]) |> round(1)


## ----gee-vs-lmm---------------------------------------------------------------
gee_exch_adj <- geeglm(biomarker ~ time_numeric * treatment + age_c + sex, id = id,
                       data = data_gee, family = gaussian, corstr = "exchangeable")

cmp <- full_join(
  tidy(model_ri, effects = "fixed") |> select(term, LMM = estimate, LMM_SE = std.error),
  tidy(gee_exch_adj) |> select(term, GEE = estimate, GEE_robust_SE = std.error),
  by = "term"
)
kable(cmp, digits = 3)


## ----km-----------------------------------------------------------------------
km_fit <- survfit(Surv(follow, death) ~ trt, data = data_wide)
km_fit


## ----fig-km, fig.height = 4.6-------------------------------------------------
km_df <- tidy(km_fit) |>
  mutate(trt = factor(sub("trt=", "", strata), levels = c("Standard", "New")))
km_start <- km_df |> distinct(trt) |>
  mutate(time = 0, estimate = 1, conf.low = 1, conf.high = 1, n.censor = 0)
km_df <- bind_rows(km_start, km_df) |> arrange(trt, time) |>
  group_by(trt) |> mutate(time_next = lead(time, default = max(time))) |> ungroup()

lr <- survdiff(Surv(follow, death) ~ trt, data = data_wide)
lr_p <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)

p_km <- ggplot(km_df, aes(time, estimate, colour = trt, fill = trt)) +
  geom_step(linewidth = 0.9) +
  geom_rect(aes(xmin = time, xmax = time_next,
                ymin = conf.low, ymax = conf.high), alpha = 0.10, colour = NA) +
  geom_point(data = filter(km_df, n.censor > 0), shape = 3, size = 1.6) +
  scale_colour_manual(values = trt_cols) + scale_fill_manual(values = trt_cols) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  annotate("text", x = 5.2, y = 0.95, label = sprintf("Log-rank p = %.3f", lr_p)) +
  labs(x = "Years since baseline", y = "Survival probability",
       colour = "Treatment", fill = "Treatment",
       title = "Kaplan–Meier survival by treatment (95% CI; + = censored)") +
  theme(legend.position = "top")
save_fig(p_km, "05_kaplan_meier.png", h = 4.6)


## ----km-table-----------------------------------------------------------------
summary(km_fit, times = c(1, 2, 3, 5))$table
summary(km_fit, times = c(1, 2, 3, 5)) |>
  (\(s) tibble(Treatment = sub("trt=", "", s$strata), Year = s$time,
               `At risk` = s$n.risk, Survival = s$surv,
               `95% CI` = sprintf("%.2f–%.2f", s$lower, s$upper)))() |>
  kable(digits = 2)
lr


## ----cox1---------------------------------------------------------------------
cox1 <- coxph(Surv(follow, death) ~ treatment + age + sex, data = data_wide)
hr_table <- function(m) {
  tidy(m, exponentiate = TRUE, conf.int = TRUE) |>
    transmute(Term = term, HR = estimate, `95% CI low` = conf.low,
              `95% CI high` = conf.high, p = p.value)
}
kable(hr_table(cox1), digits = 3)


## ----ph-test------------------------------------------------------------------
zph1 <- cox.zph(cox1)
zph1


## ----fig-schoenfeld, fig.height = 3.6, fig.width = 10-------------------------
zph_df <- map_dfr(colnames(zph1$y), \(v)
  tibble(var = v, time = zph1$x, resid = zph1$y[, v]))
coef_df <- tibble(var = names(coef(cox1)), beta = coef(cox1))
p_sch <- ggplot(zph_df, aes(time, resid)) +
  geom_hline(data = coef_df, aes(yintercept = beta), colour = "grey40", linetype = "dashed") +
  geom_point(alpha = 0.35, size = 0.9) +
  geom_smooth(method = "loess", se = TRUE, colour = "#D95F0E", fill = "#D95F0E", alpha = 0.15) +
  facet_wrap(~ var, scales = "free_y") +
  labs(title = "Scaled Schoenfeld residuals (KM-transformed time)",
       x = "Transformed time", y = "β(t) residual")
save_fig(p_sch, "06_schoenfeld.png", w = 10, h = 3.6)


## ----cox2---------------------------------------------------------------------
cox2 <- coxph(Surv(follow, death) ~ treatment + age + sex + bio_0, data = data_wide)
kable(hr_table(cox2), digits = 3)

# HR per 10-unit higher baseline biomarker (a more readable scale)
round(exp(10 * c(coef(cox2)["bio_0"], confint(cox2)["bio_0", ])), 2)

AIC(cox1, cox2)
anova(cox1, cox2)   # likelihood-ratio test for adding bio_0


## ----compare-trt--------------------------------------------------------------
bind_rows(
  `Treatment + age + sex`         = hr_table(cox1) |> filter(Term == "treatment"),
  `Treatment + age + sex + bio_0` = hr_table(cox2) |> filter(Term == "treatment"),
  .id = "Model"
) |> select(-Term) |> kable(digits = 3)


## ----session------------------------------------------------------------------
sessionInfo()

