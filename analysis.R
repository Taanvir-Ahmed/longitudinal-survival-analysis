## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 7, fig.height = 4.5,
                      fig.path = "figures/")


## ----packages-----------------------------------------------------------------
library(tidyverse)  # reading, reshaping and plotting data
library(lme4)       # linear mixed model
library(lmerTest)   # p-values for the mixed model
library(survival)   # Kaplan-Meier, log-rank test and Cox model


## ----read---------------------------------------------------------------------
data_wide <- read.table("data/extra_data.txt", header = TRUE)
head(data_wide)


## ----long---------------------------------------------------------------------
data_long <- data_wide |>
  pivot_longer(cols = starts_with("bio_"),
               names_to = "visit", values_to = "biomarker") |>
  mutate(month = as.numeric(sub("bio_", "", visit)),
         group = factor(treatment, levels = c(0, 1), labels = c("Standard", "New")))

head(data_long)
nrow(data_long)   # 240 patients x 4 visits = 960 rows


## ----trajectories-------------------------------------------------------------
group_means <- data_long |>
  group_by(group, month) |>
  summarise(biomarker = mean(biomarker))

ggplot(data_long, aes(x = month, y = biomarker, colour = group)) +
  geom_line(aes(group = id), alpha = 0.15) +
  geom_line(data = group_means, linewidth = 1.5) +
  geom_point(data = group_means, size = 3) +
  scale_x_continuous(breaks = c(0, 3, 6, 12)) +
  scale_colour_manual(values = c("Standard" = "#2C7FB8", "New" = "#D95F0E")) +
  labs(x = "Month", y = "Biomarker", colour = "Treatment",
       title = "Biomarker over 12 months") +
  theme_minimal()


## ----means-table--------------------------------------------------------------
group_means |>
  pivot_wider(names_from = month, values_from = biomarker, names_prefix = "Month ") |>
  knitr::kable(digits = 1)


## ----lmm----------------------------------------------------------------------
model <- lmer(biomarker ~ month * treatment + (1 | id), data = data_long)
summary(model)


## ----lmm-ci-------------------------------------------------------------------
confint(model, parm = "beta_", method = "Wald")


## ----km-----------------------------------------------------------------------
km <- survfit(Surv(follow, death) ~ treatment, data = data_wide)
km


## ----kaplan-meier-------------------------------------------------------------
plot(km, col = c("#2C7FB8", "#D95F0E"), lwd = 2, mark.time = TRUE,
     xlab = "Years", ylab = "Probability of survival",
     main = "Kaplan–Meier survival curves")
legend("bottomleft", legend = c("Standard", "New"),
       col = c("#2C7FB8", "#D95F0E"), lwd = 2)


## ----km-3years----------------------------------------------------------------
summary(km, times = 3)


## ----logrank------------------------------------------------------------------
survdiff(Surv(follow, death) ~ treatment, data = data_wide)


## ----cox----------------------------------------------------------------------
cox <- coxph(Surv(follow, death) ~ treatment + age + sex, data = data_wide)
summary(cox)


## ----ph-check-----------------------------------------------------------------
cox.zph(cox)

