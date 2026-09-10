library(tidyverse)
library(here)

load(here("output", "01", "riclpm_multivariate_results.Rdata"))


# Filter to only M1 (model 0) and M2 (model 5, the full model with all 5 covariates)
plot_dat <- all_params %>%
  filter(model %in% c(0, 5)) %>%
  mutate(
    model_label = ifelse(model == 0, "M1", "M2"),
    state_label = recode(state,
                         state_e = "Extraversion",
                         state_n = "Neuroticism",
                         state_o = "Openness",
                         state_a = "Agreeableness",
                         state_c = "Conscientiousness"
    )
  )

ggplot(plot_dat, aes(x = state_label, y = ri_var_est, 
                     fill = model_label, group = model_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), 
           width = 0.5) +
  geom_errorbar(aes(ymin = ri_var_ci_lower, ymax = ri_var_ci_upper),
                position = position_dodge(width = 0.6), width = 0.2) +
  scale_fill_manual(values = c("M1" = "#2C7BB6", "M2" = "#D7191C")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Random Intercept Variance",
    fill = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )


# Reshape to long format for faceting
plot_dat_long <- plot_dat %>%
  select(state_label, model_label, 
         ri_var_est, ri_var_ci_lower, ri_var_ci_upper,
         ri_trait_est, ri_trait_ci_lower, ri_trait_ci_upper) %>%
  pivot_longer(
    cols = -c(state_label, model_label),
    names_to = c("param", ".value"),
    names_pattern = "(ri_var|ri_trait)_(.*)"
  ) %>%
  mutate(param_label = ifelse(param == "ri_var", 
                              "RI Variance", 
                              "RI-Trait Correlation"))



ggplot(plot_dat_long, aes(x = model_label, y = est, 
                          fill = model_label)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2) +
  scale_fill_manual(values = c("M1" = "#2C7BB6", "M2" = "#D7191C")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = function(x) {
      if (max(x) > 1) seq(0, 1.2, by = 0.2)
      else seq(0, 0.6, by = 0.1)
    },
    limits = function(x) {
      if (max(x) > 1) c(0, 1.3)
      else c(0, 0.6)
    }
  ) +
  facet_grid(param_label ~ state_label, scales = "free_y",
             switch = "y") +
  labs(
    x = NULL,
    y = NULL,
    fill = "Model"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.3, "lines"),
    panel.spacing.y = unit(1.5, "lines")
  )

ggsave("results_figure_S1.png", 
       width = 10, height = 6, dpi = 300)
