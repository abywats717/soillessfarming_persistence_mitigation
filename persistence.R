#Date: March 22, 2026
#Author: Auja Bywater
#Project: NE SARE Persistence
# This data looks at how E. coli, Salmonella, and Listeria persist in a deep water culture
# and Kratky hydroponic system from seeding to harvest

# Questions:
#     How well do each of the species persist over time?
#     How does persistence vary between species?
#     Does hydroponic system type influence persistence?

# Load packages
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(nlme)
library(emmeans)
library(ggbreak)

# Load in the data
# nursery stage (seedling to transplant ~ 3 weeks)
nursery <- read_excel("all_persistence.xlsx", 3)

# system stage (transplant to harvest ~ 4 weeks)
system <- read_excel("all_persistence.xlsx", 4)

# harvest (samples taken at harvest)
harvest <- read_excel("all_persistence.xlsx", 5)

# environment (environmental data over the experiment)
environ <- read_excel("all_persistence.xlsx", 6)





### NURSERY DATA ###
# Convert days to factor
nursery_filtered <- nursery %>%
  mutate(day = as.factor(day))  

# Give independent replicates a unique ID
nursery_filtered$rep_id <- paste(nursery_filtered$design, nursery_filtered$replicate, sep = "_")

# Exploratory Stats
nursery_reps <- nursery %>%
  group_by(day, bacteria, replicate) %>%
  summarise(
    replicate_mean = mean(log_count, na.rm = TRUE),
    .groups = "drop"
  )

nursery_avg <- nursery_reps %>%
  group_by(day, bacteria) %>%
  summarise(
    mean_log_count = mean(replicate_mean),
    se = sd(replicate_mean) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

nursery_avg

# Fit the model
nursery_model <- lme(
  log_count ~ day  * bacteria,
  random = ~1 | rep_id/bio_rep,
  weights = varIdent(form = ~1 | day * bacteria),
  data = nursery_filtered,
  control = lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 100, opt = "optim")
)

# Make sure the model is a good fit
## Variance by day
nursery_model_homo <- update(nursery_model, weights = NULL)
anova(nursery_model, nursery_model_homo)  # LRT comparing variance structures

## Variance by bacteria only
nursery_model_bact <- update(nursery_model, weights = varIdent(form = ~1 | bacteria))

## Variance by day and bacteria combined (interaction of both grouping factors)
nursery_model_combo <- update(nursery_model, weights = varIdent(form = ~1 | day * bacteria))

## Compare all four structures by AIC/BIC
AIC(nursery_model, nursery_model_bact, nursery_model_combo, nursery_model_homo)
BIC(nursery_model, nursery_model_bact, nursery_model_combo, nursery_model_homo)

## Does day*bacteria meaningfully improve over bacteria-only?
anova(nursery_model_combo, nursery_model_bact)

## Does day*bacteria meaningfully improve over day-only?
anova(nursery_model_combo, nursery_model)

## Sample size per day x bacteria cell
table(nursery_filtered$day, nursery_filtered$bacteria)

## Check for convergence issues or wild variance estimates
summary(nursery_model_combo)$modelStruct$varStruct

# Check residuals normality 
plot(nursery_model)  

# Check residual normality
qqnorm(resid(nursery_model, type = "normalized"))
qqline(resid(nursery_model, type = "normalized"))
qqnorm(ranef(nursery_model)$rep_id[,1]); qqline(ranef(nursery_model)$rep_id[,1])

nursery_filtered$nursery_model <- resid(nursery_model, type = "normalized")




## Look at the model
anova(nursery_model)

# Check individual bacteria (sanity check only for overall trends)
emm_bacteria <- emmeans(
  nursery_model,
  ~ bacteria
)

pairs(emm_bacteria, adjust = "tukey")

# Compare bacteria at each day
emm_day <- emmeans(
  nursery_model,
  ~ bacteria | day
)

pairs(emm_day, adjust = "tukey")

# Track individual bacteria over time
emm_time <- emmeans(
  nursery_model,
  ~ day | bacteria
)

## Check consecutive days
consec_contrasts <- contrast(emm_time, method = "consec", adjust = "tukey")
consec_contrasts

## Check day 0 vs 21
day0_vs_day21 <- contrast(
  emm_time,
  method = "trt.vs.ctrl",
  ref = "day0"
)

day0_vs_day21_only <- as.data.frame(day0_vs_day21) %>%
  filter(contrast == "day21 - day0")

day0_vs_day21_only




# Prep for graphing
nursery_summary <- emmeans(
  nursery_model,
  ~ day * bacteria
)

nursery_summary <- as.data.frame(
  confint(nursery_summary)
)

# Custom labels for Bacteria (italicized) and Design
bacteria_labels <- c(
  e_coli = expression(bolditalic("E. coli")),
  listeria = expression(bolditalic("L. monocytogenes")),
  salmonella = expression(bolditalic("S. enterica"))
)

design_labels <- c(
  dwc = "Deep Water Culture",
  kr = "Kratky"
)

# Plot
# Plot model-estimated means with 95% confidence intervals
nursery_plot <- ggplot(
  nursery_summary,
  aes(
    x = as.numeric(as.character(day)),
    y = emmean,
    color = bacteria
  )
) +
  # 95% confidence interval shading
  geom_ribbon(
    aes(
      ymin = lower.CL,
      ymax = upper.CL,
      fill = bacteria,
      group = bacteria
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  # Model-estimated mean line
  geom_line(
    aes(group = bacteria),
    linewidth = 1.2
  ) +
  
  # Sampling day points
  geom_point(
    size = 3
  ) +
  
  # Colorblind-friendly palette
  scale_color_manual(
    values = c(
      e_coli = "#0072B2",
      listeria = "#E69F00",
      salmonella = "#882255"
    ),
    labels = bacteria_labels
  ) +
    scale_fill_manual(
    values = c(
      e_coli = "#0072B2",
      listeria = "#E69F00",
      salmonella = "#882255"
    ),
    guide = "none"
  ) +
  
  scale_y_continuous(
    limits = c(0,8)
  ) +
    labs(
    x = "Day",
    y = expression(bold(log[10] ~ CFU/g)),
    color = "Bacteria"
  ) +
    theme_minimal(
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.6
    ),
    axis.ticks.length = unit(0.15, "cm"),
    strip.text = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 13,
      color = "black"
    ),
    legend.text = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),
    legend.position = "right"
  )

nursery_plot

ggsave(
  "nursery_persistence_model_CI.png",
  plot = nursery_plot,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)











### SYSTEM DATA ###

# Reshape data to long format for plotting cup and corner together
system_long <- system %>%
  pivot_longer(
    cols = c(cup, corner),
    names_to = "location",
    values_to = "log_count"
  ) %>%
  filter(day != 21) %>%
  mutate(day = as.factor(day))
# day as factor
system_long <- system_long %>%
  mutate(day = as.factor(day))

# Give independent replicates a unique ID
system_long$rep_id <- paste(system_long$design, system_long$Replicate, sep = "_")

# Exploratory stats
system_reps <- system %>%
  group_by(design, day, bacteria, Replicate) %>%
  summarise(
    cup_mean = mean(cup, na.rm = TRUE),
    corner_mean = mean(corner, na.rm = TRUE),
    .groups = "drop"
  )

system_avg <- system_reps %>%
  group_by(design, day, bacteria) %>%
  summarise(
    mean_cup = mean(cup_mean),
    se_cup = sd(cup_mean) / sqrt(n()),
    mean_corner = mean(corner_mean),
    se_corner = sd(corner_mean) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

system_avg

# Fit the model with ph, ec, and ns_temp
system_model <- lme(
  log_count ~ day  * bacteria * design * location + ph + ec + ns_temp,
  random = ~1 | rep_id/bio_rep,
  weights = varIdent(form = ~1 | day * bacteria),
  data = system_long,
  control = lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 100, opt = "optim")
)

# Check residuals
plot(system_model)

# Check if ph, ec, and ns_temp are needed
full_ml <- update(system_model, method = "ML")

# Check without covariates
reduced_ml <- update(full_ml, . ~ . - ph - ec - ns_temp)
no_ph <- update(full_ml, . ~ . - ph)
no_ec <- update(full_ml, . ~ . - ec)
no_temp <- update(full_ml, . ~ . - ns_temp)

anova(reduced_ml, full_ml)
anova(no_ph, full_ml)
anova(no_ec, full_ml)
anova(no_temp, full_ml)

# Covariates don't improve fit enough to justify extra parameter

# Refit the model without covariates 
system_model <- lme(
  log_count ~ day  * bacteria * design * location,
  random = ~1 | rep_id/bio_rep,
  weights = varIdent(form = ~1 | day * bacteria),
  data = system_long,
  control = lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 100, opt = "optim")
)

# Check convergence 
system_model$apVar

# Check variance
VarCorr(system_model)

# bio_rep variance is very low - check if it's needed
sm_ml <- update(system_model, method = "ML")
sm_norep <- update(sm_ml, random = ~1 | rep_id)
anova(sm_norep, sm_ml)

# Removing bio_rep to simplify model as it isn't needed
system_model <- lme(
  log_count ~ day * bacteria * design * location,
  random = ~1 | rep_id,
  weights = varIdent(form = ~1 | day * bacteria),
  data = system_long,
  control = lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 100, opt = "optim")
)

# Confirm how the model is handling variance is appropriate
m_day        <- system_model  # current: weights = varIdent(form = ~1 | day)
m_bacteria   <- update(system_model, weights = varIdent(form = ~1 | bacteria))
m_design     <- update(system_model, weights = varIdent(form = ~1 | design))
m_location   <- update(system_model, weights = varIdent(form = ~1 | location))
m_day_bact   <- update(system_model, weights = varIdent(form = ~1 | day * bacteria))

AIC(m_day, m_bacteria, m_design, m_location, m_day_bact)

anova(m_day, m_day_bact)
anova(m_bacteria, m_day_bact)
anova(m_bacteria, m_day_bact)
table(system_long$day, system_long$bacteria)


# Final model check
system_model_sum$apVar

system_model$apVar
VarCorr(system_model)
plot(system_model)
qqnorm(resid(system_model, type = "normalized")); qqline(resid(system_model, type = "normalized"))

# Model looks good - moving on to results

joint_tests(system_model)

# Compare bacteria 
emmeans(system_model, pairwise ~ bacteria)

# Compare systems
emmeans(system_model, pairwise ~ design)

# Compare location
emmeans(system_model, pairwise ~ location)

# Compare bacteria within each combination of day/system/location
emm <- emmeans(system_model, pairwise ~ bacteria | day * design * location)
emm_means <- summary(emm$emmeans)
emm_means_df <- as.data.frame(emm_means)

emm_contrasts <- summary(emm$contrasts, adjust = "tukey")
emm_contrasts_df <- as.data.frame(emm_contrasts)

ggplot(emm_means_df, aes(x = bacteria, y = emmean, group = design, color = design)) +
  geom_point(position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.3)) +
  facet_grid(location ~ day) +
  labs(y = "Estimated marginal mean (log_count)", x = "Bacteria") +
  theme_minimal()

emm_contrasts_df[emm_contrasts_df$p.value < 0.05,]

# Compare design within each combination
emm_design <- emmeans(system_model, pairwise ~ design | bacteria * location * day, adjust = "tukey")
design_df <- as.data.frame(emm_design$contrasts)
design_df[design_df$p.value < 0.05, ]

# Compare location within each combination
emm_loc <- emmeans(system_model, pairwise ~ location | bacteria * design * day, adjust = "tukey")
loc_df <- as.data.frame(emm_loc$contrasts)
loc_df[loc_df$p.value < 0.05, ]

# prep for graph
system_summary <- emmeans(
  system_model,
  ~ day * design * bacteria * location
)

system_summary <- as.data.frame(
  confint(system_summary)
)

# Sort the data
system_summary <- system_summary %>%
  arrange(design, bacteria, location, day)

# Plot
system_plot <- ggplot(
  system_summary,
  aes(
    x = as.numeric(as.character(day)),
    y = emmean,
    color = bacteria
  )
) +
  
  # 95% confidence interval ribbons
  geom_ribbon(
    aes(
      ymin = lower.CL,
      ymax = upper.CL,
      fill = bacteria,
      group = interaction(design, bacteria, location)
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  # Model-estimated mean lines
  geom_line(
    aes(
      group = interaction(design, bacteria, location)
    ),
    linewidth = 1.2
  ) +
  
  # Sampling day points
  geom_point(
    size = 3
  ) +
  
  # Facet: hydroponic system rows, sampling location columns
  facet_grid(
    design ~ location,
    labeller = labeller(
      design = design_labels,
      location = c(
        cup = "Cup",
        corner = "Corner"
      )
    )
  ) +
  
  # Bacteria colors for lines
  scale_color_manual(
    values = c(
      e_coli = "#0072B2",
      listeria = "#E69F00",
      salmonella = "#882255"
    ),
    labels = bacteria_labels
  ) +
  
  # Bacteria colors for confidence interval ribbons
  scale_fill_manual(
    values = c(
      e_coli = "#0072B2",
      listeria = "#E69F00",
      salmonella = "#882255"
    ),
    labels = bacteria_labels,
    guide = "none"
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 6, 2)
  ) +
  
  coord_cartesian(
    ylim = c(0, 6)
  ) +
  
  labs(
    x = "Day",
    y = expression(bold(log[10] ~ CFU/mL)),
    color = "Bacteria"
  ) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.6
    ),
    axis.ticks.length = unit(
      0.15,
      "cm"
    ),
    strip.text = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    legend.title = element_text(
      face = "bold",
      size = 13,
      color = "black"
    ),
    legend.text = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),
    legend.position = "right"
  )

system_plot

ggsave(
  "system_location_persistence_model_CI.png",
  plot = system_plot,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)

# The replicates don't overlap in time - checking consistency 
system_long$date_parsed <- as.Date(system_long$date, format = "%m.%d.%y")
system_long$month <- format(system_long$date_parsed, "%Y-%m")
table(system_long$design, system_long$month) 

system_model_season <- lme(
  log_count ~ day * bacteria * design * location + month,
  random = ~1 | rep_id,
  weights = varIdent(form = ~1 | day * bacteria),
  data = system_long,
  control = lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 100, opt = "optim")
)

joint_tests(system_model_season)

system_model_season$apVar






### Harvest Data ###

# convert to long
harvest_long <- harvest %>%
  pivot_longer(
    cols = c(salmonella, e_coli, listeria),
    names_to = "bacteria",
    values_to = "count"
  )

# bring enrichment results along, matched to the correct pathogen, so
# prevalence stats below can account for enrichment-only positives
enrich_long <- harvest %>%
  select(design, replicate, date, location, bio_rep, s_enrichment, e_enrichment, l_enrichment) %>%
  rename(
    salmonella = s_enrichment,
    e_coli     = e_enrichment,
    listeria   = l_enrichment
  ) %>%
  pivot_longer(
    cols = c(salmonella, e_coli, listeria),
    names_to = "bacteria",
    values_to = "enrichment"
  )

harvest_long <- harvest_long %>%
  left_join(enrich_long, by = c("design", "replicate", "date", "location", "bio_rep", "bacteria"))

# summarize data

harvest_summary <- harvest_long %>%
  group_by(location, design, bacteria) %>%
  summarise(
    n                    = n(),
    mean_log             = mean(count, na.rm = TRUE),
    se_log               = sd(count, na.rm = TRUE) / sqrt(n()),
    n_direct_pos         = sum(count > 0),
    pct_direct_pos       = round(100 * mean(count > 0), 1),
    n_enrich_pos         = sum(count > 0 | enrichment == "+"),
    pct_enrich_pos       = round(100 * mean(count > 0 | enrichment == "+"), 1),
    mean_of_positives    = if (any(count > 0)) round(mean(count[count > 0]), 2) else NA,
    median_of_positives  = if (any(count > 0)) round(median(count[count > 0]), 2) else NA,
    min_of_positives     = if (any(count > 0)) round(min(count[count > 0]), 2) else NA,
    max_of_positives     = if (any(count > 0)) round(max(count[count > 0]), 2) else NA,
    .groups = "drop"
  )

print(harvest_summary, n = Inf)



# Fix facet titles
harvest_long$bacteria <- factor(
  harvest_long$bacteria,
  levels = c("e_coli", "listeria", "salmonella"),
  labels = c("E. coli", "L. monocytogenes", "S. enterica")
)

# Fix facet titles
harvest_long$design <- factor(
  harvest_long$design,
  levels = c("dwc", "kr"),
  labels = c("Deep Water Culture", "Kratky")
)

# Remove leaf data for ploting
harvest_long <- harvest_long %>%
  filter(location != "leaves")

# Plot
harvest_plot <- ggplot(
  harvest_long,
  aes(x = bacteria, y = count, fill = bacteria)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    aes(color = bacteria),
    width = 0,
    height = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  scale_y_continuous(
    breaks = seq(0, 5, by = 1)
  ) +
  coord_cartesian(ylim = c(0, 6)) +
  facet_grid(location ~ design) +
  scale_fill_manual(
    values = c(
      "E. coli" = "#0072B2",
      "L. monocytogenes" = "#E69F00",
      "S. enterica" = "#882255"
    ),
    labels = bacteria_labels
  ) +
  scale_color_manual(
    values = c(
      "E. coli" = "#0072B2",
      "L. monocytogenes" = "#E69F00",
      "S. enterica" = "#882255"
    ),
    labels = bacteria_labels
  ) +
  labs(
    y = expression(bold(log[10] ~ CFU/mL)),
    x = NULL,
    fill = "Bacteria",
    color = "Bacteria"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(
      face = "bold",
      size = 18,
      color = "black"
    ),
    strip.text.y = element_text(
      face = "bold",
      size = 16,
      color = "black",
      angle = 0
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 18,
      color = "black"
    ),
    axis.text.x = element_text(
      face = "bold.italic",
      size = 14,
      color = "black"
    ),
    axis.text.y = element_text(
      face = "bold.italic",
      size = 14,
      color = "black"
    ),
    panel.spacing = unit(1, "lines")
  )
harvest_plot
ggsave("harvest_plot.png", plot = harvest_plot, width = 16, height = 8, dpi = 600, bg = "white")


# Descriptive analysis

# Replicate structure check 
harvest %>%
  distinct(design, replicate, date) %>%
  arrange(design, date)

harvest %>%
  group_by(design, location) %>%
  summarize(n_replicates = n_distinct(replicate), .groups = "drop")

# Sanity check against the text's reported ranges 
harvest_summary %>% filter(location %in% c("rockwool", "roots"))
harvest_summary %>% filter(location == "leaves")

# Verify "DWC leaves: Listeria undetectable by direct plating
harvest_long %>%
  filter(design == "Deep Water Culture", location == "leaves", bacteria == "L. monocytogenes") %>%
  summarize(
    n = n(),
    n_direct_pos = sum(count > 0),
    n_enrich_pos = sum(count > 0 | enrichment == "+"),
    pct_direct   = round(100 * mean(count > 0), 1),
    pct_enrich   = round(100 * mean(count > 0 | enrichment == "+"), 1)
  )

# Design/date (season) confound check 
harvest %>%
  distinct(design, replicate, date) %>%
  mutate(date_parsed = as.Date(date, format = "%m.%d.%y")) %>%
  arrange(date_parsed)