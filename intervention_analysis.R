# Author: Auja Bywater
# Last updated: 07/08/2026

# This analysis is for the intervention testing in my hydroponic system. 
# We tested bacteriophages and bacteriocins as an intervention to reduce persistence of foodborne pathogens

# Load packages
pacman::p_load(summarytools, dplyr, ggplot2, lmerTest, readxl, tidyr, emmeans)
library(ggbreak)
library(nlme)

# Load in the data
data <- read_excel("all_intervention.xlsx", 3)

# Subset my treatment
cin <- data %>% 
  filter(intervention =="cin")
phage <- data %>% 
  filter(intervention == "phage")

## Graph phage data 

# Facetted by rep
phage_averaged <- phage %>%
  group_by(day, group, rep, bacteria) %>%
  summarize(
    mean_log_count = mean(log_cfu_g, na.rm = TRUE), 
    .groups = "drop"
  )

ggplot(phage_averaged, aes(x = factor(day), y = mean_log_count, color = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_grid(bacteria ~ rep, labeller = label_both) + 
  labs(
    title = "Phage Intervention: Bacterial Persistence by Independent Replicate",
    x = "Day",
    y = "Averaged Log Count (CFU/mL)",
    color = "Treatment Group"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Independent reps combined 
phage_summary_all <- phage %>%
  group_by(day, group, bacteria) %>%
  summarize(
    mean_log_count = mean(log_cfu_g),
    se = sd(log_cfu_g) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

ggplot(phage_summary_all,
       aes(x = factor(day),
           y = mean_log_count,
           color = group,
           group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_count - se,
                    ymax = mean_log_count + se),
                width = 0.2,
                linewidth = 0.6) +
  facet_wrap(~bacteria) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, by = 1)
  ) +
  labs(
    title = "Phage Intervention: Mean Bacterial Persistence",
    x = "Day",
    y = "Log Count (CFU/mL)",
    color = "Treatment Group"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )


# Make day equally spaced
phage$day <- factor(phage$day, levels = c(0, 1, 2, 3, 21))

# Fit the model with temp and humidity
phagein <- lme(
  log_cfu_g ~ day * group + temp + humidity,
  random = ~1 | rep/system,
  weights = varIdent(form = ~1 | day),
  data = phage,
  control = lmeControl(
    maxIter = 200,
    msMaxIter = 200,
    niterEM = 100,
    opt = "optim"
  )
)

# Check if temp and humidity are needed
full_ml_p <- update(phagein, method = "ML")

# Check without covariates
reduced_ml_p <- update(full_ml_p, . ~ . - temp - humidity)
no_temp_p <- update(full_ml_p, . ~ . - temp)
no_hum <- update(full_ml_p, . ~ . - humidity)


anova(reduced_ml_p, full_ml_p)
anova(no_temp_p, full_ml_p)
anova(no_hum, full_ml_p)

# Temp and humidity are worth keeping
plot(phagein)

# Interpret model
anova(phagein)
joint_tests(phagein)

phage_emm <- emmeans(
  phagein,
  pairwise ~ group | day
)

phage_emm

# Estimated marginal means with 95% CI
phage_summary <- emmeans(phagein, ~ group | day)

phage_summary <- as.data.frame(confint(phage_summary))

# Labels for groups
group_labels <- c(
  bacteriophage = "Bacteriophage",
  ro = "RO Control"
)

# Plot with Lmer CI
phage_plot <- ggplot(
  phage_summary,
  aes(
    x = as.numeric(as.character(day)),
    y = emmean,
    color = group
  )
) +
  
  # 95% confidence interval shading
  geom_ribbon(
    aes(
      ymin = lower.CL,
      ymax = upper.CL,
      fill = group,
      group = group
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  # Estimated mean line
  geom_line(
    aes(group = group),
    linewidth = 1.2
  ) +
  
  # Sampling day points
  geom_point(
    size = 3
  ) +
  
  # Color palette
  scale_color_manual(
    values = c(
      bacteriophage = "#D55E00",
      ro = "#009E73"
    ),
    labels = group_labels
  ) +
  
  scale_fill_manual(
    values = c(
      bacteriophage = "#D55E00",
      ro = "#009E73"
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = c(0,1,2,3,21),
    labels = c(0,1,2,3,21)
  ) +
  
  scale_y_continuous(
    limits = c(0,8)
  ) +
  
  labs(
    x = "Day",
    y = expression(bold(log[10] ~ CFU/g)),
    color = "Treatment"
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

phage_plot

ggsave(
  "bacteriophage.png",
  plot = phage_plot,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)











  


## Graph cin data
cin_averaged <- cin %>%
  group_by(day, group, rep, bacteria) %>%
  summarize(
    mean_log_count = mean(log_cfu_g, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(cin_averaged,
       aes(x = factor(day),
           y = mean_log_count,
           color = bacteria,
           linetype = group,
           group = interaction(bacteria, group))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~rep) +
  labs(
    title = "bacteriocin treatment",
    x = "Day",
    y = "Averaged Log Count (CFU/mL)",
    color = "Bacteria",
    linetype = "Treatment Group"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# bacteria facet
ggplot(
  cin_averaged,
  aes(
    x = factor(day),
    y = mean_log_count,
    color = bacteria,
    linetype = group,
    group = interaction(bacteria, group)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_grid(
    bacteria ~ rep,
    labeller = labeller(bacteria = bacteria_labels)
  ) +
  scale_color_manual(
    values = c(
      e_coli = "#0072B2",
      listeria = "#E69F00",
      salmonella = "#882255"
    ),
    labels = bacteria_labels
  ) +
  labs(
    title = "Bacteriocin Treatment",
    x = "Day",
    y = expression(log[10] ~ CFU/mL),
    linetype = "Treatment Group"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  guides(color = "none")

# Averaged independent reps
cin_summary_all <- cin %>%
  group_by(day, group, bacteria) %>%
  summarize(
    mean_log_count = mean(log_cfu_g, na.rm = TRUE),
    se = sd(log_cfu_g, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(cin_summary_all,
       aes(x = factor(day),
           y = mean_log_count,
           color = bacteria,
           linetype = group,
           group = interaction(bacteria, group))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_log_count - se,
                    ymax = mean_log_count + se),
                width = 0.2,
                linewidth = 0.6) +
  labs(
    title = "Cinnamycin Intervention: Bacterial Persistence",
    x = "Day",
    y = "Mean Log Count (CFU/g)",
    color = "Bacteria",
    linetype = "Treatment Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

# Facetted
ggplot(
  cin_summary_all,
  aes(
    x = factor(day),
    y = mean_log_count,
    color = bacteria,
    linetype = group,
    group = interaction(bacteria, group)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(
      ymin = mean_log_count - se,
      ymax = mean_log_count + se
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ bacteria,
    labeller = labeller(bacteria = bacteria_labels)
  ) +
  labs(
    x = "Day",
    y = expression(log[10] ~ CFU/mL),
    linetype = "Treatment Group"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  ) +
  guides(color = "none")





# Run stats

# Make day equally spaced
cin$day <- factor(cin$day, levels = c(0, 1, 2, 3, 21))

# Fit the model
cinin <- lme(
  log_cfu_g ~ day * group * bacteria + temp + humidity,
  random = ~1 | rep/system,
  weights = varIdent(form = ~1 | day * bacteria),
  data = cin,
  control = lmeControl(
    maxIter = 200,
    msMaxIter = 200,
    niterEM = 100,
    opt = "optim"
  )
)

# Check if temp and humidity are needed
full_ml_b <- update(cinin, method = "ML")

# Check without covariates
reduced_ml_b <- update(full_ml_b, . ~ . - temp - humidity)
no_temp_b <- update(full_ml_b, . ~ . - temp)
no_hum_b <- update(full_ml_b, . ~ . - humidity)


anova(reduced_ml_b, full_ml_b)
anova(no_temp_b, full_ml_b)
anova(no_hum_b, full_ml_b)

# Humidity is necessary, temp isn't necessarily
# For consistency with the bacteriophage treatment, I'm keeping them both

# Interpret model
anova(cinin)

cin_emm <- emmeans(
  cinin,
  pairwise ~ group | bacteria *day
)

cin_emm


# Estimated marginal means with 95% CI
cin_summary <- emmeans(cinin, ~ group | bacteria * day)

cin_summary <- as.data.frame(confint(cin_summary))

# Italicized bacteria labels
bacteria_labels <- c(
  e_coli = "bolditalic('E. coli')",
  listeria = "bolditalic('Listeria monocytogenes')",
  salmonella = "bolditalic('Salmonella')"
)

group_labels <- c(
  bacteriocin = "Bacteriocin",
  ro_only = "RO Control"
)


# Make day numerical just for graph
cin_summary$day <- as.numeric(as.character(cin_summary$day))

# Graph
cin_plot <- ggplot(
  cin_summary,
  aes(
    x = day,
    y = emmean,
    color = group
  )
) +
  
  # 95% CI shading
  geom_ribbon(
    aes(
      ymin = lower.CL,
      ymax = upper.CL,
      fill = group,
      group = group
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  # Estimated mean lines
  geom_line(
    aes(group = group),
    linewidth = 1.2
  ) +
  
  # Sampling points
  geom_point(
    size = 3
  ) +
  
  # Treatment colors
  scale_color_manual(
    values = c(
      bacteriocin = "#CC79A7",
      ro_only = "#009E73"
    ),
    labels = group_labels
  ) +
  
  scale_fill_manual(
    values = c(
      bacteriocin = "#CC79A7",
      ro_only = "#009E73"
    ),
    guide = "none"
  ) +
  
  # Days
  scale_x_continuous(
    breaks = c(0,1,2,3,21),
    labels = c(0,1,2,3,21)
  ) +
  
  scale_y_continuous(
    limits = c(0,8)
  ) +
  
  # Facet by bacteria
  facet_wrap(
    ~ bacteria,
    labeller = as_labeller(bacteria_labels, label_parsed)
  ) +
  
  labs(
    x = "Day",
    y = expression(bold(log[10] ~ CFU/g)),
    color = "Treatment"
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

cin_plot

ggsave(
  "bacteriocin.png",
  plot = cin_plot,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)


