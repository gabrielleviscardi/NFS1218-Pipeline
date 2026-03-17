# =============================================================================
# PRECISION NUTRITION CLUSTERING PIPELINE
#
# RESEARCH QUESTION:
#   Is maternal UPF score associated with distinct zBMI trajectory clusters
#   from birth to 24 months? 
#   Do children of mothers with higher UPF scores
#   belong to higher-risk zBMI growth clusters, and does maternal dietary
#   quality predict cluster membership at 24 months?
#
#
# GROUP MEMBERS: Gabrielle Viscardi, Brighid McKay, Diana Ghidanac
# DATE: April 7, 2026

#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                            Line
#TOC> -----------------------------------------------------------------
#TOC>   1        Packages                                          35
#TOC>   2        Data loading                                      50
#TOC>   3        Missingness exploration                           61    
#TOC>   4        Distribution inspection                           78
#TOC>   5        Identifying implausible values                   117
#TOC>   6        Data cleaning                                    162
#TOC>   7        Post-cleaning diagnostics                        221
#TOC>   8        Data preparation for clustering                  239
#TOC>   9        Assessing clustering tendency & optimal k        287
#TOC>   10       Clustering                                       376
#TOC>   11       Cluster validation                               486
#TOC>   12       Answering the research question                  562
#TOC>   13       Save final results                               783
#TOC> ===========================================================================
#
# =============================================================================
# SECTION 1 — PACKAGES
# =============================================================================

packages <- c("tidyverse", "naniar", "skimr", "NbClust", "factoextra",
              "cluster", "clValid", "dendextend", "gridExtra",
              "RColorBrewer", "nnet", "reshape2", "ggplot2")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# SECTION 2 — DATA LOADING 
# =============================================================================

setwd("~/Desktop")
data <- read_csv("mock_precision_growth_dataset.csv")

# Inspect structure
glimpse(data)
summary(data)

# =============================================================================
# SECTION 3 — MISSINGNESS EXPLORATION
# =============================================================================

missing_summary <- data %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(),
               names_to  = "variable",
               values_to = "percent_missing") %>%
  arrange(desc(percent_missing))

print(missing_summary)

# Visual inspection of missing data 
vis_miss(data)
gg_miss_upset(data)

# =============================================================================
# SECTION 4 — DISTRIBUTION INSPECTION
# =============================================================================

# zBMI distribution at birth, 12 months, and 24 months
if ("WHO_zBMI_birth" %in% names(data)) {
  print(ggplot(data, aes(x = WHO_zBMI_birth)) +
          geom_histogram(bins = 30) +
          theme_minimal() +
          labs(title = "Distribution of zBMI at birth")
  )
}

if ("WHO_zBMI_12m" %in% names(data)) {
  print(ggplot(data, aes(x = WHO_zBMI_12m)) +
          geom_histogram(bins = 30) +
          theme_minimal() +
          labs(title = "Distribution of zBMI at 12 months")
  )
}

if ("zBMI_24m" %in% names(data)) {
  print(ggplot(data, aes(x = zBMI_24m)) +
          geom_histogram(bins = 30) +
          theme_minimal() +
          labs(title = "Distribution of zBMI at 24 months")
  )
}

if ("Ultra_processed_score" %in% names(data)) {
  print(ggplot(data, aes(x = Ultra_processed_score)) +
          geom_histogram(bins = 25) +
          theme_minimal() +
          labs(title = "Maternal UPF Score")
  )
}

# =============================================================================
# SECTION 5 — IDENTIFYING IMPLAUSIBLE VALUES
# =============================================================================

# Implausible age values
# NOTE: Age_24m_months has values of -3, 5, and 120 and they should all be 24

if ("Age_24m_months" %in% names(data)) {
  
  implausible_age <- data %>%
    filter(Age_24m_months < 0 | Age_24m_months > 60)
  
  print(implausible_age)
}

# WHO-style plausibility cut-offs for zBMI
# WHO commonly flags z-scores < -5 or > +5 as implausible at any age
# WHO_zBMI_birth has 11 outliers < -5 and no missing values
# WHO_zBMI_12m has 1 outlier < -5 and 30 missing values
# zBMI_24m has one extreme outlier > +5 and no missing values

if ("WHO_zBMI_birth" %in% names(data)) {
  
  implausible_zBMI <- data %>%
    filter(WHO_zBMI_birth < -5 | WHO_zBMI_birth > 5)
  
  print(implausible_zBMI)
}

if ("WHO_zBMI_12m" %in% names(data)) {
  
  implausible_zBMI <- data %>%
    filter(WHO_zBMI_12m < -5 | WHO_zBMI_12m > 5)
  
  print(implausible_zBMI)
}

if ("zBMI_24m" %in% names(data)) {
  
  implausible_zBMI <- data %>%
    filter(zBMI_24m < -5 | zBMI_24m > 5)
  
  print(implausible_zBMI)
}

# =============================================================================
# SECTION 6 — DATA CLEANING
#   1. Age_24m_months: has values of -3, 5, and 120 (should all be 24)
#   2. WHO_zBMI_birth: has 11 outliers < -5
#   3. WHO_zBMI_12m: has 1 outlier < -5 and 30 missing values
#   4. zBMI_24m: one extreme outlier above +5 
# =============================================================================

clean_data <- data

# Implausible age values
if ("Age_24m_months" %in% names(data)) {
  implausible_age <- data %>%
    filter(Age_24m_months < 0 | Age_24m_months > 60)
  print(implausible_age)
}

# WHO-style plausibility cut-offs for zBMI
if ("WHO_zBMI_birth" %in% names(data)) {
  implausible_zbmi_birth <- data %>%
    filter(WHO_zBMI_birth < -5 | WHO_zBMI_birth > 5)
  print(implausible_zbmi_birth)
}

if ("WHO_zBMI_12m" %in% names(data)) {
  implausible_zbmi_12m <- data %>%
    filter(WHO_zBMI_12m < -5 | WHO_zBMI_12m > 5)
  print(implausible_zbmi_12m)
}

if ("zBMI_24m" %in% names(data)) {
  implausible_zbmi_24m <- data %>%
    filter(zBMI_24m < -5 | zBMI_24m > 5)
  print(implausible_zbmi_24m)
}

# Handle missing values

clean_data <- clean_data %>%
  filter(!is.na(zBMI_24m),
         !is.na(WHO_zBMI_12m),
         !is.na(WHO_zBMI_birth),
         !is.na(Ultra_processed_score))

# WHO_zBMI_12m has 30 missing values (10%)
# This is a significant amount of missing value that are core trajectory variables
# Dropping those values would lose too much data
# We will impute missing data using the median zBMI at 12 months

cols_to_impute <- c("WHO_zBMI_12m")
for (col in cols_to_impute) {
  if (col %in% names(clean_data)) {
    med_val <- median(clean_data[[col]], na.rm = TRUE)
    n_miss  <- sum(is.na(clean_data[[col]]))
    clean_data[[col]][is.na(clean_data[[col]])] <- med_val
    if (n_miss > 0) cat(sprintf("Imputed %d missing values in '%s' with median (%.3f)\n",
                                n_miss, col, med_val))
  }
}

# =============================================================================
# SECTION 7 — POST-CLEANING DIAGNOSTICS
# =============================================================================

# Recalculate missingness after cleaning
missing_summary_clean <- clean_data %>%
  summarise(across(everything(),
                   ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(),
               names_to = "variable",
               values_to = "percent_missing") %>%
  arrange(desc(percent_missing))

print(missing_summary_clean)

# Save cleaned dataset
write_csv(clean_data, "clean_precision_growth_dataset.csv")

# =============================================================================
# SECTION 8 — DATA PREPARATION FOR CLUSTERING
#
# We have zBMI at 3 time points (birth, 12-months, 24-months)
# =============================================================================

# Install and load required packages
packages <- c("NbClust", "factoextra", "ggplot2", "gridExtra", "cluster", 
              "RColorBrewer", "reshape2")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Select growth trajectory features

growth_vars <- c("WHO_zBMI_birth", "WHO_zBMI_12m", "zBMI_24m")

growth_matrix <- clean_data %>%
  select(all_of(growth_vars)) %>%
  as.data.frame()

rownames(growth_matrix) <- clean_data$ID

# Standardize the data 

growth_scaled <- scale(growth_matrix)
cat("\nMean of each column after scaling (should all be ~0):\n")
print(round(colMeans(growth_scaled), 3))
cat("\nSD of each column after scaling (should all be ~1):\n")
print(round(apply(growth_scaled, 2, sd), 3))

# Compute and visualize distance matrix 
dist_eucl <- dist(growth_scaled, method = "euclidean")

# Visualize: children who are close together (similar trajectories) appear
# in the same color block
fviz_dist(dist_eucl,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) +
  labs(title = "Distance Matrix: Child Growth Trajectories")

# From distance matrix, we can see some block structure forming, which suggests
# that there are some subgroups of children with similar growth trajectories
# Clustering is therefore useful in this case

# =============================================================================
# SECTION 9 — ASSESSING CLUSTERING TENDENCY & OPTIMAL K
# 
# We need to determine the optimal number of clusters k
# =============================================================================

cat("Running NbClust — this may take 1-2 minutes...\n")
nbclust_result <- NbClust(
  data     = growth_scaled,
  distance = "euclidean",
  min.nc   = 2,
  max.nc   = 6,       # max 6 clusters is reasonable for n=293
  method   = "kmeans",
  index    = "all"
)

optimal_k <- as.numeric(names(which.max(table(nbclust_result$Best.nc[1, ]))))
cat(sprintf("\nNbClust recommends k = %d clusters\n", optimal_k))

# Left plot shows how the Dindex score decreases as the number of clusters 
# increases from 2 to 6
# Right plot shows the highest peak at k = 3, which means the biggest 
# meaningful jump in cluster quality happens when going from 2 to 3 clusters
# NbClust recommends k=2 clusters

# =============================================================================
# Plot 1: Voting results 
votes    <- table(nbclust_result$Best.nc[1, ])
vote_df  <- as.data.frame(votes)
names(vote_df) <- c("k", "Votes")
vote_df$k <- as.numeric(as.character(vote_df$k))

plot_votes <- ggplot(vote_df, aes(x = k, y = Votes, fill = k == optimal_k)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.8) +
  geom_text(aes(label = Votes), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("grey70", "#E74C3C"), guide = "none") +
  scale_x_continuous(breaks = vote_df$k) +
  labs(title    = "A. NbClust Voting Results — zBMI Trajectory Clustering",
       subtitle = sprintf("Optimal k = %d (most index votes)", optimal_k),
       x = "Number of Clusters (k)", y = "Number of Indices Voting") +
  theme_classic(base_size = 14)

print(plot_votes)

# Each bar represents a possible number of clusters
# The red bar shows that k = 2 won with 11 votes

# =============================================================================
# Elbow method
wss <- sapply(1:6, function(k) {
  if (k == 1) sum(scale(growth_scaled, scale = FALSE)^2)
  else kmeans(growth_scaled, centers = k, nstart = 25)$tot.withinss
})

elbow_df <- data.frame(k = 1:6, WSS = wss)

plot_elbow <- ggplot(elbow_df, aes(x = k, y = WSS)) +
  geom_line(color = "#3498DB", linewidth = 1.5) +
  geom_point(color = "#3498DB", size = 4) +
  geom_point(data = elbow_df[elbow_df$k == optimal_k, ],
             aes(x = k, y = WSS), color = "#E74C3C", size = 6) +
  scale_x_continuous(breaks = 1:6) +
  labs(title    = "B. Elbow Method",
       subtitle = "Look for the bend — big drop then levels off",
       x = "Number of Clusters (k)",
       y = "Total Within-Cluster Sum of Squares") +
  theme_classic(base_size = 14)

print(plot_elbow)

# Here, the bend is clearly at k = 2, meaning the the elbow plot 
# agrees with NbClust

# =============================================================================
# Silhouette method 
# Silhouette width measures how well each point fits its own cluster vs others.

fviz_nbclust(growth_scaled, kmeans, method = "silhouette") +
  labs(title    = "C. Silhouette Method",
       subtitle = "Higher average silhouette = better-defined clusters") +
  theme_classic(base_size = 14)

# 0.50 = strong cluster structure
# 0.25–0.50 = weak but present cluster structure
# < 0.25 = no substantial structure

# Here, the highest peak is at k = 2 (~0.3), meaning there is a weak, but 
# present cluster structure

# =============================================================================
# SECTION 10 — CLUSTERING
#
# We apply the three methods above and compare them.
# =============================================================================

# K-means clustering

km_result <- kmeans(growth_scaled, centers = optimal_k, nstart = 25)

cat(sprintf("\nK-means with k = %d:\n", optimal_k))
print(km_result)

# Add cluster assignment back to the cleaned data
clean_data$kmeans_cluster <- factor(km_result$cluster)

# Mean of each growth variable per cluster — tells us what each cluster IS
cat("\nMean growth values per k-means cluster:\n")
print(
  aggregate(growth_matrix, by = list(Cluster = km_result$cluster), mean) %>%
    mutate(across(where(is.numeric), ~ round(., 2)))
)

#  Cluster  WHO_zBMI_birth  WHO_zBMI_12m  zBMI_24m
#     1          1.05          1.87        0.92
#     2         -1.42         -0.62       -0.68

# =============================================================================
# Visualize
# PCA projection for visualization
pca_result <- prcomp(growth_scaled)
pca_df <- data.frame(
  PC1     = pca_result$x[, 1],
  PC2     = pca_result$x[, 2],
  Cluster = factor(km_result$cluster)
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "K-Means Clusters: Child zBMI Trajectories",
    subtitle = sprintf("k = %d | Each point = one child", optimal_k),
    x        = "Principal Component 1",
    y        = "Principal Component 2"
  ) +
  theme_classic(base_size = 13)

# K-means clustering with k=2 identified two distinct zBMI trajectory groups
# Clusters show meaningful overlap in the middle (expected with a silhouette
# score of 0.3)
# Children in cluster 2 tend to have higher zBMI trajectories overall

# =============================================================================
# PAM clustering
# PAM (Partitioning Around Medoids) is more robust to outliers than k-means
# because cluster centers are actual data points (medoids), not means.

fviz_nbclust(growth_scaled, pam, method = "silhouette") +
  labs(title = "Silhouette Method for PAM") +
  theme_classic()

pam_result <- pam(growth_scaled, k = optimal_k)
clean_data$pam_cluster <- factor(pam_result$cluster)

cat("\nPAM cluster medoids (most representative child in each cluster):\n")
print(pam_result$medoids)

fviz_cluster(pam_result,
             palette      = "Set1",
             ellipse.type = "t",
             repel        = TRUE,
             ggtheme      = theme_classic()) +
  labs(title    = "PAM Clusters: Child Growth Trajectories",
       subtitle = sprintf("k = %d | Medoid-based clustering", optimal_k))

# Clusters are pretty much identical to k-means results with overlap in the middle

# =============================================================================
# Hierarchical clustering 

res_dist <- dist(growth_scaled, method = "euclidean")
res_hc   <- hclust(d = res_dist, method = "complete")

# Dendrogram
fviz_dend(res_hc, cex = 0.4,
          main = "Hierarchical Clustering Dendrogram — Child Growth Trajectories")

# Cut into optimal_k groups
hc_groups <- cutree(res_hc, k = optimal_k)
clean_data$hc_cluster <- factor(hc_groups)

cat("\nHierarchical cluster sizes:\n")
print(table(hc_groups))

# The x-axis labels are unreadable because there are 293 children 
# crammed into a small space
# The tree structure is clearly visible

# Colored dendrogram
fviz_dend(res_hc,
          k              = optimal_k,
          cex            = 0.4,
          k_colors       = RColorBrewer::brewer.pal(optimal_k, "Set2"),
          color_labels_by_k = TRUE,
          rect           = TRUE,
          main           = "Hierarchical Clusters — Child Growth Trajectories")


# =============================================================================
# SECTION 11 — CLUSTER VALIDATION
# =============================================================================

cat("\nRunning clValid — comparing hierarchical, k-means, and PAM...\n")

# Internal validation (from 4_clust_validation.R)
# Connectivity: lower = better
# Dunn index:   higher = better (compact, well-separated clusters)
# Silhouette:   higher = better

clmethods <- c("hierarchical", "kmeans", "pam")

intern_valid <- clValid(
  obj      = growth_scaled,
  nClust   = 2:5,
  clMethods = clmethods,
  validation = "internal"
)

cat("\n--- Internal Validation Summary ---\n")
print(summary(intern_valid))

# Output - Optimal Scores:
#                Score        Method       Clusters
# Connectivity  2.9290      hierarchical       2       
# Dunn          0.5246      hierarchical       2       
# Silhouette    0.6883      hierarchical       2      

# Stability validation
# APN: lower = better (consistent clusters when columns removed)
# AD:  lower = better
# ADM: lower = better
# FOM: lower = better

# nClust = 2:5; we test up to 5 for an extra layer of validation, even
# if we established k = 2

stab_valid <- clValid(
  obj       = growth_scaled,
  nClust    = 2:5,
  clMethods = clmethods,
  validation = "stability"
)

cat("\n--- Stability Validation — Optimal Scores ---\n")
print(optimalScores(stab_valid))

# Output -     Score      Method        Clusters
#   APN     0.03166109  hierarchical        2
#   AD      1.45934211           pam        5
#   ADM     0.11252769  hierarchical        2
#   FOM     0.76721028           pam        5

# Based on these results, choose the best method and k for Section 7.
# Update BEST_METHOD and BEST_K below if the validation results differ
# from the NbClust recommendation.

BEST_METHOD <- "kmeans"   # <-- update based on validation results
BEST_K      <- optimal_k  # <-- update if validation suggests different k

# Add the best cluster assignment as a dedicated column
if (BEST_METHOD == "kmeans") {
  clean_data$best_cluster <- clean_data$kmeans_cluster
} else if (BEST_METHOD == "pam") {
  clean_data$best_cluster <- clean_data$pam_cluster
} else {
  clean_data$best_cluster <- clean_data$hc_cluster
}

cat(sprintf("\nUsing '%s' with k=%d for downstream analysis.\n",
            BEST_METHOD, BEST_K))
cat("Cluster sizes:\n")
print(table(clean_data$best_cluster))


# =============================================================================
# SECTION 12 — ANSWERING THE RESEARCH QUESTION
#
# Research Question:
#   Do children of mothers with higher UPF scores belong to
#   higher-risk (higher zBMI) growth clusters?
# =============================================================================

# Profile each cluster
# Describe what each cluster looks like in terms of zBMI and weight gain.
# This tells us which clusters are "high-risk" vs "normal" growth.

cat("\n=== SECTION 7: LINKING CLUSTERS TO MATERNAL UPF SCORE ===\n")

cluster_profile <- clean_data %>%
  group_by(best_cluster) %>%
  summarise(
    n                  = n(),
    mean_zBMI_24m      = round(mean(zBMI_24m), 3),
    sd_zBMI_24m        = round(sd(zBMI_24m), 3),
    mean_birth_weight  = round(mean(Birth_weight_g), 1),
    mean_weight_6m     = round(mean(Weight_6m_kg), 2),
    mean_weight_12m    = round(mean(Weight_12m_kg), 2),
    mean_weight_24m    = round(mean(Weight_24m_kg), 2),
    mean_UPF_score     = round(mean(Ultra_processed_score), 3),
    pct_stunted        = round(mean(Stunted_24m) * 100, 1)
  )

cat("\nCluster profiles:\n")
print(cluster_profile)

# Does UPF score differ between clusters? 
# One-way ANOVA: is maternal UPF score significantly different across clusters?

upf_anova <- aov(Ultra_processed_score ~ best_cluster, data = clean_data)
cat("\nANOVA: Maternal UPF score by cluster:\n")
print(summary(upf_anova))

# If ANOVA is significant (p < 0.05), run Tukey post-hoc to see which
# pairs of clusters differ
if (summary(upf_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
  cat("\nTukey HSD post-hoc test:\n")
  print(TukeyHSD(upf_anova))
}

# Visualize UPF score by cluster (boxplot)

plot_upf_cluster <- ggplot(clean_data,
                            aes(x = best_cluster, y = Ultra_processed_score,
                                fill = best_cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  labs(
    title    = "Maternal UPF Score by zBMI Growth Cluster",
    subtitle = "Higher cluster number = higher zBMI trajectory (if ordered by zBMI)",
    x        = "Growth Cluster",
    y        = "Maternal Ultra-Processed Food Score (0-1)"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

print(plot_upf_cluster)
ggsave("plot_upf_by_cluster.png", plot_upf_cluster, width = 8, height = 5, dpi = 300)

# Maternal UPF score did not appear to differ meaningfully between the two 
# zBMI trajectory clusters, with both groups showing similar median scores and
# overlapping distributions. This suggests that in this sample, maternal UPF 
# consumption alone may not be sufficient to predict child growth cluster 
# membership.

# Visualize zBMI by cluster (boxplot) 

plot_zbmi_cluster <- ggplot(clean_data,
                             aes(x = best_cluster, y = zBMI_24m,
                                 fill = best_cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  geom_hline(yintercept = 0,  linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 2,  linetype = "dotted", color = "orange",
             linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = -2, linetype = "dotted", color = "blue",
             linewidth = 1, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "zBMI at 24 Months by Growth Cluster",
    subtitle = "Orange dotted = overweight threshold (z=2); Blue = underweight (z=-2)",
    x        = "Growth Cluster",
    y        = "zBMI at 24 Months"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

print(plot_zbmi_cluster)
ggsave("plot_zbmi_by_cluster.png", plot_zbmi_cluster, width = 8, height = 5, dpi = 300)

# The point at ~7.4 is concerning. This should have been removed 
# since it's above the WHO cut-off of +5
# Check code earlier  

# Weight trajectory plot by cluster
# This is the key visualization for "trajectory" — shows how each cluster
# grows differently from birth to 24 months.

traj_long <- clean_data %>%
  select(ID, best_cluster,
         Birth_weight_g, Weight_6m_kg, Weight_12m_kg, Weight_24m_kg) %>%
  mutate(Birth_weight_kg = Birth_weight_g / 1000) %>%
  select(-Birth_weight_g) %>%
  pivot_longer(
    cols      = c(Birth_weight_kg, Weight_6m_kg, Weight_12m_kg, Weight_24m_kg),
    names_to  = "timepoint",
    values_to = "weight_kg"
  ) %>%
  mutate(age_months = case_when(
    timepoint == "Birth_weight_kg"  ~ 0,
    timepoint == "Weight_6m_kg"     ~ 6,
    timepoint == "Weight_12m_kg"    ~ 12,
    timepoint == "Weight_24m_kg"    ~ 24
  ))

traj_summary <- traj_long %>%
  group_by(best_cluster, age_months) %>%
  summarise(mean_weight = mean(weight_kg),
            se_weight   = sd(weight_kg) / sqrt(n()),
            .groups = "drop")

plot_trajectory <- ggplot(traj_summary,
                           aes(x = age_months, y = mean_weight,
                               color = best_cluster, group = best_cluster)) +
  geom_line(linewidth = 1.8) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_weight - se_weight,
                    ymax = mean_weight + se_weight),
                width = 0.8, linewidth = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  scale_x_continuous(breaks = c(0, 6, 12, 24)) +
  labs(
    title    = "Mean Weight Trajectory by Growth Cluster (Birth to 24 Months)",
    subtitle = "Error bars = ±1 SE  |  Each cluster represents a distinct growth pattern",
    x        = "Age (months)",
    y        = "Mean Weight (kg)"
  ) +
  theme_classic(base_size = 14)

print(plot_trajectory)
ggsave("plot_weight_trajectory_by_cluster.png", plot_trajectory,
       width = 9, height = 6, dpi = 300)

# Does UPF score predict cluster membership?
# This directly answers: "does maternal dietary quality PREDICT cluster membership?"
# We use multinomial logistic regression because cluster is a categorical variable
# with more than 2 levels.
# Predictors: UPF score + covariates (Maternal_BMI, Sex, Gestational_age_weeks)

# Set reference cluster (the one with lowest mean zBMI = "normal" growth)
ref_cluster <- cluster_profile %>%
  filter(mean_zBMI_24m == min(mean_zBMI_24m)) %>%
  pull(best_cluster) %>%
  as.character()

clean_data$best_cluster <- relevel(clean_data$best_cluster, ref = ref_cluster)

cat(sprintf("\nMultinomial logistic regression: predicting cluster membership\n"))
cat(sprintf("Reference cluster: %s (lowest mean zBMI = 'normal growth')\n\n", ref_cluster))

multinom_model <- nnet::multinom(
  best_cluster ~ Ultra_processed_score + Maternal_BMI + Sex +
    Gestational_age_weeks + Household_income_index,
  data  = clean_data,
  trace = FALSE
)

cat("Model summary:\n")
print(summary(multinom_model))

# Compute z-scores and p-values (multinom doesn't give p-values directly)
z_scores <- summary(multinom_model)$coefficients /
  summary(multinom_model)$standard.errors
p_values <- (1 - pnorm(abs(z_scores))) * 2

cat("\nP-values for each predictor:\n")
print(round(p_values, 4))

cat("\nOdds ratios (exp of coefficients):\n")
print(round(exp(coef(multinom_model)), 3))

# Summary heatmap: cluster profiles
# Shows how all key variables compare across clusters

profile_long <- cluster_profile %>%
  select(best_cluster, mean_zBMI_24m, mean_UPF_score,
         mean_weight_24m, pct_stunted) %>%
  pivot_longer(-best_cluster, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(scaled_value = as.numeric(scale(value))) %>%
  ungroup()

plot_heatmap <- ggplot(profile_long,
                        aes(x = best_cluster, y = variable, fill = scaled_value)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 0, name = "Scaled\nValue") +
  scale_y_discrete(labels = c("mean_UPF_score"   = "Mean UPF Score",
                               "mean_weight_24m"  = "Mean Weight 24m (kg)",
                               "mean_zBMI_24m"    = "Mean zBMI 24m",
                               "pct_stunted"      = "% Stunted")) +
  labs(
    title    = "Cluster Profile Heatmap",
    subtitle = "Red = higher value, Blue = lower value (scaled within each variable)",
    x        = "Growth Cluster",
    y        = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(size = 12, face = "bold"))

print(plot_heatmap)
ggsave("plot_cluster_profile_heatmap.png", plot_heatmap,
       width = 9, height = 5, dpi = 300)


# =============================================================================
# SECTION 13 — SAVE FINAL RESULTS
# =============================================================================

# Save final dataset with cluster assignments
write_csv(clean_data, "final_clustered_data.csv")

# Save cluster profiles table
write_csv(cluster_profile, "cluster_profiles_summary.csv")

# Print final summary
cat("\n")
cat(paste(rep("=", 65), collapse = ""), "\n")
cat("FINAL SUMMARY\n")
cat("Research Question: Does maternal UPF score predict zBMI\n")
cat("growth cluster membership from birth to 24 months?\n")
cat(paste(rep("=", 65), collapse = ""), "\n")
cat(sprintf("Final sample size (after cleaning): %d children\n", nrow(clean_data)))
cat(sprintf("Optimal clustering method: %s\n", BEST_METHOD))
cat(sprintf("Optimal number of clusters: %d\n", BEST_K))
cat("\nCluster profiles (mean values):\n")
print(cluster_profile)
cat("\nFiles saved:\n")
cat("  clean_precision_growth_dataset.csv\n")
cat("  final_clustered_data.csv\n")
cat("  cluster_profiles_summary.csv\n")
cat("  plot_upf_by_cluster.png\n")
cat("  plot_zbmi_by_cluster.png\n")
cat("  plot_weight_trajectory_by_cluster.png\n")
cat("  plot_cluster_profile_heatmap.png\n")
cat(paste(rep("=", 65), collapse = ""), "\n")

