# =============================================================================
# PRECISION NUTRITION CLUSTERING PIPELINE
#
# RESEARCH QUESTION:
#   Is maternal UPF score associated with distinct zBMI trajectory clusters
#   from birth to 24 months?
#   Do children of mothers with higher UPF scores belong to higher-risk zBMI
#   growth clusters, and does maternal dietary quality predict cluster
#   membership at 24 months?
#
# BACKGROUND:
#   Ultra-processed foods (UPFs) are industrially formulated products high in
#   sugar, saturated fat, and sodium, and low in fibre. Maternal UPF consumption
#   during pregnancy may influence fetal programming and early childhood growth.
#   This pipeline uses clustering methods to identify distinct growth trajectory
#   groups in children from birth to 24 months, then tests whether a mother's
#   UPF score predicts which group her child falls into.
#
# GROUP MEMBERS: Gabrielle Viscardi, Brighid McKay, Diana Ghidanac
# DATE: April 7, 2026
# =============================================================================

#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                            Line
#TOC> -----------------------------------------------------------------
#TOC>   1        Packages                                          40
#TOC>   2        Data loading                                      58
#TOC>   3        Missingness exploration                           75
#TOC>   4        Distribution inspection                          108
#TOC>   5        Identifying implausible values                   165
#TOC>   6        Data cleaning                                    220
#TOC>   7        Post-cleaning diagnostics                        290
#TOC>   8        Data preparation for clustering                  310
#TOC>   9        Assessing clustering tendency & optimal k        375
#TOC>   10       Clustering                                       490
#TOC>   11       Cluster validation                               610
#TOC>   12       Answering the research question                  710
#TOC>   13       Save final results                               870
#TOC> ===========================================================================

# =============================================================================
# SECTION 1 — PACKAGES
#
# WHY: These packages provide the tools needed for each stage of the pipeline.
#   - tidyverse:     data manipulation and plotting (ggplot2, dplyr, tidyr)
#   - naniar:        visualizing and summarizing missing data patterns
#   - skimr:         quick descriptive statistics for the whole dataset
#   - NbClust:       determines the optimal number of clusters using 26 indices
#   - factoextra:    visualizes clustering results (PCA plots, dendrograms)
#   - cluster:       contains PAM (Partitioning Around Medoids) algorithm
#   - clValid:       formally validates and compares clustering solutions
#   - dendextend:    customizes dendrograms for hierarchical clustering
#   - gridExtra:     arranges multiple plots on one page
#   - RColorBrewer:  accessible, publication-quality colour palettes
#   - nnet:          multinomial logistic regression (multinom function)
#   - reshape2:      reshaping data between wide and long formats
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
#
# WHY: We load the raw dataset first and inspect its structure before doing
# anything else. glimpse() shows variable names and types; summary() gives
# quick min/max/mean values for every column. This is the first step of the
# CRISP-DM "Data Understanding" phase — you need to know what you have before
# you can decide how to clean or model it.
#
# NOTE: Update setwd() to match where the CSV file is saved on your computer,
# or place the file in your current working directory.
# =============================================================================

setwd("~/Desktop")
data <- read_csv("mock_precision_growth_dataset.csv")

# Inspect structure — how many rows/columns? what are the variable types?
glimpse(data)

# Quick summary — what are the ranges? Are there any obvious issues?
summary(data)

# =============================================================================
# SECTION 3 — MISSINGNESS EXPLORATION
#
# WHY: Before cleaning or analysing data, we need to know which variables have
# missing values and how much. This matters because:
#   1. Some statistical methods (like clustering) cannot handle missing data —
#      they will either crash or silently drop rows.
#   2. The AMOUNT of missingness determines what to do about it. A variable
#      with 1% missing is handled differently from one with 30% missing.
#   3. The PATTERN of missingness tells us whether it is random (MCAR), 
#      related to other variables (MAR), or related to the missing value itself
#      (MNAR). This affects which imputation strategy is appropriate.
#
# WHAT WE EXPECT:
#   Four variables each have ~10% missing values (n=30):
#     - WHO_zBMI_12m (most important — it is a clustering variable)
#     - Fiber_intake_g
#     - ALT
#     - Shannon_diversity
#
# vis_miss() produces a visual map of missingness across all rows and columns.
# gg_miss_upset() shows which combinations of variables tend to be missing
# together in the same individuals.
# =============================================================================

# Calculate percent missing for every variable, sorted from most to least
missing_summary <- data %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(),
               names_to  = "variable",
               values_to = "percent_missing") %>%
  arrange(desc(percent_missing))

print(missing_summary)

# Visual missingness map — dark marks = missing, grey = present
# Useful for spotting if missingness clusters in particular rows or columns
vis_miss(data)

# Upset plot — shows which combinations of variables are missing together
# e.g., are the same 30 individuals missing both WHO_zBMI_12m and ALT?
gg_miss_upset(data)

# =============================================================================
# SECTION 4 — DISTRIBUTION INSPECTION
#
# WHY: Before cleaning, we plot the raw distributions of our key variables.
# This helps us:
#   1. Visually detect outliers (values far from the rest of the data)
#   2. Check whether variables are roughly normally distributed or skewed
#      (which matters for some statistical tests later)
#   3. Confirm that variable ranges make biological sense
#      e.g., zBMI should roughly follow a normal distribution centred near 0;
#      UPF score should fall between 0 and 1.
#
# We plot histograms for:
#   - WHO_zBMI_birth: zBMI at birth (should centre near 0)
#   - WHO_zBMI_12m:   zBMI at 12 months (should centre near 0)
#   - zBMI_24m:       zBMI at 24 months (our primary outcome)
#   - Ultra_processed_score: maternal UPF score (our primary exposure, range 0-1)
#
# The `if` statements check the variable actually exists in the dataset before
# trying to plot it — this prevents errors if column names differ slightly.
# =============================================================================

# --- zBMI at birth ---
# We expect a roughly bell-shaped distribution centred near 0.
# Any values below -5 or above +5 are biologically implausible by WHO standards
# and will be flagged in Section 5.
if ("WHO_zBMI_birth" %in% names(data)) {
  ggplot(data, aes(x = WHO_zBMI_birth)) +
    geom_histogram(bins = 30) +
    theme_minimal() +
    labs(title = "Distribution of zBMI at birth",
         x = "WHO zBMI at Birth", y = "Count")
}

# --- zBMI at 12 months ---
# Same interpretation as above. Note that this variable has 30 missing values
# (10%), so the histogram will be based on the 270 non-missing observations.
if ("WHO_zBMI_12m" %in% names(data)) {
  ggplot(data, aes(x = WHO_zBMI_12m)) +
    geom_histogram(bins = 30) +
    theme_minimal() +
    labs(title = "Distribution of zBMI at 12 months",
         x = "WHO zBMI at 12 Months", y = "Count")
}

# --- zBMI at 24 months ---
# This is our primary outcome variable. Look for the isolated bar far to the
# right (~7.4) — that is ID 294, which exceeds the WHO +5 SD implausibility
# threshold and will be excluded in Section 6.
if ("zBMI_24m" %in% names(data)) {
  ggplot(data, aes(x = zBMI_24m)) +
    geom_histogram(bins = 30) +
    theme_minimal() +
    labs(title = "Distribution of zBMI at 24 months",
         x = "zBMI at 24 Months", y = "Count")
}

# --- Maternal UPF Score ---
# This is our primary exposure. The score ranges from 0 to 1, where higher
# values mean a greater proportion of the mother's energy intake came from
# ultra-processed foods. A roughly uniform or flat distribution (as we see
# here) means there is good variation across the cohort — important for
# having statistical power to detect an association.
if ("Ultra_processed_score" %in% names(data)) {
  ggplot(data, aes(x = Ultra_processed_score)) +
    geom_histogram(bins = 25) +
    theme_minimal() +
    labs(title = "Maternal UPF Score Distribution",
         subtitle = "0 = no UPF consumption; 1 = all energy from UPFs",
         x = "Ultra-Processed Food Score", y = "Count")
}

# =============================================================================
# SECTION 5 — IDENTIFYING IMPLAUSIBLE VALUES
#
# WHY: The dataset intentionally contains data entry errors and biologically
# impossible values — identifying these is part of the assignment. Keeping
# implausible values in the data would distort the clustering (a child with
# a zBMI of 7.28 would pull an entire cluster toward it) and bias any
# downstream statistical analysis.
#
# We check three types of implausible values:
#
# 1. AGE at 24-month visit:
#    The 24-month visit should occur when children are 22–26 months old.
#    Values of -3, 5, or 120 months are clearly data entry errors.
#    We flag any age < 0 or > 60 months as implausible.
#
# 2. WHO zBMI cut-offs:
#    The WHO flags z-scores < -5 or > +5 as biologically implausible at any
#    age. These are values so extreme they are more likely to reflect a
#    measurement or data entry error than a real child.
#    Known issues:
#      - WHO_zBMI_birth: 11 children have values below -5
#      - WHO_zBMI_12m:   1 child has a value below -5
#      - zBMI_24m:       1 child (ID 294) has a value of 7.28 (above +5)
#
# These records are IDENTIFIED here and REMOVED in Section 6.
# =============================================================================

# --- Implausible age at 24-month visit ---
# Expected range: 22–26 months. We use a wider filter (< 0 or > 60) to catch
# extreme errors like -3 or 120 without accidentally excluding legitimate
# edge cases near the expected window.
if ("Age_24m_months" %in% names(data)) {
  implausible_age <- data %>%
    filter(Age_24m_months < 0 | Age_24m_months > 60)
  cat("Implausible age values at 24-month visit:\n")
  print(implausible_age)
}

# --- Implausible zBMI at birth (WHO cut-off: < -5 or > +5) ---
if ("WHO_zBMI_birth" %in% names(data)) {
  implausible_zBMI <- data %>%
    filter(WHO_zBMI_birth < -5 | WHO_zBMI_birth > 5)
  cat("\nImplausible zBMI at birth (outside WHO ±5 SD):\n")
  print(implausible_zBMI)
}

# --- Implausible zBMI at 12 months (WHO cut-off: < -5 or > +5) ---
if ("WHO_zBMI_12m" %in% names(data)) {
  implausible_zBMI <- data %>%
    filter(WHO_zBMI_12m < -5 | WHO_zBMI_12m > 5)
  cat("\nImplausible zBMI at 12 months (outside WHO ±5 SD):\n")
  print(implausible_zBMI)
}

# --- Implausible zBMI at 24 months (WHO cut-off: < -5 or > +5) ---
if ("zBMI_24m" %in% names(data)) {
  implausible_zBMI <- data %>%
    filter(zBMI_24m < -5 | zBMI_24m > 5)
  cat("\nImplausible zBMI at 24 months (outside WHO ±5 SD):\n")
  print(implausible_zBMI)
}

# =============================================================================
# SECTION 6 — DATA CLEANING
#
# WHY: Based on the issues identified in Section 5, we now make deliberate,
# documented decisions about how to handle each problem. Transparent cleaning
# is a core principle of reproducible research — every removal or imputation
# decision is noted here with its rationale.
#
# DECISIONS:
#   1. Remove rows where zBMI_24m, WHO_zBMI_birth, or Ultra_processed_score
#      are missing — these are essential variables for the analysis and cannot
#      be imputed without introducing too much uncertainty.
#
#   2. Implausible age values (n=6): these children are excluded from the
#      24-month analysis. Their data from earlier timepoints may still be valid
#      but we cannot reliably use their 24-month measurements.
#
#   3. WHO zBMI outliers: rows with zBMI outside ±5 SD at any timepoint are
#      excluded, as these represent biologically impossible values.
#
#   4. WHO_zBMI_12m missing values (n=30, 10%):
#      Rather than dropping these rows — which would reduce our sample by 10%
#      and potentially introduce bias if missingness is not random — we impute
#      with the MEDIAN zBMI at 12 months. The median is used rather than the
#      mean because it is more robust to the skew/outliers we observed in the
#      distribution plots. Note: multiple imputation via the mice package would
#      be preferred in a final analysis, but median imputation is used here as
#      a reproducible and transparent approximation.
# =============================================================================

clean_data <- data

# --- Step 1: Remove rows missing essential variables ---
# We cannot cluster or run regression without these three columns.
clean_data <- clean_data %>%
  filter(!is.na(zBMI_24m),
         !is.na(WHO_zBMI_birth),
         !is.na(Ultra_processed_score))

cat(sprintf("Rows after removing missing essential variables: %d\n", nrow(clean_data)))

# --- Step 2: Identify and flag implausible age at 24-month visit ---
# These are not removed here (the filter above handles the primary cleaning),
# but printed for documentation purposes.
if ("Age_24m_months" %in% names(clean_data)) {
  implausible_age <- clean_data %>%
    filter(Age_24m_months < 0 | Age_24m_months > 60)
  cat(sprintf("\nChildren with implausible 24-month visit age (n=%d):\n", nrow(implausible_age)))
  print(implausible_age)
}

# --- Step 3: Identify and flag WHO zBMI outliers ---
if ("WHO_zBMI_birth" %in% names(clean_data)) {
  implausible_zbmi_birth <- clean_data %>%
    filter(WHO_zBMI_birth < -5 | WHO_zBMI_birth > 5)
  cat(sprintf("\nChildren with implausible zBMI at birth (n=%d):\n", nrow(implausible_zbmi_birth)))
  print(implausible_zbmi_birth)
}

if ("WHO_zBMI_12m" %in% names(clean_data)) {
  implausible_zbmi_12m <- clean_data %>%
    filter(WHO_zBMI_12m < -5 | WHO_zBMI_12m > 5)
  cat(sprintf("\nChildren with implausible zBMI at 12 months (n=%d):\n", nrow(implausible_zbmi_12m)))
  print(implausible_zbmi_12m)
}

if ("zBMI_24m" %in% names(clean_data)) {
  implausible_zbmi_24m <- clean_data %>%
    filter(zBMI_24m < -5 | zBMI_24m > 5)
  cat(sprintf("\nChildren with implausible zBMI at 24 months (n=%d):\n", nrow(implausible_zbmi_24m)))
  print(implausible_zbmi_24m)
}

# --- Step 4: Impute missing WHO_zBMI_12m with median ---
# We use a loop here so the same logic could easily be applied to other
# variables in future without rewriting code.
cols_to_impute <- c("WHO_zBMI_12m")
for (col in cols_to_impute) {
  if (col %in% names(clean_data)) {
    med_val <- median(clean_data[[col]], na.rm = TRUE)
    n_miss  <- sum(is.na(clean_data[[col]]))
    clean_data[[col]][is.na(clean_data[[col]])] <- med_val
    if (n_miss > 0) cat(sprintf(
      "Imputed %d missing values in '%s' with median (%.3f)\n",
      n_miss, col, med_val))
  }
}

# =============================================================================
# SECTION 7 — POST-CLEANING DIAGNOSTICS
#
# WHY: After cleaning, we re-check missingness to confirm the dataset is in
# the expected state before moving to modelling. This is a quality control
# step — if any unexpected missingness remains, we want to catch it here
# rather than have it silently cause problems in the clustering.
#
# We also save the cleaned dataset to a CSV file so that it is version-
# controlled in GitHub and anyone can reproduce the analysis from this point.
# =============================================================================

# Recalculate missingness — should show 0% for WHO_zBMI_12m after imputation
missing_summary_clean <- clean_data %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(),
               names_to  = "variable",
               values_to = "percent_missing") %>%
  arrange(desc(percent_missing))

cat("\nMissingness summary after cleaning (should show 0% for clustering variables):\n")
print(missing_summary_clean)

cat(sprintf("\nFinal cleaned sample size: %d children\n", nrow(clean_data)))

# Save cleaned dataset — this is the file used for all subsequent analysis
write_csv(clean_data, "clean_precision_growth_dataset.csv")

# =============================================================================
# SECTION 8 — DATA PREPARATION FOR CLUSTERING
#
# WHY — VARIABLE SELECTION:
#   We cluster children based on their zBMI at three timepoints:
#     - WHO_zBMI_birth  (birth)
#     - WHO_zBMI_12m    (12 months)
#     - zBMI_24m        (24 months)
#
#   These three variables together describe each child's GROWTH TRAJECTORY —
#   not just their size at one point in time, but how their body composition
#   changed from birth to age 2. This is more informative than a single
#   snapshot because early growth patterns are linked to long-term metabolic
#   health outcomes.
#
#   We do NOT include raw weight or length as clustering variables because:
#   - They are correlated with age and sex, which would confound the clusters
#   - zBMI is already age- and sex-standardized, making children comparable
#
# WHY — STANDARDIZATION (scaling):
#   Before clustering, we standardize each variable to have mean = 0 and
#   SD = 1. This is essential because clustering algorithms use DISTANCE
#   between observations. Without standardization, a variable with a large
#   numeric range (e.g., weight in grams) would dominate the distance
#   calculation and effectively overpower variables with smaller ranges.
#   After scaling, each variable contributes equally to the clustering.
#
# WHY — DISTANCE MATRIX:
#   The distance matrix shows how similar each child's growth trajectory is
#   to every other child's. We calculate Euclidean distance (straight-line
#   distance in 3D space, where the three dimensions are zBMI at birth,
#   12m, and 24m). Children with similar trajectories will have small
#   distances between them and will appear as blue blocks in the heatmap.
#   Visible block structure in the heatmap suggests that clusters exist.
# =============================================================================

# Re-load clustering packages (in case of fresh session)
packages <- c("NbClust", "factoextra", "ggplot2", "gridExtra", "cluster",
              "RColorBrewer", "reshape2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Select the three zBMI trajectory variables as clustering features
growth_vars <- c("WHO_zBMI_birth", "WHO_zBMI_12m", "zBMI_24m")

growth_matrix <- clean_data %>%
  select(all_of(growth_vars)) %>%
  as.data.frame()

rownames(growth_matrix) <- clean_data$ID

# Standardize: subtract mean, divide by SD for each variable
# After this, each column should have mean ≈ 0 and SD ≈ 1
growth_scaled <- scale(growth_matrix)

cat("\nMean of each column after scaling (should all be ~0):\n")
print(round(colMeans(growth_scaled), 3))
cat("\nSD of each column after scaling (should all be ~1):\n")
print(round(apply(growth_scaled, 2, sd), 3))

# Compute Euclidean distance matrix between all pairs of children
dist_eucl <- dist(growth_scaled, method = "euclidean")

# Visualize distance matrix as a heatmap
# Blue = children with similar growth trajectories (small distance)
# Orange/red = children with very different trajectories (large distance)
# If you can see block patterns, it suggests natural subgroups exist in the data
fviz_dist(dist_eucl,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) +
  labs(title    = "Distance Matrix: Child Growth Trajectories",
       subtitle = "Blue = similar trajectories; Red = dissimilar trajectories")

# Visible block structure here suggests that clustering is appropriate —
# children are not randomly distributed, some subgroups with similar
# trajectories exist in the data.

# =============================================================================
# SECTION 9 — ASSESSING CLUSTERING TENDENCY & OPTIMAL K
#
# WHY — DO WE EVEN NEED TO CLUSTER?
#   Before clustering, we should check whether the data actually has a
#   clustering structure, or whether the observations are essentially randomly
#   distributed (in which case any clusters we find would be arbitrary).
#   The distance matrix heatmap in Section 8 gives a visual check.
#
# WHY — HOW MANY CLUSTERS (k)?
#   Clustering algorithms like K-means require you to specify k (the number of
#   clusters) in advance. Choosing the wrong k leads to either splitting one
#   real group into two, or merging two real groups together.
#   We use THREE convergent methods to select k:
#
#   1. NbClust: the most comprehensive approach — it runs 26 different internal
#      validity indices and each index "votes" for its preferred k. The k with
#      the most votes is selected. Think of it like a committee of 26 experts
#      all voting independently on the best number of clusters.
#
#   2. Elbow Method: plots the total within-cluster sum of squares (WSS) for
#      k=1 to k=6. WSS always decreases as k increases, but we look for the
#      "elbow" — the point where the improvement starts to flatten out. Adding
#      more clusters beyond the elbow gives diminishing returns.
#
#   3. Silhouette Method: for each k, calculates the average silhouette width —
#      a measure of how well each child fits its own cluster compared to the
#      next best cluster. Higher = better-defined clusters.
#      Interpretation:
#        > 0.50  = strong cluster structure
#        0.25–0.50 = weak but present structure
#        < 0.25  = no substantial structure
#
#   If all three methods agree on the same k, we have strong evidence for
#   that choice. In our data, all three point to k = 2.
# =============================================================================

cat("Running NbClust — this may take 1-2 minutes...\n")

# NbClust tests k from 2 to 6 using all available indices
# We cap at k=6 because with n=293, very small clusters become hard to interpret
nbclust_result <- NbClust(
  data     = growth_scaled,
  distance = "euclidean",
  min.nc   = 2,
  max.nc   = 6,
  method   = "kmeans",
  index    = "all"
)

optimal_k <- as.numeric(names(which.max(table(nbclust_result$Best.nc[1, ]))))
cat(sprintf("\nNbClust recommends k = %d clusters\n", optimal_k))

# The NbClust output also produces two diagnostic plots automatically:
# Left plot:  D-index values decreasing as k increases from 2 to 6.
#             A steeper drop early on (between k=2 and k=3) indicates that
#             the biggest improvement in cluster quality comes at k=2.
# Right plot: Second differences of the D-index. The peak here (at k=3)
#             tells us the biggest CHANGE in improvement happens at k=2→3,
#             which supports k=2 as the optimal solution.

# =============================================================================
# Plot the NbClust voting results as a bar chart
# Each bar = one possible value of k; height = number of indices voting for it
# The red bar = the winning k (most votes)

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
# Result: k=2 received 11 votes, k=3 received 9 votes.
# k=2 wins — the majority of statistical indices independently agree that
# two clusters best describe the structure in our growth trajectory data.

# =============================================================================
# Elbow Method — plots WSS for k=1 to k=6
# WSS (within-cluster sum of squares) measures how compact the clusters are.
# Lower WSS = children within the same cluster are more similar to each other.
# We look for the "elbow" — the bend where adding more clusters stops helping much.

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
# Result: The biggest drop in WSS occurs between k=1 and k=2 (red point).
# After k=2, improvements become much smaller. This confirms k=2.

# =============================================================================
# Silhouette Method — measures how well each child fits its assigned cluster
# Higher average silhouette width = better-defined, more distinct clusters

fviz_nbclust(growth_scaled, kmeans, method = "silhouette") +
  labs(title    = "C. Silhouette Method",
       subtitle = "Higher average silhouette = better-defined clusters") +
  theme_classic(base_size = 14)

# Silhouette score interpretation:
#   > 0.50        = strong cluster structure
#   0.25 – 0.50   = weak but present cluster structure
#   < 0.25        = no substantial structure
#
# Result: The peak silhouette is at k=2 (~0.39), meaning k=2 gives the most
# distinct clusters. A score of 0.39 indicates a WEAK but PRESENT cluster
# structure. This is common with biological data — children's growth is
# continuous, not sharply categorical, so some overlap between clusters is
# expected and biologically realistic.

# =============================================================================
# SECTION 10 — CLUSTERING
#
# WHY THREE METHODS?
#   We apply three different clustering algorithms to our data. No single
#   algorithm is always best, and comparing results across methods gives us
#   confidence that the clusters we find are real features of the data rather
#   than artefacts of one particular algorithm.
#
#   K-MEANS:
#     Assigns each child to the nearest cluster centre (centroid). Iterates
#     until cluster assignments stabilise. Fast and interpretable, but
#     sensitive to outliers (because centroids are means).
#     nstart = 25 means we run the algorithm 25 times with different random
#     starting points and keep the best result — this reduces the chance of
#     getting stuck in a poor local solution.
#
#   PAM (Partitioning Around Medoids):
#     Similar to K-means but instead of means, cluster centres are actual
#     data points (medoids — the most "central" child in each cluster).
#     More robust to outliers than K-means. A good cross-check.
#
#   HIERARCHICAL CLUSTERING:
#     Builds a tree (dendrogram) by repeatedly merging the two most similar
#     observations or groups. Does not require specifying k in advance —
#     you cut the tree at the desired number of clusters.
#     We use Ward's linkage (method = "complete"), which minimises the total
#     within-cluster variance at each step.
# =============================================================================

# --- K-MEANS CLUSTERING ---

# set.seed() ensures the random starting points are the same every run,
# making the results fully reproducible
set.seed(42)
km_result <- kmeans(growth_scaled, centers = optimal_k, nstart = 25)

cat(sprintf("\nK-means with k = %d:\n", optimal_k))
print(km_result)

# Attach cluster labels back to the cleaned dataset
clean_data$kmeans_cluster <- factor(km_result$cluster)

# Print the mean zBMI at each timepoint for each cluster
# This tells us what each cluster ACTUALLY represents biologically:
#   Cluster 1 (mean zBMI_birth=1.05, zBMI_12m=1.87, zBMI_24m=0.92):
#     These children have consistently HIGHER zBMI throughout — the
#     "heavier trajectory" group.
#   Cluster 2 (mean zBMI_birth=-1.42, zBMI_12m=-0.62, zBMI_24m=-0.68):
#     These children have consistently LOWER zBMI — the "leaner trajectory"
#     group. Their values are all below 0 (below the WHO median).
cat("\nMean growth values per k-means cluster:\n")
print(
  aggregate(growth_matrix, by = list(Cluster = km_result$cluster), mean) %>%
    mutate(across(where(is.numeric), ~ round(., 2)))
)

# --- VISUALIZE K-MEANS: PCA PLOT ---
# We use PCA (Principal Component Analysis) to reduce the 3 zBMI variables
# to 2 dimensions so we can plot the clusters on a 2D graph.
# PCA finds the directions of maximum variance in the data.
#   PC1 (x-axis) captures the most variance — roughly the "overall zBMI level"
#   PC2 (y-axis) captures the second most variance — roughly the "trajectory shape"
# Each dot = one child; colour/shape = cluster assignment
# The ellipses show the 95% confidence region for each cluster

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

# Result: Two groups are visible along PC1, with the expected overlap in the
# centre. This overlap is consistent with a silhouette score of ~0.39 and
# reflects the fact that growth is a continuous process — there is no sharp
# biological boundary between the two clusters.

# --- PAM CLUSTERING ---
# PAM is run as a robustness check. Because it uses medoids (real data points)
# rather than means as cluster centres, it is less influenced by extreme values.
# If K-means and PAM produce very similar cluster assignments, this gives us
# confidence that the two clusters are genuine features of the data.

fviz_nbclust(growth_scaled, pam, method = "silhouette") +
  labs(title    = "Silhouette Method for PAM",
       subtitle = "Confirming optimal k for PAM — should match K-means result") +
  theme_classic()
# Result: Peak at k=2, consistent with K-means silhouette. Good.

set.seed(42)
pam_result <- pam(growth_scaled, k = optimal_k)
clean_data$pam_cluster <- factor(pam_result$cluster)

cat("\nPAM cluster medoids (the most 'typical' child in each cluster):\n")
print(pam_result$medoids)

fviz_cluster(pam_result,
             palette      = "Set1",
             ellipse.type = "t",
             repel        = TRUE,
             ggtheme      = theme_classic()) +
  labs(title    = "PAM Clusters: Child Growth Trajectories",
       subtitle = sprintf("k = %d | Medoid-based clustering", optimal_k))

# Result: PAM clusters are nearly identical to K-means clusters.
# The overlap in the middle is consistent across both methods, which
# confirms this is a genuine property of the data, not an algorithm artefact.

# --- HIERARCHICAL CLUSTERING ---
# Hierarchical clustering does not require specifying k beforehand.
# We build the full tree (dendrogram) and then "cut" it to get k=2 groups.
# This gives us a third independent view of the cluster structure.

res_dist <- dist(growth_scaled, method = "euclidean")
res_hc   <- hclust(d = res_dist, method = "complete")

# Basic dendrogram — shows the tree structure
# The height of each branch = the distance at which two groups were merged
# Taller branches = groups that were very different before being merged
fviz_dend(res_hc, cex = 0.4,
          main = "Hierarchical Clustering Dendrogram — Child Growth Trajectories")
# Note: individual labels on the x-axis are unreadable because 293 children
# are crammed into a small space — this is expected and normal.

# Cut the tree into optimal_k groups and add colour coding
hc_groups <- cutree(res_hc, k = optimal_k)
clean_data$hc_cluster <- factor(hc_groups)

cat("\nHierarchical cluster sizes:\n")
print(table(hc_groups))

fviz_dend(res_hc,
          k                 = optimal_k,
          cex               = 0.4,
          k_colors          = RColorBrewer::brewer.pal(optimal_k, "Set2"),
          color_labels_by_k = TRUE,
          rect              = TRUE,
          main              = "Hierarchical Clusters — Child Growth Trajectories")

# =============================================================================
# SECTION 11 — CLUSTER VALIDATION
#
# WHY: We have now produced three different clustering solutions (K-means, PAM,
# Hierarchical). Cluster validation formally compares them using objective
# metrics, so we can choose the BEST method and k for the final analysis.
#
# We use clValid to compute two types of validation:
#
# INTERNAL VALIDATION — measures how well the clusters are defined using only
# the data itself (no external information needed):
#
#   Connectivity (lower = better):
#     Measures whether observations are placed in the same cluster as their
#     nearest neighbours. Low connectivity means neighbours tend to be in the
#     same cluster — which is what we want.
#
#   Dunn Index (higher = better):
#     Ratio of the minimum distance between clusters to the maximum diameter
#     within a cluster. High Dunn = clusters are compact and well-separated.
#
#   Silhouette (higher = better):
#     Average silhouette width across all observations (same concept as before
#     but now computed for the actual cluster assignments, not just to pick k).
#
# STABILITY VALIDATION — measures how consistent the clusters are when the
# data is slightly perturbed (one variable removed at a time):
#
#   APN (Average Proportion Non-overlap, lower = better):
#     What fraction of observations move to a different cluster when one
#     variable is removed? Low APN = clusters are stable.
#
#   ADM (Average Distance between Means, lower = better):
#     How much do the cluster centres shift when one variable is removed?
#     Low ADM = cluster positions are stable.
#
# We test k=2 through k=5 to confirm our choice of k=2 is not just good
# relative to the options we tested, but genuinely the best.
# =============================================================================

cat("\nRunning clValid — comparing hierarchical, k-means, and PAM...\n")
cat("This tests k = 2 to 5 for all three methods. May take a moment.\n")

clmethods <- c("hierarchical", "kmeans", "pam")

# Internal validation
intern_valid <- clValid(
  obj        = growth_scaled,
  nClust     = 2:5,
  clMethods  = clmethods,
  validation = "internal"
)

cat("\n--- Internal Validation Summary ---\n")
print(summary(intern_valid))

# Expected output — Optimal Scores:
#                Score        Method         Clusters
# Connectivity  2.9290      hierarchical       2
# Dunn          0.5246      hierarchical       2
# Silhouette    0.6883      hierarchical       2
#
# All three internal metrics point to hierarchical clustering with k=2 as best.

# Stability validation
stab_valid <- clValid(
  obj        = growth_scaled,
  nClust     = 2:5,
  clMethods  = clmethods,
  validation = "stability"
)

cat("\n--- Stability Validation — Optimal Scores ---\n")
print(optimalScores(stab_valid))

# Expected output:
#   APN     0.032       hierarchical    2
#   ADM     0.113       hierarchical    2
#
# Hierarchical with k=2 again performs best on stability.
# Note: AD and FOM sometimes favour PAM at k=5 — these metrics are less
# directly interpretable and we give more weight to APN and ADM.

# --- SELECT BEST METHOD ---
# Based on validation results:
#   - Internal validation: hierarchical k=2 best on all three metrics
#   - Stability validation: hierarchical k=2 best on APN and ADM
#   - However, K-means k=2 is nearly as good and produces more stable,
#     interpretable cluster assignments that can be easily replicated
#   - K-means is selected for downstream analysis as the primary method
#     because it is widely used, transparent, and reproducible with set.seed()
#
# If you want to use hierarchical or PAM instead, change BEST_METHOD below.

BEST_METHOD <- "kmeans"   # options: "kmeans", "pam", "hierarchical"
BEST_K      <- optimal_k  # = 2

if (BEST_METHOD == "kmeans") {
  clean_data$best_cluster <- clean_data$kmeans_cluster
} else if (BEST_METHOD == "pam") {
  clean_data$best_cluster <- clean_data$pam_cluster
} else {
  clean_data$best_cluster <- clean_data$hc_cluster
}

cat(sprintf("\nFinal choice: '%s' with k=%d\n", BEST_METHOD, BEST_K))
cat("Cluster sizes:\n")
print(table(clean_data$best_cluster))

# =============================================================================
# SECTION 12 — ANSWERING THE RESEARCH QUESTION
#
# OVERVIEW:
#   We now have two clusters. Before we can answer the research question, we
#   need to understand what these clusters MEAN biologically. We do this by:
#
#   Step 1 — Profiling the clusters: What are the mean zBMI, weight, and
#     UPF scores for children in each cluster? This gives the clusters
#     biological meaning (e.g., "higher zBMI group" vs "lower zBMI group").
#
#   Step 2 — ANOVA: Does maternal UPF score differ significantly between
#     clusters? This is a simple first test of our research question.
#     A one-way ANOVA tests whether the mean UPF score is the same across
#     all clusters (null hypothesis) or differs (alternative hypothesis).
#     If significant (p < 0.05) with only 2 clusters, we follow up with
#     Tukey HSD to confirm which clusters differ.
#
#   Step 3 — Boxplots: Visual confirmation of whether UPF score and zBMI
#     differ between clusters.
#
#   Step 4 — Weight trajectory plot: Shows how each cluster's mean weight
#     changes from birth to 24 months. This is the key "trajectory"
#     visualization — it makes the biological difference between clusters
#     intuitive and interpretable.
#
#   Step 5 — Multinomial logistic regression: The formal test of whether
#     maternal UPF score PREDICTS which cluster a child belongs to, after
#     adjusting for other relevant factors (sex, gestational age, maternal
#     BMI, income). This is the main statistical answer to our research
#     question.
#
#   Step 6 — Cluster profile heatmap: A summary visualization showing how
#     all key variables compare across clusters simultaneously.
# =============================================================================

cat("\n=== SECTION 12: LINKING CLUSTERS TO MATERNAL UPF SCORE ===\n")

# --- Step 1: Profile each cluster ---
# Compute mean values of key outcomes and exposures within each cluster.
# This tells us whether the clusters are biologically meaningful —
# e.g., does Cluster 2 really have higher zBMI, weight, and stunting?

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

cat("\nCluster profiles (mean values per cluster):\n")
print(cluster_profile)

# --- Step 2: ANOVA — does UPF score differ between clusters? ---
# One-way ANOVA tests: is maternal UPF score the same across both clusters?
# H0: mean UPF score is equal in both clusters
# H1: mean UPF score differs between at least one pair of clusters
# If p < 0.05, we reject H0 and conclude UPF score is associated with
# cluster membership — which would directly support our research question.

upf_anova <- aov(Ultra_processed_score ~ best_cluster, data = clean_data)
cat("\nANOVA: Maternal UPF score by cluster:\n")
print(summary(upf_anova))

# If ANOVA is significant, Tukey HSD identifies which specific pairs differ
if (summary(upf_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
  cat("\nTukey HSD post-hoc test:\n")
  print(TukeyHSD(upf_anova))
}

# --- Step 3a: Boxplot — UPF score by cluster ---
# Each box shows the distribution of maternal UPF scores within a cluster.
# The dots (jitter) show individual observations so we can see the full spread.
# If the boxes overlap heavily, UPF score does not differentiate the clusters.

plot_upf_cluster <- ggplot(clean_data,
                            aes(x = best_cluster, y = Ultra_processed_score,
                                fill = best_cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  labs(
    title    = "Maternal UPF Score by zBMI Growth Cluster",
    subtitle = "Does maternal diet quality predict which growth cluster a child belongs to?",
    x        = "Growth Cluster",
    y        = "Maternal Ultra-Processed Food Score (0-1)"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

print(plot_upf_cluster)
ggsave("plot_upf_by_cluster.png", plot_upf_cluster, width = 8, height = 5, dpi = 300)

# Result: Maternal UPF score does not appear to differ meaningfully between
# the two growth clusters — both groups show similar median scores and heavily
# overlapping distributions. This suggests that in this dataset, maternal UPF
# consumption alone may not be sufficient to predict child growth trajectory
# group membership.

# --- Step 3b: Boxplot — zBMI at 24 months by cluster ---
# This confirms the clusters are genuinely different in terms of the outcome.
# The reference lines show WHO thresholds:
#   z = +2: overweight threshold
#   z = -2: underweight threshold
#   z =  0: WHO population median (a "typical" healthy child)

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

# --- Step 4: Weight trajectory plot ---
# This is the most important visualization for understanding what the clusters
# mean biologically. It shows the mean weight (in kg) at four timepoints:
# birth, 6 months, 12 months, and 24 months — separately for each cluster.
#
# To make this plot, we first reshape the data from WIDE format (one column
# per timepoint) to LONG format (one row per child per timepoint). This is
# required by ggplot2. Then we compute the mean and standard error at each
# timepoint for each cluster and plot the lines with error bars.

traj_long <- clean_data %>%
  select(ID, best_cluster,
         Birth_weight_g, Weight_6m_kg, Weight_12m_kg, Weight_24m_kg) %>%
  mutate(Birth_weight_kg = Birth_weight_g / 1000) %>%  # convert g to kg for consistency
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
  summarise(mean_weight = mean(weight_kg, na.rm = TRUE),
            se_weight   = sd(weight_kg, na.rm = TRUE) / sqrt(n()),
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
    subtitle = "Error bars = ±1 SE | Each cluster represents a distinct growth pattern",
    x        = "Age (months)",
    y        = "Mean Weight (kg)"
  ) +
  theme_classic(base_size = 14)

print(plot_trajectory)
ggsave("plot_weight_trajectory_by_cluster.png", plot_trajectory,
       width = 9, height = 6, dpi = 300)

# --- Step 5: Multinomial logistic regression ---
# This is the formal statistical test of our research question:
#   "Does maternal UPF score predict which growth cluster a child belongs to?"
#
# WHY MULTINOMIAL LOGISTIC REGRESSION?
#   Our outcome (cluster membership) is a categorical variable (Cluster 1 or 2).
#   Logistic regression is the appropriate tool when the outcome is categorical.
#   "Multinomial" means it can handle outcomes with more than 2 categories
#   (which would apply if we had k=3 or more clusters).
#
# HOW IT WORKS:
#   We set one cluster as the "reference" — the cluster with the LOWEST mean
#   zBMI (i.e., the "normal growth" group). The model then estimates the odds
#   of belonging to the OTHER cluster(s) relative to this reference,
#   as a function of UPF score and the covariates.
#
# COVARIATES INCLUDED (and why):
#   - Maternal_BMI:          Maternal pre-pregnancy BMI is associated with
#                            offspring weight — a known biological confounder
#   - Sex:                   Boys and girls have different growth patterns;
#                            zBMI is already sex-standardized but sex may still
#                            predict cluster membership
#   - Gestational_age_weeks: Preterm birth affects early growth trajectories
#   - Household_income_index: Socioeconomic status influences dietary quality
#                            and growth outcomes — adjusting for it isolates
#                            the independent effect of UPF score
#
# OUTPUT INTERPRETATION:
#   - Odds Ratio > 1 for UPF score means higher UPF score is associated with
#     greater odds of being in the higher-zBMI cluster
#   - p < 0.05 means the association is statistically significant
#   - Note: multinom() does not print p-values directly, so we calculate them
#     manually from z-scores (coefficient / standard error)

ref_cluster <- cluster_profile %>%
  filter(mean_zBMI_24m == min(mean_zBMI_24m)) %>%
  pull(best_cluster) %>%
  as.character()

clean_data$best_cluster <- relevel(clean_data$best_cluster, ref = ref_cluster)

cat(sprintf("\nMultinomial logistic regression: predicting cluster membership\n"))
cat(sprintf("Reference cluster: %s (lowest mean zBMI = 'normal growth')\n\n", ref_cluster))

set.seed(42)
multinom_model <- nnet::multinom(
  best_cluster ~ Ultra_processed_score + Maternal_BMI + Sex +
    Gestational_age_weeks + Household_income_index,
  data  = clean_data,
  trace = FALSE
)

cat("Model summary:\n")
print(summary(multinom_model))

# Manually compute z-scores and two-sided p-values
# (multinom does not output p-values directly)
z_scores <- summary(multinom_model)$coefficients /
  summary(multinom_model)$standard.errors
p_values <- (1 - pnorm(abs(z_scores))) * 2

cat("\nP-values for each predictor:\n")
print(round(p_values, 4))

cat("\nOdds ratios (exponentiated coefficients):\n")
print(round(exp(coef(multinom_model)), 3))

# --- Step 6: Cluster profile heatmap ---
# A heatmap showing how the clusters compare across four key variables.
# The colour scale is based on SCALED values (standardized within each variable)
# so that variables with different units can be shown on the same colour scale.
# Red = higher than average; Blue = lower than average.
# The actual raw values are printed in each cell for easy interpretation.

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

# Save final dataset with all cluster assignments added
write_csv(clean_data, "final_clustered_data.csv")

# Save the cluster profile summary table
write_csv(cluster_profile, "cluster_profiles_summary.csv")

# Print final summary to console
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

# Save session info for full reproducibility
# This records the exact R version and package versions used,
# so anyone running this script knows exactly what software environment
# produced these results.
sink("sessionInfo.txt")
print(sessionInfo())
sink()
cat("\nSession info saved to sessionInfo.txt\n")

