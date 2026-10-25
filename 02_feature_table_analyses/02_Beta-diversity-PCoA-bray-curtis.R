# BETA DIVERSITY
# PCoA 
# Bray-Curtis matrix distance

# 1. Load packages-----------

library(dplyr)
library(tidyr)
library(stringr)
library(vegan)
library(ggplot2)


# 1.1. Check package versions-----------

R.version.string
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("stringr")
packageVersion("vegan")
packageVersion("ggplot2")


# 2. Define input and output paths-----------

# set up working directory
# Input: ASV-level feature table with taxonomic annotations.
input_file <- "feature_table.csv"

# Metadata containing sample identifiers and experimental conditions.
metadata_file <- "sample_metadata.csv"

# Output directory.
output_dir <- "06_beta_diversity_results"


similarity_threshold <- 60 # may be changed

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# 3. Load input files-----------

if (!file.exists(input_file)) {
  stop("Feature table not found. Check the input file path.")
}

if (!file.exists(metadata_file)) {
  stop("Sample metadata not found. Check the input file path.")
}

feature_table <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata_raw <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(feature_table) <- trimws(names(feature_table))
names(metadata_raw) <- trimws(names(metadata_raw))


# 4. Validate input structure-----------

# The feature table must contain genus-level taxonomy.

if (!"Genus" %in% names(feature_table)) {
  stop("The feature table must contain a Genus column.")
}

# The metadata must contain sample identifiers and conditions.

required_metadata_columns <- c(
  "sample_name",
  "location"
)

if (!all(required_metadata_columns %in% names(metadata_raw))) {
  stop(
    "The metadata must contain sample_name and location columns."
  )
}

if (anyDuplicated(names(feature_table))) {
  stop("The feature table contains duplicated column names.")
}

# Standardise the metadata fields used in this analysis.

sample_metadata <- metadata_raw %>%
  dplyr::transmute(
    SampleID = trimws(as.character(.data$sample_name)),
    Condition = trimws(as.character(.data$location))
  )

if (
  anyNA(sample_metadata$SampleID) ||
  any(sample_metadata$SampleID == "") ||
  anyDuplicated(sample_metadata$SampleID)
) {
  stop("Sample identifiers must be present and unique.")
}

if (
  anyNA(sample_metadata$Condition) ||
  any(sample_metadata$Condition == "")
) {
  stop("Every sample must have a condition assigned.")
}


# 5. Identify sample columns automatically-----------

# Sample columns are identified by matching the metadata identifiers
# to the feature table headers. No sample names are hard-coded.

missing_samples <- setdiff(
  sample_metadata$SampleID,
  names(feature_table)
)

if (length(missing_samples) > 0) {
  stop(
    "Some metadata samples are missing from the feature table. ",
    "Check sample_name against the feature table column names."
  )
}

sample_columns <- names(feature_table)[
  names(feature_table) %in% sample_metadata$SampleID
]

if (length(sample_columns) < 3) {
  stop("At least three samples are required for a two-axis PCoA.")
}

# Align metadata rows to the feature table column order.

sample_metadata <- sample_metadata[
  match(sample_columns, sample_metadata$SampleID),
  ,
  drop = FALSE
]

# Generate anonymous sample labels for analysis and output.
# Original sample identifiers are used internally for matching only.

sample_metadata$AnalysisSample <- sprintf(
  "Sample_%03d",
  seq_len(nrow(sample_metadata))
)

condition_order <- unique(sample_metadata$Condition)

if (length(condition_order) < 2) {
  stop("At least two conditions are required for this analysis.")
}


# 6. Validate abundance data-----------

# Convert sample abundance columns to numeric.
# Empty cells are interpreted as zero counts.

for (column in sample_columns) {
  
  values <- trimws(
    as.character(feature_table[[column]])
  )
  
  values[
    !is.na(values) & values == ""
  ] <- NA_character_
  
  numeric_values <- suppressWarnings(
    as.numeric(values)
  )
  
  if (any(!is.na(values) & is.na(numeric_values))) {
    stop("A sample column contains non-numeric abundance values.")
  }
  
  numeric_values[is.na(numeric_values)] <- 0
  
  if (
    any(!is.finite(numeric_values)) ||
    any(numeric_values < 0)
  ) {
    stop("Abundance values must be finite and non-negative.")
  }
  
  feature_table[[column]] <- numeric_values
}

original_totals <- vapply(
  sample_columns,
  function(column) {
    sum(feature_table[[column]], na.rm = TRUE)
  },
  numeric(1)
)

if (
  any(!is.finite(original_totals)) ||
  any(original_totals <= 0)
) {
  stop("At least one sample has an invalid or zero read total.")
}


# 7. Prepare the ASV-level count matrix-----------

# Beta diversity is calculated from the original ASV counts.
# Each row represents an ASV and each sample is a column.
# Do not aggregate ASVs by genus before calculating Bray-Curtis.

count_matrix <- as.matrix(
  feature_table[, sample_columns, drop = FALSE]
)

storage.mode(count_matrix) <- "numeric"

# Remove features with zero abundance in every sample.

keep_features <- rowSums(count_matrix, na.rm = TRUE) > 0

count_matrix_filtered <- count_matrix[
  keep_features,
  ,
  drop = FALSE
]

if (nrow(count_matrix_filtered) == 0) {
  stop("No features with positive abundance remain.")
}

# Save the number of features retained for quality control.

feature_filter_summary <- data.frame(
  Initial_features = nrow(count_matrix),
  Features_retained = nrow(count_matrix_filtered),
  Features_removed = sum(!keep_features)
)

write.csv(
  feature_filter_summary,
  file.path(output_dir, "feature_filter_summary.csv"),
  row.names = FALSE
)


# 8. Create the grouped genus table-----------

# Clean genus names without replacing missing genera with higher
# taxonomic levels.

clean_genus_name <- function(x) {
  
  x <- stringr::str_squish(
    stringr::str_remove(
      as.character(x),
      stringr::regex("^g__", ignore_case = TRUE)
    )
  )
  
  unclassified <- (
    is.na(x) |
      x == "" |
      stringr::str_detect(
        x,
        stringr::regex(
          "incertae.*sedis|unclassified|uncultured|unknown|unassigned|unidentified|norank",
          ignore_case = TRUE
        )
      )
  )
  
  x[unclassified] <- NA_character_
  
  x
}

# Aggregate original feature counts by genus.

genus_counts <- feature_table %>%
  dplyr::mutate(
    Genus_clean = clean_genus_name(.data$Genus),
    Genus_clean = dplyr::if_else(
      is.na(.data$Genus_clean),
      "Unclassified at genus level",
      .data$Genus_clean
    )
  ) %>%
  dplyr::group_by(.data$Genus_clean) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Calculate relative abundance separately within each sample.

genus_relative_abundance <- genus_counts %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ .x / sum(.x, na.rm = TRUE)
    )
  )

# Calculate mean relative abundance by condition.

genus_relative_long <- genus_relative_abundance %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_columns),
    names_to = "SampleID",
    values_to = "Relative_abundance"
  ) %>%
  dplyr::left_join(
    sample_metadata %>%
      dplyr::select(SampleID, Condition),
    by = "SampleID"
  )

genus_mean_abundance <- genus_relative_long %>%
  dplyr::group_by(
    .data$Genus_clean,
    .data$Condition
  ) %>%
  dplyr::summarise(
    Mean_relative_abundance = mean(
      .data$Relative_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Create a wide table with one column per condition.

otu_table_grouped <- genus_mean_abundance %>%
  tidyr::pivot_wider(
    names_from = "Condition",
    values_from = "Mean_relative_abundance",
    values_fill = 0
  ) %>%
  dplyr::rename(
    Genus = Genus_clean
  ) %>%
  dplyr::select(
    Genus,
    dplyr::all_of(condition_order)
  ) %>%
  dplyr::arrange(.data$Genus)

write.csv(
  otu_table_grouped,
  file.path(output_dir, "otu_table_grouped-fungil.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 9. Calculate relative abundance for beta diversity-----------

# Normalise each feature column by its sample total.
# Features are rows; samples are columns.

relative_matrix <- sweep(
  count_matrix_filtered,
  2,
  colSums(count_matrix_filtered),
  FUN = "/"
)

# Replace original sample names with anonymous analysis labels.

colnames(relative_matrix) <- sample_metadata$AnalysisSample

# vegan expects samples in rows and features in columns.

community_matrix <- t(relative_matrix)

rownames(community_matrix) <- sample_metadata$AnalysisSample

if (any(!is.finite(community_matrix))) {
  stop("The community matrix contains non-finite values.")
}


# 10. Calculate Bray-Curtis dissimilarity-----------

bray_dist <- vegan::vegdist(
  community_matrix,
  method = "bray"
)

bray_matrix <- as.matrix(bray_dist)

# Export the distance matrix using anonymous sample identifiers.

bray_distance_export <- data.frame(
  AnalysisSample = rownames(bray_matrix),
  round(bray_matrix, 3),
  check.names = FALSE
)

write.csv(
  bray_distance_export,
  file.path(output_dir, "bray_curtis_distance_matrix.csv"),
  row.names = FALSE
)


# 11. Run PCoA-----------

pcoa_result <- stats::cmdscale(
  bray_dist,
  k = 2,
  eig = TRUE,
  add = TRUE
)

pcoa_points <- pcoa_result$points

positive_eigenvalues <- pcoa_result$eig[
  pcoa_result$eig > 0
]

if (length(positive_eigenvalues) == 0) {
  stop("PCoA did not return positive eigenvalues.")
}

positive_eigenvalue_sum <- sum(positive_eigenvalues)

pc1_pct <- 100 * pcoa_result$eig[1] /
  positive_eigenvalue_sum

pc2_pct <- 100 * pcoa_result$eig[2] /
  positive_eigenvalue_sum

pcoa_df <- data.frame(
  AnalysisSample = rownames(pcoa_points),
  PCoA1 = pcoa_points[, 1],
  PCoA2 = pcoa_points[, 2],
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    sample_metadata %>%
      dplyr::select(AnalysisSample, Condition),
    by = "AnalysisSample"
  )

write.csv(
  pcoa_df,
  file.path(output_dir, "pcoa_coordinates.csv"),
  row.names = FALSE
)


# 12. Plot PCoA-----------

# Assign colours dynamically to the conditions.

condition_colours <- stats::setNames(
  grDevices::hcl.colors(
    length(condition_order),
    palette = "Set 2"
  ),
  condition_order
)

p_pcoa <- ggplot2::ggplot(
  pcoa_df,
  ggplot2::aes(
    x = .data$PCoA1,
    y = .data$PCoA2,
    fill = .data$Condition
  )
) +
  ggplot2::geom_point(
    shape = 21,
    color = "grey30",
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = condition_colours
  ) +
  ggplot2::labs(
    title = "Principal Coordinates Analysis (PCoA)",
    subtitle = "Bray-Curtis dissimilarity based on feature-level relative abundance",
    x = sprintf("PCoA1 (%.1f%%)", pc1_pct),
    y = sprintf("PCoA2 (%.1f%%)", pc2_pct),
    fill = "Condition"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom"
  )


# Display the PCoA plot in RStudio for inspection.

print(p_pcoa)

# Export the PCoA plot.

ggplot2::ggsave(
  file.path(output_dir, "pcoa_bray_curtis.pdf"),
  plot = p_pcoa,
  width = 8,
  height = 6,
  units = "in"
)

ggplot2::ggsave(
  file.path(output_dir, "pcoa_bray_curtis.pdf"),
  plot = p_pcoa,
  width = 8,
  height = 6,
  units = "in"
)

ggplot2::ggsave(
  file.path(output_dir, "pcoa_bray_curtis.png"),
  plot = p_pcoa,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)


# 13. Calculate pairwise Bray-Curtis similarities-----------

# Convert dissimilarity to percentage similarity.
# Similarity (%) = 100 * (1 - Bray-Curtis dissimilarity).

similarity_matrix <- 100 * (1 - bray_matrix)

pair_indices <- utils::combn(
  seq_len(nrow(similarity_matrix)),
  2
)

pair_conditions <- sample_metadata$Condition

pair_similarity <- data.frame(
  Sample_1 = sample_metadata$AnalysisSample[
    pair_indices[1, ]
  ],
  Sample_2 = sample_metadata$AnalysisSample[
    pair_indices[2, ]
  ],
  Condition_1 = pair_conditions[
    pair_indices[1, ]
  ],
  Condition_2 = pair_conditions[
    pair_indices[2, ]
  ],
  Bray_Curtis_dissimilarity = bray_matrix[
    cbind(pair_indices[1, ], pair_indices[2, ])
  ],
  Similarity_pct = similarity_matrix[
    cbind(pair_indices[1, ], pair_indices[2, ])
  ],
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(
    dplyr::desc(.data$Similarity_pct)
  )

write.csv(
  pair_similarity,
  file.path(output_dir, "bray_curtis_pairwise_similarity.csv"),
  row.names = FALSE
)

# Export only sample pairs above the chosen similarity threshold.

similar_pairs <- pair_similarity %>%
  dplyr::filter(
    .data$Similarity_pct >= similarity_threshold
  )

write.csv(
  similar_pairs,
  file.path(output_dir, "bray_curtis_pairs_above_threshold.csv"),
  row.names = FALSE
)



# 15. Output-----------

# All results are saved in beta_diversity_results/.
#
# Input:
# - merged_feature_taxonomy_fungil.csv
# - sample_metadata.csv
#
# Grouped genus table for reuse:
# - otu_table_grouped-fungil.csv
#
# Quality control:
# - feature_filter_summary.csv
#
# Beta diversity:
# - bray_curtis_distance_matrix.csv
# - pcoa_coordinates.csv
# - bray_curtis_pairwise_similarity.csv
# - bray_curtis_pairs_above_threshold.csv
#
# Figures:
# - pcoa_bray_curtis.pdf
# - pcoa_bray_curtis.png

# Bray-Curtis distances are calculated from relative abundances
# at the original feature level.
# The grouped genus table is generated separately for reuse.
# Original sample identifiers are not exported.