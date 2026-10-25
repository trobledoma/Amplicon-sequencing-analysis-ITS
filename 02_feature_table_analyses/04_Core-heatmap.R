
# HEATMAP OF CORE FUNGAL COMMUNITY
# Genus-level shared taxa in Z-score
# Heatmap with shared taxa 


# 1. Load packages-----------

library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(grid)


# 1.1. Check package versions-----------

R.version.string
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("stringr")
packageVersion("pheatmap")
packageVersion("grid")


# 2. Define input and output paths-----------

# set up working directory
input_file <- "feature_table.csv"
output_dir <- "08_heatmap_results"

# Set to TRUE to include the unclassified category in the heatmap.

include_unclassified_in_heatmap <- FALSE

# Keep the original visual contrast range.

visual_limit <- 1

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# 3. Load input file-----------

if (!file.exists(input_file)) {
  stop("Input feature table not found. Check the file path.")
}

feature_table <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(feature_table) <- trimws(names(feature_table))

if (anyDuplicated(names(feature_table))) {
  stop("The input table contains duplicated column names.")
}

if (!"Genus" %in% names(feature_table)) {
  stop("The input feature table must contain a Genus column.")
}


# 4. Identify sample columns and conditions-----------

# Sample columns are detected automatically.
# Supported examples: its-LC1, its-LRS2, its-LL3, LC1 and LRS1.
# Column names must follow a condition label plus replicate number.

sample_pattern <- stringr::regex(
  "^(?:its[-_.]?)?([A-Za-z][A-Za-z_-]*?)[0-9]+$",
  ignore_case = TRUE
)

sample_matches <- stringr::str_match(
  names(feature_table),
  sample_pattern
)

is_sample_column <- !is.na(sample_matches[, 1])

sample_columns <- names(feature_table)[is_sample_column]

if (length(sample_columns) == 0) {
  stop(
    "No sample columns were detected. ",
    "Check that sample headers contain a condition label and replicate number."
  )
}

condition_labels <- toupper(
  sample_matches[is_sample_column, 2]
)

condition_labels <- stringr::str_remove(
  condition_labels,
  "[-_.]+$"
)

# Create an internal sample metadata table.

sample_metadata <- data.frame(
  SampleID = sample_columns,
  Condition = condition_labels,
  stringsAsFactors = FALSE
)

# Preserve the order in which conditions first appear in the input table.

condition_order <- unique(sample_metadata$Condition)

if (length(condition_order) < 2) {
  stop("At least two conditions are required for a heatmap.")
}


# 5. Validate abundance data-----------

# Convert abundance columns to numeric.
# Empty abundance cells are interpreted as zero counts.

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


# 6. Clean genus names-----------

# Remove the g__ prefix.
# Unclassified records are retained as a separate category.

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

feature_table <- feature_table %>%
  dplyr::mutate(
    Genus_clean = clean_genus_name(.data$Genus),
    Genus_clean = dplyr::if_else(
      is.na(.data$Genus_clean),
      "Unclassified at genus level",
      .data$Genus_clean
    )
  )


# 7. Aggregate features at genus level-----------

# Sum counts from features assigned to the same genus.

genus_counts <- feature_table %>%
  dplyr::group_by(.data$Genus_clean) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Verify that aggregation preserves the original read totals.

aggregated_totals <- vapply(
  sample_columns,
  function(column) {
    sum(genus_counts[[column]], na.rm = TRUE)
  },
  numeric(1)
)

if (!isTRUE(
  all.equal(
    unname(original_totals),
    unname(aggregated_totals),
    tolerance = 1e-8
  )
)) {
  stop("Read totals changed during genus aggregation.")
}


# 8. Calculate relative abundance-----------

# Calculate relative abundance separately within each sample.

relative_abundance <- genus_counts %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ .x / sum(.x, na.rm = TRUE)
    )
  )

relative_abundance_sums <- vapply(
  sample_columns,
  function(column) {
    sum(relative_abundance[[column]], na.rm = TRUE)
  },
  numeric(1)
)

if (any(abs(relative_abundance_sums - 1) > 1e-8)) {
  stop("Relative abundances do not sum to one in every sample.")
}


# 9. Reshape relative abundance and join metadata-----------

relative_abundance_long <- relative_abundance %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_columns),
    names_to = "SampleID",
    values_to = "Relative_abundance"
  ) %>%
  dplyr::left_join(
    sample_metadata,
    by = "SampleID"
  )


# 10. Calculate mean relative abundance by condition-----------

# Average relative abundance across replicates within each condition.

genus_mean_abundance <- relative_abundance_long %>%
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


# 11. Create and export the grouped genus table-----------

# Create one row per genus and one column per condition.

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

# Verify that mean relative abundance sums to one per condition.

condition_totals <- vapply(
  condition_order,
  function(condition) {
    sum(otu_table_grouped[[condition]], na.rm = TRUE)
  },
  numeric(1)
)

if (any(abs(condition_totals - 1) > 1e-8)) {
  stop("Mean relative abundances do not sum to one per condition.")
}

write.csv(
  otu_table_grouped,
  file.path(output_dir, "otu_table_grouped-fungil.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 12. Select genera present in all conditions-----------

# Presence is defined as mean relative abundance greater than zero
# in every condition.

otu_shared <- otu_table_grouped %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(condition_order),
      ~ is.finite(.x) & .x > 0
    )
  )

# Exclude unclassified records by default.

if (!include_unclassified_in_heatmap) {
  
  otu_shared <- otu_shared %>%
    dplyr::filter(
      .data$Genus != "Unclassified at genus level"
    )
}

if (nrow(otu_shared) == 0) {
  stop(
    "No identified genera have positive mean relative abundance ",
    "in every condition."
  )
}

# Export the genera selected for the heatmap.

write.csv(
  otu_shared,
  file.path(output_dir, "genera_shared_across_conditions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 13. Prepare the numeric heatmap matrix-----------

heat_data <- as.matrix(
  otu_shared[, condition_order, drop = FALSE]
)

rownames(heat_data) <- otu_shared$Genus

storage.mode(heat_data) <- "numeric"

if (any(!is.finite(heat_data))) {
  stop("The heatmap matrix contains missing or non-finite values.")
}

if (anyDuplicated(rownames(heat_data))) {
  stop("Duplicate genus names remain in the heatmap matrix.")
}


# 14. Scale rows using Z-scores-----------

# Remove genera with no variation across conditions.

row_sd <- apply(heat_data, 1, stats::sd)

valid_rows <- is.finite(row_sd) & row_sd > 0

heat_data <- heat_data[
  valid_rows,
  ,
  drop = FALSE
]

if (nrow(heat_data) == 0) {
  stop(
    "No shared genera show variation across conditions; ",
    "row-wise Z-scores cannot be calculated."
  )
}

# Standardise each genus across conditions.

scaled_rows <- t(scale(t(heat_data)))

# Limit extreme Z-scores for visualisation.

scaled_rows[scaled_rows > visual_limit] <- visual_limit
scaled_rows[scaled_rows < -visual_limit] <- -visual_limit

# Export the scaled matrix for reproducibility.

scaled_table <- data.frame(
  Genus = rownames(scaled_rows),
  scaled_rows,
  check.names = FALSE
)

write.csv(
  scaled_table,
  file.path(output_dir, "heatmap_row_zscores.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 15. Define the pastel colour palette-----------

pastel_contrast <- grDevices::colorRampPalette(
  c("#C8A2C8", "#FFFFFF", "#8FD19E")
)(100)

heatmap_breaks <- seq(
  -visual_limit,
  visual_limit,
  length.out = 101
)


# 16. Generate heatmap with italic genus labels-----------

# Scale figure height according to the number of genera.

plot_height <- max(
  10,
  0.20 * nrow(scaled_rows) + 3
)

heatmap_plot <- pheatmap::pheatmap(
  scaled_rows,
  scale = "none",
  color = pastel_contrast,
  breaks = heatmap_breaks,
  border_color = NA,
  cellwidth = 55,
  fontsize_row = 12,
  fontsize_col = 12,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Genera Shared Across Conditions",
  silent = TRUE
)

# Format genus names in italics.

row_label_index <- which(
  heatmap_plot$gtable$layout$name == "row_names"
)

if (length(row_label_index) == 1) {
  
  row_label_grob <- heatmap_plot$gtable$grobs[[
    row_label_index
  ]]
  
  row_label_grob$gp$fontface <- "italic"
  
  heatmap_plot$gtable$grobs[[
    row_label_index
  ]] <- row_label_grob
}


# 17. Display and export the heatmap-----------

# Display the heatmap in RStudio.

grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)

# Save the heatmap as PDF.

grDevices::pdf(
  file.path(output_dir, "shared_genera_heatmap.pdf"),
  width = 10,
  height = plot_height,
  useDingbats = FALSE
)

grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)

grDevices::dev.off()

# Save the heatmap as high-resolution PNG.

grDevices::png(
  file.path(output_dir, "shared_genera_heatmap.png"),
  width = 2400,
  height = ceiling(plot_height * 180),
  res = 180
)

grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)

grDevices::dev.off()


# 18. Output-----------

# All output files are saved in heatmap_results/.
#
# Grouped genus table:
# - otu_table_grouped-fungil.csv
#
# Shared genera:
# - genera_shared_across_conditions.csv
#
# Heatmap data:
# - heatmap_row_zscores.csv
#
# Figures:
# - shared_genera_heatmap.pdf
# - shared_genera_heatmap.png
#
# The heatmap includes genera with positive mean relative abundance
# in every detected condition.
# Unclassified records are excluded by default.
# Row labels are displayed in italics.

