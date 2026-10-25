### RELATIVE ABUNDANCE PLOT
### CLASS LEVE
### TOP 20

# 1. Load packages-----------

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# 1.1. Check package versions-----------

R.version.string
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("stringr")
packageVersion("ggplot2")


# 2. Define input and output paths-----------

input_file <- "feature_table.csv"
metadata_file <- "sample_metadata.csv"
output_dir <- "05_top20_class_abundance_results"

number_of_classes <- 20

# Set to TRUE to include unclassified records among the candidates.

include_unclassified <- FALSE

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

sample_metadata_raw <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(feature_table) <- trimws(names(feature_table))
names(sample_metadata_raw) <- trimws(names(sample_metadata_raw))


# 4. Validate input structure-----------

if (!"Class" %in% names(feature_table)) {
  stop("The feature table must contain a Class column.")
}

if (!all(c("sample_name", "location") %in%
         names(sample_metadata_raw))) {
  stop(
    "The metadata must contain sample_name and location columns."
  )
}

# Prepare metadata using generic column names.

sample_metadata <- sample_metadata_raw %>%
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


# 5. Identify and validate sample columns-----------

# Match sample identifiers from the metadata to feature table columns.
# No individual sample names are hard-coded.

sample_columns <- sample_metadata$SampleID

if (!all(sample_columns %in% names(feature_table))) {
  stop(
    "Metadata sample identifiers must match feature table columns."
  )
}

if (nrow(sample_metadata) < 2) {
  stop("At least two samples are required.")
}

# Convert abundances to numeric.
# Empty cells are treated as zero counts.

for (column in sample_columns) {
  
  values <- trimws(as.character(feature_table[[column]]))
  
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


# 6. Clean class names-----------

# Remove the c__ prefix.
# Missing or unclassified taxonomy is grouped separately.

clean_class_name <- function(x) {
  
  x <- stringr::str_squish(
    stringr::str_remove(
      as.character(x),
      stringr::regex("^c__", ignore_case = TRUE)
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
    Class_clean = clean_class_name(.data$Class),
    Class_clean = dplyr::if_else(
      is.na(.data$Class_clean),
      "Unclassified at class level",
      .data$Class_clean
    )
  )


# 7. Aggregate ASV counts by class-----------

# Combine ASVs assigned to the same fungal class.

class_counts <- feature_table %>%
  dplyr::group_by(.data$Class_clean) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Verify that class aggregation preserves read totals.

aggregated_totals <- vapply(
  sample_columns,
  function(column) {
    sum(class_counts[[column]], na.rm = TRUE)
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
  stop("Read totals changed during class aggregation.")
}


# 8. Calculate relative abundance by sample-----------

# Relative abundance is calculated using the total reads
# in each original sample as the denominator.

class_relative_abundance <- class_counts %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ .x / sum(.x, na.rm = TRUE)
    )
  )


# 9. Reshape abundance data and join metadata-----------

class_relative_long <- class_relative_abundance %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_columns),
    names_to = "SampleID",
    values_to = "Relative_abundance"
  ) %>%
  dplyr::left_join(
    sample_metadata,
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    Relative_abundance_pct =
      100 * .data$Relative_abundance
  )


# 10. Identify the 20 most abundant classes-----------

# Rank classes by mean relative abundance across all samples.

class_ranking_data <- class_relative_long

if (!include_unclassified) {
  class_ranking_data <- class_ranking_data %>%
    dplyr::filter(
      .data$Class_clean != "Unclassified at class level"
    )
}

top_classes <- class_ranking_data %>%
  dplyr::group_by(.data$Class_clean) %>%
  dplyr::summarise(
    Mean_relative_abundance_pct = mean(
      .data$Relative_abundance_pct,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(.data$Mean_relative_abundance_pct)
  ) %>%
  dplyr::slice_head(n = number_of_classes) %>%
  dplyr::mutate(
    Rank = dplyr::row_number()
  ) %>%
  dplyr::select(
    Rank,
    Class = Class_clean,
    Mean_relative_abundance_pct
  )

if (nrow(top_classes) == 0) {
  stop("No classified fungal classes are available for plotting.")
}

write.csv(
  top_classes,
  file.path(output_dir, "top20_fungal_classes.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 11. Prepare plotting data-----------

# Retain individual samples for the selected classes.
# Sample identifiers are used internally but are not displayed.

plot_data <- class_relative_long %>%
  dplyr::filter(
    .data$Class_clean %in% top_classes$Class
  ) %>%
  dplyr::mutate(
    Feature_label = factor(
      .data$Class_clean,
      levels = rev(top_classes$Class)
    ),
    Condition = factor(
      .data$Condition,
      levels = unique(sample_metadata$Condition)
    )
  )

# Export a summary by class and condition without sample identifiers.

class_condition_summary <- plot_data %>%
  dplyr::group_by(
    .data$Class_clean,
    .data$Condition
  ) %>%
  dplyr::summarise(
    Mean_relative_abundance_pct = mean(
      .data$Relative_abundance_pct,
      na.rm = TRUE
    ),
    SD_relative_abundance_pct = stats::sd(
      .data$Relative_abundance_pct,
      na.rm = TRUE
    ),
    Number_of_samples = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    .data$Class_clean,
    .data$Condition
  )

write.csv(
  class_condition_summary,
  file.path(output_dir, "top20_class_abundance_by_condition.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 12. Define pastel colours for experimental conditions-----------

condition_palette <- c(
  "LC"  = "#98F5FF",
  "LRS" = "#EE6AA7",
  "LL"  = "lightgoldenrod1"
)

condition_colours <- condition_palette[condition_order]

if (anyNA(condition_colours)) {
  stop("Some conditions do not have an assigned colour.")
}

# 13. Calculate mean relative abundance across replicates-----------

# Retain the 20 most abundant fungal classes.

plot_data <- class_relative_long %>%
  dplyr::filter(
    .data$Class_clean %in% top_classes$Class
  ) %>%
  dplyr::group_by(
    .data$Class_clean,
    .data$Condition
  ) %>%
  dplyr::summarise(
    Mean_relative_abundance_pct = mean(
      .data$Relative_abundance_pct,
      na.rm = TRUE
    ),
    SD_relative_abundance_pct = stats::sd(
      .data$Relative_abundance_pct,
      na.rm = TRUE
    ),
    n_replicates = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Feature_label = factor(
      .data$Class_clean,
      levels = rev(top_classes$Class)
    ),
    Condition = factor(
      .data$Condition,
      levels = names(condition_colours)
    )
  )


# Export mean relative abundance by class and condition.

write.csv(
  plot_data %>%
    dplyr::select(
      Class = Class_clean,
      Condition,
      Mean_relative_abundance_pct,
      SD_relative_abundance_pct,
      n_replicates
    ),
  file.path(
    output_dir,
    "top20_class_mean_abundance_by_condition.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 14. Plot mean relative abundance in a single panel-----------

p_top_classes <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = Feature_label,
    y = Mean_relative_abundance_pct,
    fill = Condition
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.72
  ) +
  ggplot2::scale_fill_manual(
    values = condition_colours
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Relative Abundance of the 20 Most Abundant Fungal Classes",
    subtitle = "Bars represent mean relative abundance across replicates",
    x = "Fungal class",
    y = "Mean relative abundance (%)",
    fill = "Condition"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      size = 9,
      face = "italic"
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold")
  )


# 15. Display and export the figure-----------

# Display the plot in RStudio.

print(p_top_classes)

# Export as PDF.

ggplot2::ggsave(
  file.path(
    output_dir,
    "top20_fungal_classes_mean_abundance.pdf"
  ),
  plot = p_top_classes,
  width = 12,
  height = 9,
  units = "in"
)

# Export as high-resolution PNG.

ggplot2::ggsave(
  file.path(
    output_dir,
    "top20_fungal_classes_mean_abundance.png"
  ),
  plot = p_top_classes,
  width = 12,
  height = 9,
  units = "in",
  dpi = 300
)

