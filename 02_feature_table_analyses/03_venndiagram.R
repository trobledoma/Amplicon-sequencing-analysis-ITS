# VENN DIAGRAMA FOR CORE FUNGAL COMMUNITY
# Genus-level shared taxa in presence absence
# Venn diagrama with shared taxa in number and percentages


# 1. Load packages-----------

library(dplyr)
library(tidyr)
library(stringr)
library(VennDiagram)
library(grid)


# 1.1. Check package versions-----------

R.version.string
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("stringr")
packageVersion("VennDiagram")
packageVersion("grid")


# 2. Define input and output paths-----------

# set up working directory

input_file <- "feature_table.csv"

output_dir <- "07_venn_results"


# Set to TRUE to include the unclassified category in the Venn diagram.
# FALSE excludes it from the genus comparison.

include_unclassified_in_venn <- FALSE

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

# Detect sample names containing a condition label followed by replicate numbers.
# Examples: its-LC1, its-LRS2, its-LL3.

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
  stop("No sample columns could be identified from the column names.")
}

condition_labels <- sample_matches[is_sample_column, 2]

# Remove trailing separators from condition labels, if present.

condition_labels <- stringr::str_remove(
  condition_labels,
  "[-_.]+$"
)

sample_metadata <- data.frame(
  SampleID = sample_columns,
  Condition = condition_labels,
  stringsAsFactors = FALSE
)

# Automatically define condition order from the input table.

condition_order <- unique(sample_metadata$Condition)


# 5. Validate abundance data-----------

# Convert sample abundances to numeric.
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


# 6. Clean genus names-----------

# Remove the g__ prefix and standardise unclassified records.
# Missing genus labels are not replaced by family or order labels.

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


# 7. Aggregate feature counts by genus-----------

# Sum reads from features assigned to the same genus.

genus_counts <- feature_table %>%
  dplyr::group_by(.data$Genus_clean) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Verify that aggregation preserves read totals.

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


# 8. Calculate relative abundance per sample-----------

# Convert genus counts to relative abundance within each sample.

genus_relative_abundance <- genus_counts %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ .x / sum(.x, na.rm = TRUE)
    )
  )

relative_abundance_sums <- vapply(
  sample_columns,
  function(column) {
    sum(genus_relative_abundance[[column]], na.rm = TRUE)
  },
  numeric(1)
)

if (any(abs(relative_abundance_sums - 1) > 1e-8)) {
  stop("Relative abundances do not sum to one in every sample.")
}


# 9. Calculate mean relative abundance by condition-----------

# Reshape the data and attach the condition inferred from sample names.

relative_abundance_long <- genus_relative_abundance %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_columns),
    names_to = "SampleID",
    values_to = "Relative_abundance"
  ) %>%
  dplyr::left_join(
    sample_metadata,
    by = "SampleID"
  )

# Calculate the mean relative abundance of each genus
# across samples belonging to the same condition.

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


# 10. Create the grouped genus table-----------

# Create one row per genus and one column per condition.

feature_table_grouped <- genus_mean_abundance %>%
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
  dplyr::arrange(Genus)

# Export the grouped table for reuse in downstream analyses.

write.csv(
  feature_table_grouped,
  file.path(output_dir, "feature_table_grouped-fungil.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 11. Convert abundance to presence and absence-----------

# A genus is present when its mean relative abundance is greater than zero.

feature_presence_absence <- feature_table_grouped %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(condition_order),
      ~ as.integer(.x > 0)
    )
  )

# Exclude unclassified records from the genus-level comparison by default.

if (include_unclassified_in_venn) {
  
  venn_taxa <- feature_presence_absence
  
} else {
  
  venn_taxa <- feature_presence_absence %>%
    dplyr::filter(
      Genus != "Unclassified at genus level"
    )
}


# 12. Summarise genus presence by condition-----------

group_presence_summary <- data.frame(
  Condition = condition_order,
  Number_of_genera_present = vapply(
    condition_order,
    function(condition) {
      sum(venn_taxa[[condition]], na.rm = TRUE)
    },
    integer(1)
  )
)

write.csv(
  group_presence_summary,
  file.path(output_dir, "genus_presence_summary_by_condition.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 13. Identify genera shared across all conditions-----------

shared_genera_presence <- venn_taxa %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(condition_order),
      ~ .x == 1
    )
  )

# Export the names of genera shared across all conditions.

shared_genera_names <- shared_genera_presence %>%
  dplyr::select(Genus)

write.csv(
  shared_genera_names,
  file.path(output_dir, "genera_shared_across_all_conditions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Export mean relative abundances for shared genera.

shared_genera_abundance <- feature_table_grouped %>%
  dplyr::semi_join(
    shared_genera_names,
    by = "Genus"
  )

write.csv(
  shared_genera_abundance,
  file.path(output_dir, "shared_genera_relative_abundance.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# 14. Create genus lists for the Venn diagrams-----------

venn_lists <- stats::setNames(
  lapply(
    condition_order,
    function(condition) {
      unique(
        venn_taxa$Genus[
          venn_taxa[[condition]] == 1
        ]
      )
    }
  ),
  condition_order
)

if (any(lengths(venn_lists) == 0)) {
  stop(
    "At least one condition has no genera available for the Venn diagram."
  )
}

pastel_colors <- c(
  "#6F9FC6",
  "#D8B2C6",
  "#EACF65"
)


# 15. Create Venn diagram: genus counts-----------

venn_plot_counts <- VennDiagram::venn.diagram(
  x = venn_lists,
  filename = NULL,
  fill = pastel_colors,
  alpha = 0.6,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.3,
  cat.fontfamily = "sans",
  main = "Genera Present in Each Condition",
  main.cex = 1.5,
  main.fontfamily = "sans"
)

# Display the count-based diagram in RStudio.

grid::grid.newpage()
grid::grid.draw(venn_plot_counts)

# Save the diagram as PDF.

grDevices::pdf(
  file.path(output_dir, "venn_genera_counts.pdf"),
  width = 8,
  height = 8
)

grid::grid.newpage()
grid::grid.draw(venn_plot_counts)
grDevices::dev.off()

# Save the diagram as a high-resolution PNG.

grDevices::png(
  file.path(output_dir, "venn_genera_counts.png"),
  width = 2400,
  height = 2400,
  res = 300
)

grid::grid.newpage()
grid::grid.draw(venn_plot_counts)
grDevices::dev.off()


# 16. Create Venn diagram: counts and percentages-----------

venn_plot_percentages <- VennDiagram::venn.diagram(
  x = venn_lists,
  filename = NULL,
  fill = pastel_colors,
  alpha = 0.6,
  
  # Numbers inside the circles.
  
  cex = 1.5,
  fontfamily = "sans",
  
  # Condition labels.
  
  cat.cex = 1.3,
  cat.fontfamily = "sans",
  
  # Main title.
  
  main = "Genera Shared Across Conditions",
  main.cex = 1.5,
  main.fontfamily = "sans",
  
  # Display counts and percentages.
  
  print.mode = c("raw", "percent"),
  sigdigs = 3
)

# Display the count-and-percentage diagram in RStudio.

grid::grid.newpage()
grid::grid.draw(venn_plot_percentages)

# Save the diagram as PDF.

grDevices::pdf(
  file.path(output_dir, "venn_genera_counts_percentages.pdf"),
  width = 8,
  height = 8
)

grid::grid.newpage()
grid::grid.draw(venn_plot_percentages)
grDevices::dev.off()

# Save the diagram as a high-resolution PNG.

grDevices::png(
  file.path(output_dir, "venn_genera_counts_percentages.png"),
  width = 2400,
  height = 2400,
  res = 300
)

grid::grid.newpage()
grid::grid.draw(venn_plot_percentages)
grDevices::dev.off()


# 17. Output-----------

# All output files are saved in venn_results/.
#
# Input:
# - merged_feature_taxonomy_fungil.csv
#
# Intermediate grouped table:
# - feature_table_grouped-fungil.csv
#
# Summary tables:
# - genus_presence_summary_by_condition.csv
# - genera_shared_across_all_conditions.csv
# - shared_genera_relative_abundance.csv
#
# Figures:
# - venn_genera_counts.pdf
# - venn_genera_counts.png
# - venn_genera_counts_percentages.pdf
# - venn_genera_counts_percentages.png
#
# The Venn diagrams use genus-level presence/absence.
# Unclassified records are excluded by default.
# The grouped table retains the unclassified abundance category.

