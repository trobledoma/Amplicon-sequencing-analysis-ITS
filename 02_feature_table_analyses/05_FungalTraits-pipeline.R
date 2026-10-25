# FUNGALTRAITS ANALYSIS
# Genus-level trophic modes and relative abundance
# Chord diagram for trophic mode and guid


# 1. Load packages-----------

library(fungaltraits)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(scales)
library(circlize)


# 1.1. Check package versions-----------

R.version.string
packageVersion("fungaltraits")
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("stringr")
packageVersion("ggplot2")
packageVersion("scales")
packageVersion("circlize")


# 2. Define input and output paths-----------
# set up working directory
input_file <- "feature_table.csv"
metadata_file <- "sample_metadata.csv"
output_dir <- "09_fungaltraits_results"

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

sample_metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(feature_table) <- trimws(names(feature_table))
names(sample_metadata) <- trimws(names(sample_metadata))


# Helper function to export CSV files

write_output_csv <- function(data, filename) {
  utils::write.csv(
    data,
    file = file.path(output_dir, filename),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}


# 4. Validate input structure-----------

# The feature table must contain a Genus column.

if (!"Genus" %in% names(feature_table)) {
  stop("The feature table must contain a Genus column.")
}

# The metadata must contain sample identifiers and locations.

required_metadata_columns <- c(
  "sample_name",
  "location"
)

if (!all(required_metadata_columns %in% names(sample_metadata))) {
  stop(
    "The metadata must contain sample_name and location columns."
  )
}

# Standardise metadata column names for downstream analysis.

sample_metadata <- sample_metadata %>%
  dplyr::rename(
    SampleID = sample_name,
    Condition = location
  ) %>%
  dplyr::mutate(
    SampleID = trimws(as.character(.data$SampleID)),
    Condition = trimws(as.character(.data$Condition))
  )

# Validate sample identifiers and conditions.

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

# Sample columns are identified through the metadata file.
# Individual sample identifiers are not hard-coded.

sample_columns <- sample_metadata$SampleID

if (!all(sample_columns %in% names(feature_table))) {
  stop(
    "Metadata sample identifiers must match feature table columns."
  )
}

# Convert abundances to numeric.
# Empty abundance cells are interpreted as zero counts.

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


# 6. Check original read totals-----------

# Save read totals before aggregating features by genus.

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


# 7. Aggregate features at genus level-----------

# Standardise genus names and group unclassified features.

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

# Sum feature counts for each genus.

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

relative_abundance <- genus_counts %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_columns),
      ~ .x / sum(.x, na.rm = TRUE)
    )
  )

# Verify that relative abundances sum to one per sample.

relative_abundance_sums <- vapply(
  sample_columns,
  function(column) {
    sum(relative_abundance[[column]], na.rm = TRUE)
  },
  numeric(1)
)

if (any(abs(relative_abundance_sums - 1) > 1e-8)) {
  stop("Relative abundances do not sum to one.")
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


# 10. Summarise mean relative abundance-----------

# Calculate mean relative abundance by genus and condition.

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
  ) %>%
  dplyr::rename(
    Genus = Genus_clean
  )

write_output_csv(
  genus_mean_abundance,
  "genus_mean_relative_abundance_by_condition.csv"
)


# 11. Load FungalTraits-----------

traits_raw <- fungaltraits::fungal_traits()

if (!all(c("Genus", "trophic_mode_fg") %in% names(traits_raw))) {
  stop("The FungalTraits database lacks required columns.")
}

# Standardise trophic mode labels.

normalize_trophic_mode <- function(x) {
  
  x <- stringr::str_squish(as.character(x))
  x[x == ""] <- NA_character_
  
  vapply(
    x,
    function(mode) {
      
      if (is.na(mode)) {
        return(NA_character_)
      }
      
      components <- strsplit(
        mode,
        "-",
        fixed = TRUE
      )[[1]]
      
      components <- sort(unique(trimws(components)))
      
      paste(components, collapse = "-")
    },
    character(1),
    USE.NAMES = FALSE
  )
}

traits <- traits_raw %>%
  dplyr::mutate(
    Genus_clean = clean_genus_name(.data$Genus),
    Genus_key = stringr::str_to_lower(.data$Genus_clean),
    Trophic_mode = normalize_trophic_mode(
      .data$trophic_mode_fg
    )
  )

# Export frequency tables for database annotations.

export_frequency_table <- function(
    data,
    column_name,
    output_column_name,
    filename
) {
  
  if (!column_name %in% names(data)) {
    return(invisible(NULL))
  }
  
  values <- as.character(data[[column_name]])
  
  values[
    is.na(values) | stringr::str_trim(values) == ""
  ] <- "Missing"
  
  frequency_table <- as.data.frame(
    table(values),
    stringsAsFactors = FALSE
  )
  
  names(frequency_table) <- c(
    output_column_name,
    "n_records"
  )
  
  frequency_table <- frequency_table[
    order(frequency_table$n_records, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  rownames(frequency_table) <- NULL
  
  write_output_csv(frequency_table, filename)
  
  invisible(frequency_table)
}

export_frequency_table(
  traits_raw,
  "trophic_mode_fg",
  "Trophic_mode",
  "database_frequency_trophic_modes.csv"
)

export_frequency_table(
  traits_raw,
  "guild_fg",
  "Guild",
  "database_frequency_guilds.csv"
)

export_frequency_table(
  traits_raw,
  "confidence_fg",
  "Confidence",
  "database_frequency_confidence.csv"
)

export_frequency_table(
  traits_raw,
  "taxonomic_level_fg",
  "Taxonomic_level",
  "database_frequency_taxonomic_levels.csv"
)


# 12. Resolve trophic modes by genus-----------

# Identify all trophic modes available in FungalTraits.

available_trophic_modes <- sort(unique(
  traits$Trophic_mode[
    !is.na(traits$Trophic_mode) &
      traits$Trophic_mode != ""
  ]
))

# Resolve multiple modes only when a compatible
# combined mode already exists in FungalTraits.

resolve_combined_mode <- function(modes) {
  
  modes <- unique(modes[!is.na(modes) & modes != ""])
  
  if (length(modes) == 0) {
    return(NA_character_)
  }
  
  if (length(modes) == 1) {
    return(modes[[1]])
  }
  
  components <- strsplit(
    modes,
    "-",
    fixed = TRUE
  )
  
  all_components <- unique(
    unlist(components, use.names = FALSE)
  )
  
  available_components <- strsplit(
    available_trophic_modes,
    "-",
    fixed = TRUE
  )
  
  compatible_modes <- available_trophic_modes[
    vapply(
      available_components,
      function(x) all(all_components %in% x),
      logical(1)
    )
  ]
  
  if (length(compatible_modes) == 0) {
    return(NA_character_)
  }
  
  # Prefer the smallest existing combination that covers all components.
  
  n_components <- lengths(
    strsplit(
      compatible_modes,
      "-",
      fixed = TRUE
    )
  )
  
  best_modes <- compatible_modes[
    n_components == min(n_components)
  ]
  
  if (length(best_modes) == 1) {
    return(best_modes[[1]])
  }
  
  NA_character_
}

# Summarise database records and unique modes by genus.

traits_by_genus <- traits %>%
  dplyr::filter(
    !is.na(.data$Genus_key),
    .data$Genus_key != ""
  ) %>%
  dplyr::group_by(.data$Genus_key) %>%
  dplyr::summarise(
    Database_genus = dplyr::first(.data$Genus_clean),
    n_database_records = dplyr::n(),
    
    Available_modes = list(
      sort(unique(
        .data$Trophic_mode[
          !is.na(.data$Trophic_mode) &
            .data$Trophic_mode != ""
        ]
      ))
    ),
    
    n_trophic_modes = dplyr::n_distinct(
      .data$Trophic_mode,
      na.rm = TRUE
    ),
    
    Trophic_modes = {
      modes <- sort(unique(
        .data$Trophic_mode[
          !is.na(.data$Trophic_mode) &
            .data$Trophic_mode != ""
        ]
      ))
      
      if (length(modes) == 0) {
        NA_character_
      } else {
        paste(modes, collapse = "; ")
      }
    },
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Resolved_mode = purrr::map_chr(
      .data$Available_modes,
      resolve_combined_mode
    )
  )

# Export genus-level trophic mode assignments.

write_output_csv(
  traits_by_genus %>%
    dplyr::select(
      Database_genus,
      n_database_records,
      n_trophic_modes,
      Trophic_modes,
      Resolved_mode
    ),
  "fungaltraits_trophic_modes_by_genus.csv"
)


# 13. Review trophic mode assignments-----------

mode_review <- traits_by_genus %>%
  dplyr::mutate(
    Decision = dplyr::case_when(
      .data$n_trophic_modes == 0 ~
        "No trophic mode assigned",
      
      .data$n_trophic_modes == 1 ~
        "Single trophic mode",
      
      .data$n_trophic_modes > 1 &
        !is.na(.data$Resolved_mode) ~
        "Resolved using an existing combined mode",
      
      TRUE ~ "Pending manual review"
    )
  ) %>%
  dplyr::transmute(
    Genus = .data$Database_genus,
    n_trophic_modes = .data$n_trophic_modes,
    Original_modes = .data$Trophic_modes,
    Resolved_mode = .data$Resolved_mode,
    Decision = .data$Decision
  )

write_output_csv(
  mode_review,
  "trophic_mode_assignment_review.csv"
)

# Summarise review decisions.

mode_review_summary <- mode_review %>%
  dplyr::count(
    .data$Decision,
    name = "n_genera"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$n_genera))

write_output_csv(
  mode_review_summary,
  "trophic_mode_review_summary.csv"
)

# Export genera requiring manual review.

pending_mode_review <- mode_review %>%
  dplyr::filter(
    .data$Decision == "Pending manual review"
  )

write_output_csv(
  pending_mode_review,
  "pending_trophic_mode_review.csv"
)


# 14. Match observed genera with FungalTraits-----------

# Preserve taxa that cannot be identified to genus level.

observed_genera <- genus_mean_abundance %>%
  dplyr::distinct(Genus) %>%
  dplyr::mutate(
    Genus_key = dplyr::if_else(
      .data$Genus == "Unclassified at genus level",
      NA_character_,
      stringr::str_to_lower(
        clean_genus_name(.data$Genus)
      )
    )
  )

genus_assignment <- observed_genera %>%
  dplyr::left_join(
    traits_by_genus %>%
      dplyr::select(
        Genus_key,
        n_database_records,
        n_trophic_modes,
        Trophic_modes,
        Resolved_mode
      ),
    by = "Genus_key"
  ) %>%
  dplyr::mutate(
    Assignment_status = dplyr::case_when(
      .data$Genus == "Unclassified at genus level" ~
        "Unclassified at genus level",
      
      is.na(.data$n_database_records) ~
        "Genus not found in FungalTraits",
      
      .data$n_trophic_modes == 0 ~
        "No trophic mode assigned",
      
      .data$n_trophic_modes == 1 ~
        "Single trophic mode",
      
      .data$n_trophic_modes > 1 &
        !is.na(.data$Resolved_mode) ~
        "Resolved using combined mode",
      
      .data$n_trophic_modes > 1 ~
        "Multiple modes / unresolved",
      
      TRUE ~ "Unknown trophic mode"
    ),
    
    Trophic_mode_final = dplyr::case_when(
      .data$Genus == "Unclassified at genus level" ~
        "Unclassified at genus level",
      
      is.na(.data$n_database_records) |
        .data$n_trophic_modes == 0 ~
        "Unknown trophic mode",
      
      .data$n_trophic_modes == 1 |
        (
          .data$n_trophic_modes > 1 &
            !is.na(.data$Resolved_mode)
        ) ~
        .data$Resolved_mode,
      
      .data$n_trophic_modes > 1 ~
        "Multiple / Unresolved",
      
      TRUE ~ "Unknown trophic mode"
    )
  )

write_output_csv(
  genus_assignment %>%
    dplyr::select(
      Genus,
      Assignment_status,
      Trophic_mode_final,
      Trophic_modes
    ),
  "observed_genus_trophic_mode_assignment.csv"
)


# 15. FungalTraits coverage by condition-----------

# Join genus-level assignments to mean abundance by condition.
# Individual sample identifiers are not exported.

genus_condition_assignment <- genus_mean_abundance %>%
  dplyr::left_join(
    genus_assignment %>%
      dplyr::select(
        Genus,
        Assignment_status,
        Trophic_mode_final
      ),
    by = "Genus"
  )

write_output_csv(
  genus_condition_assignment,
  "genus_abundance_and_trophic_assignment.csv"
)

# Calculate the percentage of total community abundance
# represented by each assignment status in each condition.

coverage_by_condition <- genus_condition_assignment %>%
  dplyr::group_by(
    .data$Condition,
    .data$Assignment_status
  ) %>%
  dplyr::summarise(
    Relative_abundance = sum(
      .data$Mean_relative_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(.data$Condition) %>%
  dplyr::mutate(
    Relative_abundance_pct =
      100 * .data$Relative_abundance /
      sum(.data$Relative_abundance, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    .data$Condition,
    .data$Assignment_status
  )

write_output_csv(
  coverage_by_condition,
  "fungaltraits_coverage_by_condition.csv"
)


# 16. Figure: FungalTraits coverage by condition-----------

coverage_plot <- ggplot2::ggplot(
  coverage_by_condition,
  ggplot2::aes(
    x = .data$Assignment_status,
    y = .data$Relative_abundance_pct,
    fill = .data$Condition
  )
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "FungalTraits Coverage by Condition",
    subtitle = "Genus-level trophic mode assignment status",
    x = "Assignment status",
    y = "Relative abundance (%)",
    fill = "Condition"
  ) +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave(
  file.path(output_dir, "fungaltraits_coverage_by_condition.png"),
  plot = coverage_plot,
  width = 11,
  height = 7,
  dpi = 300
)

ggplot2::ggsave(
  file.path(output_dir, "fungaltraits_coverage_by_condition.pdf"),
  plot = coverage_plot,
  width = 11,
  height = 7
)


# 17. Relative abundance by trophic mode-----------

# Aggregate mean relative abundance by trophic mode and condition.
# These categories include unclassified and unresolved assignments.

abundance_by_mode_condition <- genus_condition_assignment %>%
  dplyr::group_by(
    .data$Condition,
    .data$Trophic_mode_final
  ) %>%
  dplyr::summarise(
    Mean_relative_abundance = sum(
      .data$Mean_relative_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Relative_abundance_pct =
      100 * .data$Mean_relative_abundance
  ) %>%
  dplyr::arrange(
    .data$Condition,
    .data$Trophic_mode_final
  )

write_output_csv(
  abundance_by_mode_condition,
  "trophic_mode_abundance_by_condition.csv"
)

# Export a wide table with one column per condition.

abundance_by_mode_wide <- abundance_by_mode_condition %>%
  dplyr::select(
    Trophic_mode_final,
    Condition,
    Relative_abundance_pct
  ) %>%
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = Relative_abundance_pct,
    values_fill = 0
  )

write_output_csv(
  abundance_by_mode_wide,
  "trophic_mode_relative_abundance_pct_by_condition.csv"
)

# Plot relative abundance by trophic mode.

trophic_mode_plot <- ggplot2::ggplot(
  abundance_by_mode_condition,
  ggplot2::aes(
    x = .data$Trophic_mode_final,
    y = .data$Relative_abundance_pct,
    fill = .data$Condition
  )
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Relative Abundance by Trophic Mode",
    subtitle = "Percentages are relative to total community abundance",
    x = "Trophic mode",
    y = "Relative abundance (%)",
    fill = "Condition"
  ) +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave(
  file.path(output_dir, "trophic_mode_abundance.png"),
  plot = trophic_mode_plot,
  width = 11,
  height = 7,
  dpi = 300
)

ggplot2::ggsave(
  file.path(output_dir, "trophic_mode_abundance.pdf"),
  plot = trophic_mode_plot,
  width = 11,
  height = 7
)


# 18. Genera assigned to each trophic mode-----------

# Unknown, unclassified and unresolved genera are excluded
# from this assigned-genus summary only.

excluded_categories <- c(
  "Unknown trophic mode",
  "Unclassified at genus level",
  "Multiple / Unresolved"
)

genera_by_mode <- genus_assignment %>%
  dplyr::filter(
    !.data$Trophic_mode_final %in% excluded_categories
  ) %>%
  dplyr::distinct(
    Trophic_mode_final,
    Genus
  ) %>%
  dplyr::arrange(
    .data$Trophic_mode_final,
    .data$Genus
  )

write_output_csv(
  genera_by_mode,
  "genera_by_trophic_mode.csv"
)

genera_summary_by_mode <- genera_by_mode %>%
  dplyr::group_by(.data$Trophic_mode_final) %>%
  dplyr::summarise(
    n_genera = dplyr::n_distinct(.data$Genus),
    Genera = paste(
      sort(unique(.data$Genus)),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data$Trophic_mode_final)

write_output_csv(
  genera_summary_by_mode,
  "genera_summary_by_trophic_mode.csv"
)


# 19. Chord diagram: conditions and trophic modes-----------

# This diagram represents trophic modes only, not guilds.
#
# Excluded from the diagram:
# - Unknown trophic mode
# - Unclassified at genus level
# - Multiple / Unresolved
#
# These categories remain in the exported CSV tables.
# The percentages use total community abundance as denominator.
# Therefore, visible links may sum to less than 100%.

condition_order <- unique(
  as.character(sample_metadata$Condition)
)

chord_data <- abundance_by_mode_condition %>%
  dplyr::filter(
    !.data$Trophic_mode_final %in% excluded_categories,
    .data$Mean_relative_abundance > 0
  )

chord_mode_levels <- sort(unique(
  chord_data$Trophic_mode_final
))

if (length(chord_mode_levels) > 0) {
  
  chord_long <- tidyr::expand_grid(
    Condition = condition_order,
    Trophic_mode_final = chord_mode_levels
  ) %>%
    dplyr::left_join(
      chord_data %>%
        dplyr::transmute(
          Condition = .data$Condition,
          Trophic_mode_final = .data$Trophic_mode_final,
          Percent = .data$Relative_abundance_pct
        ),
      by = c("Condition", "Trophic_mode_final")
    ) %>%
    dplyr::mutate(
      Percent = tidyr::replace_na(.data$Percent, 0)
    )
  
  # Prepare the matrix: rows are conditions, columns are trophic modes.
  
  chord_wide <- chord_long %>%
    tidyr::pivot_wider(
      names_from = Trophic_mode_final,
      values_from = Percent,
      values_fill = 0
    ) %>%
    dplyr::select(
      Condition,
      dplyr::all_of(chord_mode_levels)
    ) %>%
    dplyr::arrange(
      match(.data$Condition, condition_order)
    )
  
  chord_matrix <- as.matrix(
    chord_wide[, chord_mode_levels, drop = FALSE]
  )
  
  rownames(chord_matrix) <- paste0(
    "Condition: ",
    chord_wide$Condition
  )
  
  colnames(chord_matrix) <- paste0(
    "Mode: ",
    chord_mode_levels
  )
  
  condition_sector_names <- rownames(chord_matrix)
  mode_sector_names <- colnames(chord_matrix)
  
  # Assign colours automatically to conditions and trophic modes.
  
  condition_colours <- stats::setNames(
    grDevices::hcl.colors(
      length(condition_sector_names),
      palette = "Set 2"
    ),
    condition_sector_names
  )
  
  mode_colours <- stats::setNames(
    grDevices::hcl.colors(
      length(mode_sector_names),
      palette = "Dynamic"
    ),
    mode_sector_names
  )
  
  sector_colours <- c(
    condition_colours,
    mode_colours
  )
  
  # Scale condition sectors to 100%.
  # Mode sectors scale to their total displayed abundance.
  
  xmax <- c(
    stats::setNames(
      rep(100, length(condition_sector_names)),
      condition_sector_names
    ),
    stats::setNames(
      colSums(chord_matrix),
      mode_sector_names
    )
  )
  
  link_colours <- matrix(
    rep(
      unname(mode_colours[colnames(chord_matrix)]),
      each = nrow(chord_matrix)
    ),
    nrow = nrow(chord_matrix),
    dimnames = dimnames(chord_matrix)
  )
  
  # Function to draw the diagram.
  
  draw_trophic_chord <- function() {
    
    circlize::circos.clear()
    
    on.exit(
      circlize::circos.clear(),
      add = TRUE
    )
    
    circlize::circos.par(
      start.degree = 90,
      gap.after = c(
        rep(8, length(condition_sector_names)),
        rep(2, length(mode_sector_names))
      )
    )
    
    circlize::chordDiagram(
      chord_matrix,
      xmax = xmax,
      grid.col = sector_colours,
      col = link_colours,
      directional = 0,
      transparency = 0.25,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.14),
      reduce = 0
    )
    
    graphics::title(
      main = "Trophic Modes and Relative Abundance by Condition (%)"
    )
    
    # Label the sectors.
    
    circlize::circos.trackPlotRegion(
      track.index = 1,
      bg.border = NA,
      panel.fun = function(x, y) {
        
        sector <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        
        circlize::circos.text(
          mean(xlim),
          ylim[1] + circlize::mm_y(2),
          sector,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = 0.65
        )
      }
    )
    
    # Percentage axes for condition sectors.
    
    for (sector in condition_sector_names) {
      
      circlize::circos.axis(
        h = "top",
        sector.index = sector,
        track.index = 2,
        major.at = seq(0, 100, by = 20),
        labels = TRUE,
        labels.cex = 0.45,
        major.tick.length = circlize::mm_y(1)
      )
    }
  }
  
  # Display the diagram in RStudio.
  
  draw_trophic_chord()
  
  # Save the diagram in PDF format.
  
  grDevices::pdf(
    file.path(output_dir, "trophic_mode_chord_diagram.pdf"),
    width = 12,
    height = 12
  )
  draw_trophic_chord()
  grDevices::dev.off()
  
  # Save the diagram in SVG format.
  
  grDevices::svg(
    file.path(output_dir, "trophic_mode_chord_diagram.svg"),
    width = 12,
    height = 12
  )
  draw_trophic_chord()
  grDevices::dev.off()
  
  # Save the diagram in PNG format.
  
  grDevices::png(
    file.path(output_dir, "trophic_mode_chord_diagram.png"),
    width = 3600,
    height = 3600,
    res = 300
  )
  draw_trophic_chord()
  grDevices::dev.off()
}


# 20. Guild analysis and visualization-----------

# Guilds are analysed independently from trophic modes.
# The chord diagram above represents trophic modes only.
#
# Guild categories are not restricted to a predefined list.
# All distinct guild labels found in FungalTraits are considered.
#
# Guild plots exclude:
# - Unclassified at genus level
# - Not in FungalTraits
# - No guild assigned
# - Multiple / unresolved
# - Unknown guild
#
# These categories remain in the exported CSV files.
# Relative abundance percentages use total community abundance
# as denominator; displayed guilds may therefore sum to less
# than 100%.


if ("guild_fg" %in% names(traits_raw)) {
  
  # Prepare genus and guild annotations.
  
  traits_guilds <- traits_raw %>%
    dplyr::mutate(
      Genus_clean = clean_genus_name(.data$Genus),
      Genus_key = stringr::str_to_lower(.data$Genus_clean),
      
      Guild_clean = stringr::str_squish(
        as.character(.data$guild_fg)
      ),
      
      Guild_clean = dplyr::na_if(
        .data$Guild_clean,
        ""
      ),
      
      Guild_clean = dplyr::case_when(
        is.na(.data$Guild_clean) |
          toupper(.data$Guild_clean) == "NULL" ~
          NA_character_,
        TRUE ~ .data$Guild_clean
      )
    )
  
  # Summarise guild annotations by genus.
  
  guilds_by_genus <- traits_guilds %>%
    dplyr::filter(
      !is.na(.data$Genus_key),
      .data$Genus_key != ""
    ) %>%
    dplyr::group_by(.data$Genus_key) %>%
    dplyr::summarise(
      Database_genus = dplyr::first(.data$Genus_clean),
      n_database_records = dplyr::n(),
      
      n_guilds = dplyr::n_distinct(
        .data$Guild_clean,
        na.rm = TRUE
      ),
      
      Guilds = {
        available_guilds <- sort(unique(
          .data$Guild_clean[!is.na(.data$Guild_clean)]
        ))
        
        if (length(available_guilds) == 0) {
          NA_character_
        } else {
          paste(available_guilds, collapse = "; ")
        }
      },
      
      .groups = "drop"
    )
  
  write_output_csv(
    guilds_by_genus %>%
      dplyr::select(
        Database_genus,
        n_database_records,
        n_guilds,
        Guilds
      ),
    "fungaltraits_guilds_by_genus.csv"
  )
  
  # Assign guild categories conservatively.
  # Genera with multiple distinct guilds remain unresolved.
  
  genus_guild_assignment <- observed_genera %>%
    dplyr::left_join(
      guilds_by_genus,
      by = "Genus_key"
    ) %>%
    dplyr::mutate(
      Assigned_guild = dplyr::case_when(
        .data$Genus == "Unclassified at genus level" ~
          "Unclassified at genus level",
        
        is.na(.data$n_database_records) ~
          "Not in FungalTraits",
        
        .data$n_guilds == 0 ~
          "No guild assigned",
        
        .data$n_guilds == 1 ~
          .data$Guilds,
        
        .data$n_guilds > 1 ~
          "Multiple / unresolved",
        
        TRUE ~ "Unknown guild"
      )
    )
  
  write_output_csv(
    genus_guild_assignment %>%
      dplyr::select(
        Genus,
        Assigned_guild,
        n_guilds,
        Guilds
      ),
    "observed_genus_guild_assignment.csv"
  )
  
  # Summarise relative abundance by guild and condition.
  
  guild_abundance_by_condition <- genus_mean_abundance %>%
    dplyr::left_join(
      genus_guild_assignment %>%
        dplyr::select(
          Genus,
          Assigned_guild
        ),
      by = "Genus"
    ) %>%
    dplyr::group_by(
      .data$Condition,
      .data$Assigned_guild
    ) %>%
    dplyr::summarise(
      Mean_relative_abundance = sum(
        .data$Mean_relative_abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Relative_abundance_pct =
        100 * .data$Mean_relative_abundance
    ) %>%
    dplyr::arrange(
      .data$Condition,
      .data$Assigned_guild
    )
  
  write_output_csv(
    guild_abundance_by_condition,
    "guild_abundance_by_condition.csv"
  )
  
  # Export guild abundance in wide format.
  
  guild_abundance_wide <- guild_abundance_by_condition %>%
    dplyr::select(
      Assigned_guild,
      Condition,
      Relative_abundance_pct
    ) %>%
    tidyr::pivot_wider(
      names_from = Condition,
      values_from = Relative_abundance_pct,
      values_fill = 0
    )
  
  write_output_csv(
    guild_abundance_wide,
    "guild_relative_abundance_pct_by_condition.csv"
  )
  
  # Exclude unassigned and unresolved categories from the figure only.
  
  excluded_guild_categories <- c(
    "Unclassified at genus level",
    "Not in FungalTraits",
    "No guild assigned",
    "Multiple / unresolved",
    "Unknown guild"
  )
  
  guild_plot_data <- guild_abundance_by_condition %>%
    dplyr::filter(
      !.data$Assigned_guild %in% excluded_guild_categories
    )
  
  if (nrow(guild_plot_data) > 0) {
    
    guild_plot <- ggplot2::ggplot(
      guild_plot_data,
      ggplot2::aes(
        x = .data$Assigned_guild,
        y = .data$Relative_abundance_pct,
        fill = .data$Condition
      )
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Relative Abundance by Fungal Guild",
        subtitle = "Guild annotations from FungalTraits",
        x = "Guild",
        y = "Relative abundance (%)",
        fill = "Condition"
      ) +
      ggplot2::theme_bw(base_size = 11)
    
    plot_height <- max(
      6,
      0.3 * dplyr::n_distinct(guild_plot_data$Assigned_guild) + 2
    )
    
    ggplot2::ggsave(
      file.path(output_dir, "guild_abundance_by_condition.png"),
      plot = guild_plot,
      width = 11,
      height = plot_height,
      dpi = 300
    )
    
    ggplot2::ggsave(
      file.path(output_dir, "guild_abundance_by_condition.pdf"),
      plot = guild_plot,
      width = 11,
      height = plot_height
    )
  }
}


# 21. Output-----------

# All output files are saved in 09_fungaltraits_results/.
#
# Main trophic mode outputs:
# - fungaltraits_trophic_modes_by_genus.csv
# - trophic_mode_assignment_review.csv
# - pending_trophic_mode_review.csv
# - observed_genus_trophic_mode_assignment.csv
# - fungaltraits_coverage_by_condition.csv
# - trophic_mode_abundance_by_condition.csv
# - trophic_mode_relative_abundance_pct_by_condition.csv
# - genera_by_trophic_mode.csv
# - genera_summary_by_trophic_mode.csv
# - trophic_mode_chord_diagram.pdf
# - trophic_mode_chord_diagram.svg
# - trophic_mode_chord_diagram.png
#
# Main guild outputs:
# - fungaltraits_guilds_by_genus.csv
# - observed_genus_guild_assignment.csv
# - guild_abundance_by_condition.csv
# - guild_relative_abundance_pct_by_condition.csv
# - guild_abundance_by_condition.pdf
# - guild_abundance_by_condition.png
#
# Individual sample identifiers are not exported.
# Abundance tables are summarised by genus, mode, guild or condition.

