#' insilico
#'
#' Generates in silico hybrid genotypes from inbred parents and writes them to a file.
#' File is written incrementally to disk as opposed to being stored in RAM. This reduces the speed of the operation but is required.
#'
#' The total memory needed is: Memory = (number of females)×(number of males)×(number of markers)×8bytes. This shows why writing to
#' disk is required as RAM requirements scales quickly (2570*1525*2490*8 ~ 72.7 Gb).
#'
#' @param geno_file_female File path to the female parent inbreds' genotype data.
#' @param geno_file_male File path to the male parent inbreds' genotype data.
#' @param output_file File path to save the hybrid genotypes (CSV).
#' @param coding_type Original genotype coding format (0,0.5,1/-1,0,1/0,1,2/0,2) with default "auto".
#' @return NULL (writes output to file instead of returning a data frame).
#'
#' @importFrom data.table fread
#' @export

insilico <- function(geno_file_female, geno_file_male, output_file = "insilico_hybrids.csv", coding_type = "auto") {
  # Function to detect input coding scheme
  detect_coding <- function(x) {
    unique_vals <- unique(na.omit(as.numeric(x)))
    if (all(unique_vals %in% c(0, 0.5, 1))) return("0_0.5_1")
    if (all(unique_vals %in% c(-1, 0, 1))) return("-1_0_1")
    if (all(unique_vals %in% c(0, 1, 2))) return("0_1_2")
    if (all(unique_vals %in% c(0, 2))) return("0_2")
    return("unknown")
  }

  # Load and process female genotypes
  markers_female <- as.data.frame(data.table::fread(geno_file_female))
  marker_info <- markers_female[, c("marker", "chrom", "position")]
  markers_female <- markers_female[, !names(markers_female) %in% c("marker", "chrom", "position"), drop = FALSE]
  markers_female <- t(markers_female)

  if (coding_type == "auto") coding_type <- detect_coding(markers_female)
  markers_female <- switch(coding_type,
                           "0_0.5_1" = markers_female * 2,
                           "-1_0_1" = markers_female + 1,
                           "0_1_2" = markers_female,
                           "0_2" = markers_female,
                           stop("Unknown or unsupported genotype coding scheme")
  )
  markers_female <- as.matrix(markers_female)
  class(markers_female) <- "integer"
  markers_female[!markers_female %in% c(0, 2) & !is.na(markers_female)] <- NA

  # Load and process male genotypes
  markers_male <- as.data.frame(data.table::fread(geno_file_male))
  markers_male <- markers_male[, !names(markers_male) %in% c("marker", "chrom", "position"), drop = FALSE]
  markers_male <- t(markers_male)

  if (coding_type == "auto") coding_type <- detect_coding(markers_male)
  markers_male <- switch(coding_type,
                         "0_0.5_1" = markers_male * 2,
                         "-1_0_1" = markers_male + 1,
                         "0_1_2" = markers_male,
                         "0_2" = markers_male,
                         stop("Unknown or unsupported genotype coding scheme")
  )
  markers_male <- as.matrix(markers_male)
  class(markers_male) <- "integer"
  markers_male[!markers_male %in% c(0, 2) & !is.na(markers_male)] <- NA

  # Check marker consistency
  if (!identical(colnames(markers_female), colnames(markers_male))) {
    stop("Marker names do not match between female and male datasets")
  }

  # Generate hybrid names
  n_female <- nrow(markers_female)
  n_male <- nrow(markers_male)
  hybrid_names <- outer(rownames(markers_female), rownames(markers_male), paste, sep = " X ")
  hybrid_names <- as.vector(hybrid_names)

  # Open file connection
  con <- file(output_file, "w")
  on.exit(close(con)) # Ensure file closes on exit

  # Write header
  header <- paste(c("marker", "chrom", "position", hybrid_names), collapse = ",")
  writeLines(header, con)

  # Process each marker and write row by row
  for (m in 1:ncol(markers_female)) {
    # Extract genotypes for this marker
    female_geno <- markers_female[, m]
    male_geno <- markers_male[, m]

    # Compute hybrid genotypes efficiently
    hybrid_geno <- (rep(female_geno, each = n_male) + rep(male_geno, times = n_female)) / 2

    # Format row
    marker_row <- c(
      marker_info$marker[m],
      marker_info$chrom[m],
      marker_info$position[m],
      format(hybrid_geno, scientific = FALSE, trim = TRUE)
    )
    row_text <- paste(marker_row, collapse = ",")

    # Write row
    writeLines(row_text, con)
  }

  invisible(NULL) # Return nothing; output is written to file
}
