#' insilico_transposed
#'
#' Generates in silico hybrid genotypes from inbred parents and writes them to a file in a transposed format,
#' with hybrids as rows and markers as columns. The reason the in silico hybrids are being regenerated as opposed
#' to just using the existing in silico hybrids file is that the prediction stage is set up to use the in silico 
#' information in a transposed manner. Making the hybrids the rows was a choice as it is assumed that the user will 
#' have many more hybrids than markers. This makes the predictions much more efficient (long vs wide format). 
#' For large datasets, RAM limitations prevent simply taking the transpose of the existing in silico file. 
#' The file is written incrementally to disk to handle large datasets efficiently, avoiding high memory usage.
#'
#' The total memory requirement for the full dataset would be approximately:
#' Memory = (number of females) × (number of males) × (number of markers) × 8 bytes.
#' For example, 2570 females × 1525 males × 2490 markers × 8 bytes ~ 72.7 GB, making incremental writing essential.
#'
#' @param geno_file_female File path to the female parent inbreds' genotype data.
#' @param geno_file_male File path to the male parent inbreds' genotype data.
#' @param output_file File path to save the hybrid genotypes (CSV), defaults to "insilico_hybrids_transposed.csv".
#' @param coding_type Original genotype coding format (0,0.5,1/-1,0,1/0,1,2/0,2) with default "auto".
#' @return NULL (writes output to file instead of returning a data frame).
#'
#' @importFrom data.table fread
#' @export
insilico_transposed <- function(geno_file_female, geno_file_male, output_file = "insilico_hybrids_transposed.csv", coding_type = "auto") {
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
  
  # Get marker names for header
  marker_names <- marker_info$marker
  
  # Get parent names
  females <- rownames(markers_female)
  males <- rownames(markers_male)
  
  # Open file connection
  con <- file(output_file, "w")
  on.exit(close(con)) # Ensure file closes on exit
  
  # Write header: "hybrid" followed by marker names
  header <- paste(c("hybrid", marker_names), collapse = ",")
  writeLines(header, con)
  
  # Generate and write each hybrid row
  for (female in females) {
    female_geno <- markers_female[female, ]
    for (male in males) {
      hybrid_name <- paste(female, " X ", male, sep = "")
      male_geno <- markers_male[male, ]
      # Vectorized computation of hybrid genotypes
      hybrid_geno <- (female_geno + male_geno) / 2
      # Format to avoid scientific notation and ensure readability
      hybrid_geno_str <- format(hybrid_geno, scientific = FALSE, trim = TRUE)
      row_text <- paste(c(hybrid_name, hybrid_geno_str), collapse = ",")
      writeLines(row_text, con)
    }
  }
  
  invisible(NULL) # Return nothing; output is written to file
}
