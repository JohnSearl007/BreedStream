#' format_markers_usef_add
#'
#' Formats genotype data coding so that SimpleMating::getUsefA() can be ran.
#'
#' @param geno_file File path to genotype data.
#' @param coding_type Original genotype coding format (0,0.5,1/-1,0,1/0,1,2/0,2) with the default of "auto" for automatic detection.
#' @return A matrix object with the genotype data formatted as 0,2,NA.
#'
#' @importFrom data.table fread
#' @export

format_markers_usef_add <- function(geno_file, coding_type = "auto") {
  # Read and initial processing
  markers <- as.data.frame(data.table::fread(geno_file))
  rownames(markers) <- markers$marker
  markers <- markers[, !names(markers) %in% c("marker", "chrom", "position")]
  markers <- t(markers)

  # Function to detect input coding scheme
  detect_coding <- function(x) {
    unique_vals <- unique(na.omit(as.numeric(x)))
    if (all(unique_vals %in% c(0, 0.5, 1))) return("0_0.5_1")
    if (all(unique_vals %in% c(-1, 0, 1))) return("-1_0_1")
    if (all(unique_vals %in% c(0, 1, 2))) return("0_1_2")
    if (all(unique_vals %in% c(0, 2))) return("0_2")
    return("unknown")
  }

  # Automatic detection if not specified
  if (coding_type == "auto") {
    coding_type <- detect_coding(markers)
  }

  # Convert to 0,2,NA based on detected or specified coding
  markers <- switch(coding_type,
                    "0_0.5_1" = markers * 2,
                    "-1_0_1" = markers + 1,
                    "0_1_2" = markers,
                    "0_2" = markers,
                    stop("Unknown or unsupported genotype coding scheme")
  )

  # Ensure only 0,2,NA values and convert to matrix with integer class
  markers <- as.matrix(markers)
  class(markers) <- "integer"
  markers[!markers %in% c(0, 2) & !is.na(markers)] <- NA

  return(markers)
}
