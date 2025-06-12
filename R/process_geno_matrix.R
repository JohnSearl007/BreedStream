#' Process Genotype Data into Centered Coefficient Matrices
#'
#' Reads genotype data from a CSV file (markers as rows, individuals as columns),
#' internally transposes it to individuals as rows and markers as columns using
#' disk-based operations, filters markers based on minor allele count, centers
#' genotypes around mean allele frequencies, and optionally computes dominance
#' coefficients. Results are stored in disk-based \code{big.matrix} objects.
#' Chunks are processed in parallel across multiple cores for efficiency.
#'
#' @param filename Character string specifying the path to the CSV file containing
#'   genotype data (markers as rows, individuals as columns). The file must include
#'   a header with at least one column for markers and one or more for genotypes.
#'   If \code{map = TRUE}, it must include columns for marker, chromosome, and position.
#' @param ploidy Integer specifying the ploidy level (e.g., 2 for diploid, 4 for tetraploid).
#'   Must be a positive even number.
#' @param map Logical indicating whether the CSV includes map information (columns for
#'   marker, chromosome, and position). If \code{TRUE}, the first three columns are
#'   assumed to be marker, chromosome, and position, with genotypes starting from the
#'   fourth column. If \code{FALSE}, the first column is markers, and genotypes start
#'   from the second column. Defaults to \code{FALSE}.
#' @param dominance Logical indicating whether to compute dominance coefficients in
#'   addition to additive coefficients. If \code{TRUE}, a second \code{big.matrix} is
#'   created for dominance effects. Defaults to \code{FALSE}.
#' @param chunk_size Integer specifying the number of individuals to process in each chunk
#'   after transposition. Adjust this to balance memory usage and performance. Defaults to 125000.
#' @param min_minor_allele Integer specifying the minimum number of individuals carrying
#'   the minor allele required for a marker to be retained after filtering. Defaults to 5.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{coeff}: A \code{big.matrix.descriptor} pointing to the disk-based
#'       matrix of centered genotype coefficients (rows = individuals, columns = markers).
#'     \item \code{coeff.D}: A \code{big.matrix.descriptor} for dominance coefficients
#'       (only included if \code{dominance = TRUE}).
#'     \item \code{id}: A character vector of individual IDs.
#'     \item \code{markers}: A character vector of kept marker names.
#'     \item \code{map}: A data.table with columns \code{marker}, \code{chr}, \code{pos}
#'       (only included if \code{map = TRUE}).
#'   }
#'
#' @details
#' The function internally transposes the input CSV (markers as rows, individuals as columns)
#' to a temporary file (individuals as rows, markers as columns) using a Bash script
#' with csvtool to avoid loading the entire dataset into memory. It assumes marker names
#' are unique to prevent errors during transposition or processing. It then processes the
#' transposed data in two passes, parallelizing chunks across multiple CPU cores:
#' \enumerate{
#'   \item \strong{First pass}: Reads the transposed CSV in chunks of individuals to compute
#'     per-marker statistics (sum, count, minor allele counts) and filters markers based
#'     on \code{min_minor_allele}.
#'   \item \strong{Second pass}: Reads the transposed CSV again to compute centered genotype
#'     coefficients (and dominance coefficients if requested) for kept markers, storing
#'     results in disk-based \code{big.matrix} objects.
#' }
#' Genotypes are centered by subtracting \code{ploidy * p2}, where \code{p2} is the mean
#' genotype divided by ploidy. Missing values are imputed with the marker's mean genotype.
#' Dominance coefficients are computed as:
#' \code{-2 * choose(ploidy, 2) * p^2 + 2 * (ploidy - 1) * p * geno - geno * (geno - 1)}.
#'
#' @examples
#' \dontrun{
#' # Example with a diploid dataset
#' result <- process_geno_matrix(
#'   filename = "genotypes.csv",
#'   ploidy = 2,
#'   map = TRUE,
#'   dominance = TRUE,
#'   chunk_size = 5000,
#'   min_minor_allele = 10
#' )
#' print(result$id)
#' print(result$map)
#' }
#'
#' @importFrom data.table fread setnames setDTthreads
#' @importFrom bigmemory big.matrix describe
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @export
process_geno_matrix <- function(filename, ploidy, map = FALSE, dominance = FALSE, chunk_size = 125000, min_minor_allele = 5) {
  # Ensure required packages are loaded
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package 'foreach' is required but not installed. Please install it.")
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' is required but not installed. Please install it.")
  }

  # Check for csvtool dependency
  if (system("csvtool --help", ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
    stop("csvtool is required but not installed. Please install csvkit.")
  }

  message("Starting process_geno_matrix with filename: ", filename, ", ploidy: ", ploidy,
          ", map: ", map, ", dominance: ", dominance, ", chunk_size: ", chunk_size,
          ", min_minor_allele: ", min_minor_allele)

  # Configure data.table threading
  data.table::setDTthreads(threads = parallel::detectCores() - 1)

  # Initialize parallel backend
  message("Initializing parallel backend with ", parallel::detectCores() - 1, " cores...")
  doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  on.exit(doParallel::stopImplicitCluster()) # Ensure parallel backend is stopped

  # Validate inputs
  message("Validating inputs...")
  if (!file.exists(filename)) stop("File does not exist: ", filename)
  if (!is.numeric(ploidy) || ploidy < 2 || ploidy %% 2 != 0) stop("Ploidy must be an even positive number")
  if (!is.logical(map)) stop("map must be TRUE or FALSE")
  if (!is.logical(dominance)) stop("dominance must be TRUE or FALSE")
  if (!is.numeric(chunk_size) || chunk_size < 1 || chunk_size %% 1 != 0) stop("chunk_size must be a positive integer")
  if (!is.numeric(min_minor_allele) || min_minor_allele < 0) stop("min_minor_allele must be non-negative")
  message("Input validation completed.")

  # Preprocess input CSV to ensure consistent encoding and line endings
  message("Preprocessing input CSV for consistent formatting...")
  temp_input <- tempfile(fileext = ".csv")
  system(paste("cat", shQuote(filename), "| iconv -t UTF-8 | tr -d '\r' >", shQuote(temp_input)))
  filename <- temp_input

  # Read header and determine columns
  message("Reading CSV header...")
  header <- data.table::fread(filename, nrows = 0)
  col_names <- names(header)
  if (map) {
    if (ncol(header) < 4) stop("With map = TRUE, CSV must have at least 4 columns (marker, chr, pos, genotypes)")
    map_cols <- 1:3
    geno_cols <- 4:ncol(header)
  } else {
    if (ncol(header) < 2) stop("With map = FALSE, CSV must have at least 2 columns (marker, genotypes)")
    map_cols <- 1
    geno_cols <- 2:ncol(header)
  }
  id <- col_names[geno_cols]
  N <- length(id)
  if (N == 0) stop("No genotype columns found in CSV")
  message("Header processed. Number of individuals (N): ", N)

  # Transpose CSV to temporary file using Bash with csvtool
  message("Transposing CSV to temporary file...")
  transposed_file <- tempfile(fileext = ".csv")
  temp_marker_file <- tempfile(fileext = ".txt")
  system(paste("tail -n +2", shQuote(filename), "| cut -d, -f1 >", shQuote(temp_marker_file)), intern = FALSE)
  bash_script <- c(
    "#!/bin/bash",
    "set -e # Exit on error",
    "input_csv=\"$1\"",
    "markers_file=\"$2\"",
    "output_csv=\"$3\"",
    "tmp_dir=\"$4\"",
    "",
    "# Ensure temporary directory exists",
    "mkdir -p \"$tmp_dir\" || { echo \"Error: Failed to create tmp_dir $tmp_dir\" >&2; exit 1; }",
    "",
    "# Extract column names (header) and markers",
    "header=$(head -n 1 \"$input_csv\") || { echo \"Error: Failed to read header from $input_csv\" >&2; exit 1; }",
    "IFS=',' read -r -a header_array <<< \"$header\"",
    "markers=$(cat \"$markers_file\" | tr '\\n' ',' | sed 's/,$//') || { echo \"Error: Failed to read markers from $markers_file\" >&2; exit 1; }",
    "",
    "# Determine genotype column indices",
    "geno_start=$([ $(csvtool col 2 \"$input_csv\" >/dev/null 2>&1; echo $?) -eq 0 ] && echo 2 || echo 4)",
    "",
    "# Get number of columns",
    "num_cols=$(csvtool width \"$input_csv\") || { echo \"Error: Failed to determine number of columns\" >&2; exit 1; }",
    "",
    "# Create temporary files for each individual (genotype columns)",
    "for ((i=$geno_start; i<=$num_cols; i++)); do",
    "  col_idx=$i",
    "  tmp_file=\"$tmp_dir/transpose_tmp.$((i-geno_start+1))\"",
    "  if [ \"$geno_start\" -eq 4 ]; then",
    "    csvtool col $col_idx \"$input_csv\" | tail -n +4 > \"$tmp_file\" || { echo \"Error: Failed to extract column $col_idx to $tmp_file\" >&2; exit 1; }",
    "  else",
    "    csvtool col $col_idx \"$input_csv\" | tail -n +2 > \"$tmp_file\" || { echo \"Error: Failed to extract column $col_idx to $tmp_file\" >&2; exit 1; }",
    "  fi",
    "done",
    "",
    "# Create output CSV header",
    "echo \"individual,$markers\" > \"$output_csv\" || { echo \"Error: Failed to write header to $output_csv\" >&2; exit 1; }",
    "",
    "# Process each individual",
    "num_individuals=$((num_cols - geno_start + 1))",
    "for ((i=1; i<=$num_individuals; i++)); do",
    "  individual_idx=$i",
    "  individual_name=\"${header_array[$((i + geno_start - 2))]}\"",
    "  tmp_file=\"$tmp_dir/transpose_tmp.$individual_idx\"",
    "  # Read genotype values",
    "  geno_values=$(cat \"$tmp_file\" | tr '\\n' ',' | sed 's/,$//') || { echo \"Error: Failed to read genotypes from $tmp_file\" >&2; exit 1; }",
    "  # Write row to output CSV",
    "  echo \"$individual_name,$geno_values\" >> \"$output_csv\"",
    "done",
    "",
    "# Clean up temporary files",
    "for ((i=$geno_start; i<=$num_cols; i++)); do",
    "  rm \"$tmp_dir/transpose_tmp.$((i-geno_start+1))\" || { echo \"Error: Failed to remove $tmp_dir/transpose_tmp.$((i-geno_start+1))\" >&2; exit 1; }",
    "done"
  )
  bash_script_file <- tempfile(fileext = ".sh")
  writeLines(bash_script, bash_script_file)
  system(paste("chmod +x", shQuote(bash_script_file)), intern = FALSE)
  system(paste(shQuote(bash_script_file), shQuote(filename), shQuote(temp_marker_file), shQuote(transposed_file), shQuote(dirname(transposed_file))), intern = TRUE)
  message("Checking transposed CSV: ", transposed_file)
  if (file.exists(transposed_file)) {
    file_info <- file.info(transposed_file)
    message("Transposed file size: ", file_info$size, " bytes")
    if (map) {
      message("Removing 'chrom' and 'position' rows from transposed CSV...")
      temp_output_file <- tempfile(fileext = ".csv")
      system(paste("awk 'NR!=2 && NR!=3' ", shQuote(transposed_file), " > ", shQuote(temp_output_file)), intern = FALSE)
      if (file.exists(temp_output_file) && file.info(temp_output_file)$size > 0) {
        system(paste("mv", shQuote(temp_output_file), shQuote(transposed_file)), intern = FALSE)
        message("Transposed file updated, rows 2 and 3 removed")
      } else {
        stop("Failed to create or write to temporary output CSV: ", temp_output_file)
      }
    } else {
      message("map = FALSE, skipping removal of 'chrom' and 'position' rows")
    }
    file_info_updated <- file.info(transposed_file)
    message("Transposed file size after processing: ", file_info_updated$size, " bytes")
  } else {
    stop("Transposed CSV does not exist: ", transposed_file)
  }
  message("Debug: First 5 rows and 5 columns of transposed CSV:")
  system(paste("head -n 5", shQuote(transposed_file), "| csvtool -u TAB col 1-5 -"))
  file.remove(temp_marker_file)

  # Get marker names from transposed CSV header
  message("Reading marker names from transposed CSV header...")
  transposed_header <- data.table::fread(transposed_file, nrows = 0)
  if (!"individual" %in% names(transposed_header)) {
    stop("Transposed CSV does not have 'individual' as the first column. Header: ", paste(names(transposed_header), collapse = ", "))
  }
  all_markers <- names(transposed_header)[-1]  # Exclude "individual"
  M <- length(all_markers)
  if (M == 0) stop("No markers found in transposed CSV")
  message("Marker names read. Total number of markers (M): ", M)

  # First pass: Compute per-marker statistics in chunks of individuals
  message("Starting first pass: computing per-marker statistics...")
  sum_geno <- numeric(M)
  count_geno <- numeric(M)
  count_above_01 <- numeric(M)
  count_below_09 <- numeric(M)
  keep_markers <- logical(M)
  names(sum_geno) <- names(count_geno) <- names(count_above_01) <- names(count_below_09) <- names(keep_markers) <- all_markers
  num_chunks <- ceiling(N / chunk_size)
  message("Number of chunks for first pass: ", num_chunks)

  # Parallelize first pass
  first_pass_results <- foreach::foreach(i = 1:num_chunks, .packages = c("data.table"), .combine = c, .errorhandling = "stop") %dopar% {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, N)
    chunk_size_actual <- end_idx - start_idx + 1

    # Read chunk with numeric genotypes
    col_types <- c("character", rep("double", length(all_markers)))
    names(col_types) <- c("individual", all_markers)
    tryCatch({
      geno_chunk <- data.table::fread(transposed_file,
                                      skip = start_idx,
                                      nrows = chunk_size_actual,
                                      header = TRUE,
                                      colClasses = col_types)
      data.table::setnames(geno_chunk, c("individual", all_markers))
    }, error = function(e) {
      stop("Error reading chunk ", i, " from transposed CSV: ", conditionMessage(e), "\nCheck transposed CSV at: ", transposed_file)
    })

    # Check for markers with all NA values
    geno_cols_dt <- all_markers
    na_cols <- geno_cols_dt[sapply(geno_chunk[, .SD, .SDcols = geno_cols_dt], function(x) all(is.na(x)))]
    if (length(na_cols) > 0) {
      warning("Chunk ", i, ": Markers ", paste(na_cols, collapse = ", "), " contain non-numeric data; treating as NA")
    }

    # Compute sums, counts, and minor allele counts
    chunk_sum <- colSums(geno_chunk[, .SD, .SDcols = geno_cols_dt], na.rm = TRUE)
    chunk_count <- colSums(!is.na(geno_chunk[, .SD, .SDcols = geno_cols_dt]))
    chunk_count_above_01 <- colSums(geno_chunk[, .SD, .SDcols = geno_cols_dt] / ploidy > 0.1, na.rm = TRUE)
    chunk_count_below_09 <- colSums(geno_chunk[, .SD, .SDcols = geno_cols_dt] / ploidy < 0.9, na.rm = TRUE)

    # Garbage collection
    gc()

    # Return results as a named list
    list(list(sum = chunk_sum, count = chunk_count, count_above_01 = chunk_count_above_01, count_below_09 = chunk_count_below_09))
  }

  # Aggregate first pass results
  message("Aggregating first pass results...")
  for (res in first_pass_results) {
    sum_geno <- sum_geno + res$sum
    count_geno <- count_geno + res$count
    count_above_01 <- count_above_01 + res$count_above_01
    count_below_09 <- count_below_09 + res$count_below_09
  }

  # Compute allele frequencies and filter markers
  message("Computing allele frequencies and filtering markers...")
  mean_geno <- sum_geno / count_geno
  mean_geno[is.nan(mean_geno)] <- 0
  p2 <- mean_geno / ploidy
  n_minor <- ifelse(p2 <= 0.5, count_above_01, count_below_09)
  keep_markers <- n_minor >= min_minor_allele & count_geno > 0
  all_markers <- all_markers[keep_markers]
  sum_geno <- sum_geno[keep_markers]
  count_geno <- count_geno[keep_markers]
  mean_geno <- mean_geno[keep_markers]
  p2 <- p2[keep_markers]
  M_kept <- length(all_markers)
  if (M_kept == 0) stop("No markers meet min_minor_allele threshold")
  message("Number of markers after filtering: ", M_kept)

  # Initialize disk-based matrices
  message("Initializing disk-based matrices...")
  coeff <- bigmemory::big.matrix(nrow = N, ncol = M_kept, type = "double", init = 0,
                                 backingfile = "coeff.bin", descriptorfile = "coeff.desc")
  message("Additive coefficient matrix initialized. Dimensions: ", N, " rows x ", M_kept, " cols")
  if (dominance) {
    coeff.D <- bigmemory::big.matrix(nrow = N, ncol = M_kept, type = "double", init = 0,
                                     backingfile = "coeff.D.bin", descriptorfile = "coeff.D.desc")
    message("Dominance coefficient matrix initialized. Dimensions: ", N, " rows x ", M_kept, " cols")
  }

  # Second pass: Compute coefficients in chunks of individuals
  message("Starting second pass: filling coefficient matrices...")
  kept_marker_cols <- which(keep_markers)
  num_chunks <- ceiling(N / chunk_size)
  message("Number of chunks for second pass: ", num_chunks)

  # Parallelize second pass
  foreach::foreach(i = 1:num_chunks, .packages = c("data.table", "bigmemory")) %dopar% {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, N)
    chunk_size_actual <- end_idx - start_idx + 1
    message("Chunk ", i, ": Processing rows ", start_idx, " to ", end_idx, " (", chunk_size_actual, " rows)")

    # Calculate column indices to select
    col_indices <- c(1, 1 + kept_marker_cols)

    # Read chunk with header=FALSE, skipping header line appropriately
    geno_chunk <- data.table::fread(transposed_file,
                                    skip = 1 + (i - 1) * chunk_size,
                                    nrows = chunk_size_actual,
                                    header = FALSE)
    geno_chunk <- geno_chunk[, col_indices, with = FALSE]
    data.table::setnames(geno_chunk, c("individual", all_markers))
    geno_chunk[, (all_markers) := lapply(.SD, as.numeric), .SDcols = all_markers]

    message("Chunk ", i, ": fread() completed")

    # Impute missing values and compute centered coefficients
    geno_cols_dt <- all_markers
    centered_chunk <- geno_chunk[, .SD, .SDcols = geno_cols_dt]
    centered_chunk[, (geno_cols_dt) := lapply(.SD, function(x, marker) {
      ifelse(is.na(x), mean_geno[match(marker, all_markers)], x)
    }, marker = geno_cols_dt), .SDcols = geno_cols_dt]

    message("Chunk ", i, ": Imputation completed for centered_chunk")

    centered_chunk[, (geno_cols_dt) := lapply(seq_along(.SD), function(j) {
      .SD[[j]] - ploidy * p2[geno_cols_dt[j]]
    }), .SDcols = geno_cols_dt]

    message("Chunk ", i, ": Centering completed")

    # Compute dominance coefficients
    if (dominance) {
      dom_chunk <- geno_chunk[, .SD, .SDcols = geno_cols_dt]
      dom_chunk[, (geno_cols_dt) := lapply(seq_along(.SD), function(j) {
        x <- .SD[[j]]
        marker <- geno_cols_dt[j]
        geno <- ifelse(is.na(x), mean_geno[match(marker, all_markers)], x)
        p <- p2[match(marker, all_markers)]
        result <- -2 * choose(ploidy, 2) * p^2 + 2 * (ploidy - 1) * p * geno - geno * (geno - 1)
        if (length(result) != length(x)) {
          stop("Length mismatch in dominance for marker ", marker, ": expected ", length(x), ", got ", length(result))
        }
        result
      }), .SDcols = geno_cols_dt]
      message("Chunk ", i, ": Dominance coefficients computed")
    }

    # Assign to big.matrix
    coeff[start_idx:end_idx, ] <- as.matrix(centered_chunk)
    if (dominance) coeff.D[start_idx:end_idx, ] <- as.matrix(dom_chunk)

    # Garbage collection
    gc()

    # Log completion
    message("Chunk ", i, ": Completed processing ", chunk_size_actual, " individuals.")

    NULL
  }

  # Clean up temporary files
  message("Cleaning up temporary files...")
  # file.remove(transposed_file) # Keep to prevent wait in future iterations
  file.remove(temp_input)

  # Return descriptors for disk-based matrices
  message("Returning results...")
  result <- list(coeff = bigmemory::describe(coeff), id = id, markers = all_markers)
  if (dominance) result$coeff.D <- bigmemory::describe(coeff.D)
  message("process_geno_matrix completed successfully.")
  return(result)
}
