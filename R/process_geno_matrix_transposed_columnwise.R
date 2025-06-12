#' Process Genotype Data into Centered Coefficient Matrices (Columnwise)
#'
#' Processes genotype data from a CSV file (hybrids as rows, markers as columns)
#' into centered coefficient matrices using a columnwise approach. Filters markers
#' based on minor allele count, centers genotypes around mean allele frequencies,
#' and optionally computes dominance coefficients. Results are stored in disk-based
#' `big.matrix` objects. Processing is done sequentially in batches to manage memory.
#'
#' @param filename Character string specifying the path to the CSV file containing
#'   genotype data (hybrids as rows, markers as columns). The file must include a
#'   header with "hybrid" as the first column followed by marker columns.
#' @param ploidy Integer specifying the ploidy level (e.g., 2 for diploid, 4 for tetraploid).
#'   Must be a positive even number. Defaults to 2.
#' @param dominance Logical indicating whether to compute dominance coefficients in
#'   addition to additive coefficients. If `TRUE`, a second `big.matrix` is created
#'   for dominance effects. Defaults to `TRUE`.
#' @param min_minor_allele Integer specifying the minimum number of hybrids carrying
#'   the minor allele required for a marker to be retained after filtering. Defaults to 1.
#' @param batch_size Integer specifying the number of markers to process in each batch.
#'   Adjust to balance memory usage and performance. Defaults to 35.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `coeff`: A `big.matrix.descriptor` pointing to the disk-based matrix of
#'       centered genotype coefficients (rows = hybrids, columns = markers).
#'     \item `coeff.D`: A `big.matrix.descriptor` for dominance coefficients (only
#'       included if `dominance = TRUE`).
#'     \item `id`: A character vector of hybrid IDs.
#'     \item `markers`: A character vector of kept marker names.
#'   }
#'
#' @details
#' The function processes the input CSV in two passes, using sequential batch processing:
#' \enumerate{
#'   \item \strong{First pass}: Reads the CSV in batches of markers to compute statistics
#'     (sum, count, minor allele counts) and filters markers based on `min_minor_allele`.
#'   \item \strong{Second pass}: Processes the CSV in batches of markers to compute
#'     centered genotype coefficients (and dominance coefficients if requested) for
#'     kept markers, storing results in disk-based `big.matrix` objects.
#' }
#' Genotypes are centered by subtracting `ploidy * p2`, where `p2` is the mean genotype
#' divided by ploidy. Missing values are imputed with the marker's mean genotype.
#' Dominance coefficients are computed as:
#' `-2 * choose(ploidy, 2) * p^2 + 2 * (ploidy - 1) * p * geno - geno * (geno - 1)`.
#'
#' @examples
#' \dontrun{
#' # Example with a diploid dataset
#' result <- process_geno_matrix_transposed_columnwise(
#'   filename = "insilico_hybrids_transposed.csv",
#'   ploidy = 2,
#'   dominance = TRUE,
#'   min_minor_allele = 5,
#'   batch_size = 50
#' )
#' print(result$id)
#' print(result$markers)
#' }
#'
#' @importFrom data.table fread setDTthreads
#' @importFrom bigmemory big.matrix describe
#' @importFrom foreach foreach %do%
#' @export
process_geno_matrix_transposed_columnwise <- function(filename, ploidy = 2, dominance = TRUE, min_minor_allele = 1, batch_size = 175) {
  # Ensure required packages are loaded
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")
  if (!requireNamespace("bigmemory", quietly = TRUE)) stop("bigmemory required")
  if (!requireNamespace("foreach", quietly = TRUE)) stop("foreach required")

  # Input validation
  if (!file.exists(filename)) stop("File does not exist: ", filename)
  if (!is.numeric(ploidy) || ploidy < 2 || ploidy %% 2 != 0) stop("Ploidy must be an even positive number")
  if (!is.logical(dominance)) stop("dominance must be TRUE or FALSE")
  if (!is.numeric(batch_size) || batch_size < 1 || batch_size %% 1 != 0) stop("batch_size must be a positive integer")
  if (!is.numeric(min_minor_allele) || min_minor_allele < 0) stop("min_minor_allele must be non-negative")

  # Precompute constant
  choose_ploidy_2 <- choose(ploidy, 2)

  # Read header to get marker names
  header <- data.table::fread(filename, nrows = 0)
  all_markers <- names(header)[-1]  # Exclude 'hybrid' column
  M <- length(all_markers)
  if (M == 0) stop("No markers found in CSV")

  # First pass: Compute statistics for filtering sequentially
  sum_geno <- numeric(M)
  count_geno <- integer(M)
  count_above_01 <- integer(M)
  count_below_09 <- integer(M)
  names(sum_geno) <- all_markers
  names(count_geno) <- all_markers
  names(count_above_01) <- all_markers
  names(count_below_09) <- all_markers

  num_batches <- ceiling(M / batch_size)
  batch_indices <- lapply(1:num_batches, function(b) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, M)
    start_idx:end_idx
  })

  first_pass_results <- foreach::foreach(
    b = 1:num_batches,
    .packages = "data.table",
    .combine = function(...) {
      results <- list(...)
      list(
        sum = do.call(c, lapply(results, function(x) x$sum)),
        count = do.call(c, lapply(results, function(x) x$count)),
        count_above_01 = do.call(c, lapply(results, function(x) x$count_above_01)),
        count_below_09 = do.call(c, lapply(results, function(x) x$count_below_09))
      )
    },
    .errorhandling = "stop"
  ) %do% {
    indices <- batch_indices[[b]]
    markers_batch <- all_markers[indices]
    geno_batch <- data.table::fread(filename, select = c("hybrid", markers_batch))

    # Convert to matrix for vectorized operations
    geno_matrix <- as.matrix(geno_batch[, -1, with = FALSE])
    colnames(geno_matrix) <- markers_batch

    # Vectorized computations
    batch_sum <- colSums(geno_matrix, na.rm = TRUE)
    batch_count <- colSums(!is.na(geno_matrix))
    pmat <- geno_matrix / ploidy
    batch_count_above_01 <- colSums(pmat > 0.1, na.rm = TRUE)
    batch_count_below_09 <- colSums(pmat < 0.9, na.rm = TRUE)

    list(
      sum = batch_sum,
      count = batch_count,
      count_above_01 = batch_count_above_01,
      count_below_09 = batch_count_below_09
    )
  }

  # Aggregate first pass results
  sum_geno[names(first_pass_results$sum)] <- first_pass_results$sum
  count_geno[names(first_pass_results$count)] <- first_pass_results$count
  count_above_01[names(first_pass_results$count_above_01)] <- first_pass_results$count_above_01
  count_below_09[names(first_pass_results$count_below_09)] <- first_pass_results$count_below_09

  # Compute allele frequency (p2) and filter markers
  mean_geno <- sum_geno / count_geno
  mean_geno[is.nan(mean_geno)] <- 0  # Handle division by zero
  p2 <- mean_geno / ploidy
  n_minor <- ifelse(p2 <= 0.5, count_above_01, count_below_09)
  keep_markers <- which(count_geno > 0 & n_minor >= min_minor_allele)
  filtered_markers <- all_markers[keep_markers]
  p2_filtered <- p2[keep_markers]
  M_filt <- length(filtered_markers)
  if (M_filt == 0) stop("No markers meet min_minor_allele threshold")

  # Read hybrid IDs
  id <- data.table::fread(filename, select = "hybrid", header = TRUE)$hybrid
  N <- length(id)

  # Initialize big.matrix for coeff and optionally coeff.D
  coeff <- bigmemory::big.matrix(nrow = N, ncol = M_filt, type = "double", init = 0,
                                 backingfile = "coeff.bin", descriptorfile = "coeff.desc")
  if (dominance) {
    coeff.D <- bigmemory::big.matrix(nrow = N, ncol = M_filt, type = "double", init = 0,
                                     backingfile = "coeff.D.bin", descriptorfile = "coeff.D.desc")
  }

  # Second pass: Compute coeff and coeff.D sequentially
  num_batches_filt <- ceiling(M_filt / batch_size)
  batch_indices_filt <- lapply(1:num_batches_filt, function(b) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, M_filt)
    list(indices = start_idx:end_idx, col_start = start_idx)
  })

  foreach::foreach(
    b = 1:num_batches_filt,
    .packages = c("data.table", "bigmemory"),
    .errorhandling = "stop"
  ) %do% {
    batch_info <- batch_indices_filt[[b]]
    indices <- batch_info$indices
    col_start <- batch_info$col_start
    markers_batch <- filtered_markers[indices]
    p2_batch <- p2_filtered[indices]
    geno_batch <- data.table::fread(filename, select = c("hybrid", markers_batch))

    # Convert to matrix for vectorized operations
    geno_matrix <- as.matrix(geno_batch[, -1, with = FALSE])
    colnames(geno_matrix) <- markers_batch
    na_mask <- is.na(geno_matrix)

    # Compute additive coefficients (coeff)
    p2_matrix <- matrix(p2_batch, nrow = N, ncol = length(p2_batch), byrow = TRUE)
    coeff_batch <- geno_matrix - ploidy * p2_matrix
    coeff_batch[na_mask] <- 0
    coeff[, col_start:(col_start + length(markers_batch) - 1)] <- coeff_batch

    # Compute dominance coefficients (coeff.D) if requested
    if (dominance) {
      geno_non_na <- geno_matrix
      geno_non_na[na_mask] <- ploidy * p2_matrix[na_mask]
      coeff.D_batch <- -2 * choose_ploidy_2 * p2_matrix^2 +
        2 * (ploidy - 1) * p2_matrix * geno_non_na -
        geno_non_na * (geno_non_na - 1)
      coeff.D_batch[na_mask] <- 0
      coeff.D[, col_start:(col_start + length(markers_batch) - 1)] <- coeff.D_batch
    }

    NULL
  }

  # Return results
  result <- list(coeff = bigmemory::describe(coeff), id = id, markers = filtered_markers)
  if (dominance) result$coeff.D <- bigmemory::describe(coeff.D)
  return(result)
}
