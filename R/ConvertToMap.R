#' ConvertToMap
#'
#' Extracts physical and genetic position information for a single chromosome.
#'
#' @param data A data.frame containing physical and genetic position data, e.g., from a file like "https://maizegdb.org/map_text?id=1160762". If using a downloaded file directly, use ConvertToMap(read.delim("~/PATH/FILE.txt", skip = 1)).
#' @param physical_pos Character string specifying the column name for physical position (default: "Zm.B73.REFERENCE.NAM.5.0_start").
#' @param genetic_pos Character string specifying the column name for genetic position in centimorgans (default: "Coordinate").
#' @return A data.frame with columns `Mb` (physical position in megabases) and `cM` (genetic position in centimorgans).
#'
#' @importFrom dplyr %>% mutate select filter arrange
#' @importFrom rlang sym
#' @export
#'
ConvertToMap <- function(data, physical_pos = "Zm.B73.REFERENCE.NAM.5.0_start", genetic_pos = "Coordinate") {
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame")
  }
  if (!physical_pos %in% names(data)) {
    stop("Column specified by `physical_pos` (", physical_pos, ") not found in data")
  }
  if (!genetic_pos %in% names(data)) {
    stop("Column specified by `genetic_pos` (", genetic_pos, ") not found in data")
  }
  if (!is.numeric(data[[physical_pos]]) && !is.integer(data[[physical_pos]])) {
    warning("`physical_pos` column (", physical_pos, ") is not numeric; attempting to convert")
  }

  # Convert character strings to symbols for dplyr
  phys_col <- rlang::sym(physical_pos)
  gen_col <- rlang::sym(genetic_pos)

  # Process data
  data %>%
    mutate(Mb = !!phys_col / 1000000) %>%
    select(Mb, cM = !!gen_col) %>%
    filter(!is.na(Mb)) %>%
    mutate(Mb = as.numeric(Mb)) %>%
    arrange(Mb)
}
