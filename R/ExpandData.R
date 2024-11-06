#' ExpandData
#'
#' Expands input data to create a complete grid of rows and ranges, filling missing plots for spatial correction.
#'
#' @param data A data frame representing a field trial layout, with potential missing plots. Must contain columns named "Row" and "Range".
#' @return A data frame with additional filler plots for a complete grid of rows and ranges.
#'
#' @importFrom dplyr left_join %>%
#' @export

ExpandData <- function(data) {

  # Calculate the maximum row and range values from the input data
  max_row <- base::max(as.numeric(as.character(data$Row)), na.rm = TRUE)
  max_range <- base::max(as.numeric(as.character(data$Range)), na.rm = TRUE)

  # Generate a full grid of rows and ranges, then merge with input data to add missing plots
  dataexpanded <- base::expand.grid(Row = factor(1:max_row), Range = factor(1:max_range)) %>%
    dplyr::left_join(data, by = c("Row", "Range"))

  return(dataexpanded)  # Return the expanded data with filler plots added
}
