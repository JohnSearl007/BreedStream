#' ExpandData
#'
#' Makes a full grid so that spatial correction is possible. This is necessary as the ranges and rows of a field trial are not always perfectly filled with plots. In order for spatial correction to work with ASReml R, these ranges and rows must be complete so "filler" plots are created to satisfy this.
#' @param data Input data that has plots missing.
#' @return Input data with filler plots added where range/row were incomplete.
#' @export
ExpandData <- function(data) {
  max_row <- max(as.numeric(as.character(data$Row)), na.rm = TRUE)
  max_range <- max(as.numeric(as.character(data$Range)), na.rm = TRUE)

  dataexpanded <- expand.grid(Row = factor(1:max_row), Range = factor(1:max_range)) %>%
    left_join(data, by = c("Row", "Range"))

  return(dataexpanded)
}
