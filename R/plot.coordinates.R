#' plot.coordinates
#'
#' This function can help to part in create a Design File as it takes a CSV Map of Plots laid out in field order (cell A1 is range:row coordinate 1,1) and generates a dataframe with the plot and its' corresponding range and row field coordinates.
#' @param Field.Map CSV file of plots laid out as the field is.
#' @return Dataframe of the plots and their range and row coordinates.
#' @export
plot.coordinates <- function(Field.Map) {
  # Read the CSV file
  data <- fread(Field.Map, header = FALSE)

  # Transpose Data
  data = as.matrix(data)

  # Determine the measure columns (assumes all columns are measure vars)
  measure_vars <- names(data)

  # Melt the data specifying the measure vars
  melted_data <- melt(data, measure.vars = measure_vars, variable.name = "Column", value.name = "Value", na.rm = TRUE)

  # Clarify Columns
  colnames(melted_data) = c("range", "row", "plot")

  # Remove "V" from row
  melted_data$row = gsub("V","",melted_data$row)

  # Remove Border and Empty Plots
  melted_data = melted_data[!melted_data$plot == "B" & !melted_data$plot == "",]

  # Formatting of Data
  melted_data$range = as.integer(melted_data$range)
  melted_data$row = as.integer(melted_data$row)
  melted_data$plot = as.character(melted_data$plot)

  return(melted_data)
}
