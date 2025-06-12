#' coord.lookup
#'
#' This function assists in creating a Design File by taking a CSV map of plots arranged in field order (with cell A1 representing the range:row coordinate 1,1).
#' It generates a dataframe containing each plot with its corresponding range and row field coordinates.
#'
#' @param Field.Map CSV file with plots laid out in field order.
#' @return Dataframe containing plot numbers and their respective range and row coordinates.
#'
#' @importFrom data.table fread
#' @importFrom reshape2 melt
#' @export

coord.lookup <- function(Field.Map) {
  # Read the CSV file without headers as the data layout represents coordinates
  data <- data.table::fread(Field.Map, header = TRUE)

  # Convert the data to a matrix to facilitate melting
  data_matrix <- as.matrix(data)

  # Melt the matrix to long format with range and row coordinates
  melted_data <- reshape2::melt(data_matrix, varnames = c("range", "row"), value.name = "plot", na.rm = TRUE)

  # Remove any non-numeric characters from the row values if needed
  melted_data$row <- gsub("V", "", melted_data$row)

  # Filter out border plots ("B") and empty cells
  melted_data <- melted_data[!melted_data$plot %in% c("B", ""), ]

  # Convert columns to appropriate data types
  melted_data$range <- as.integer(melted_data$range)
  melted_data$row <- as.integer(melted_data$row)
  melted_data$plot <- as.character(melted_data$plot)

  return(melted_data)
}
