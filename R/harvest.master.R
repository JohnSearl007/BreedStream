#' harvest.master
#'
#' Processes raw output from an H2 Harvest Master system with Mirus Software, applying formatting and calculations to
#' prepare the data for use in StageWise. Requires a design file with columns "plot", "range", and "row" that links
#' the Harvest Master data to experimental design features such as blocks. The location (environment) of the data
#' collection must be specified.
#'
#' @param Location Character; the name of the location (environment) for the data.
#' @param Design.File Filename (or path) for the design file containing experimental design details.
#' @param Combine.File Filename (or path) for the raw Harvest Master output file.
#' @param Fixed Character vector for fixed terms; "id_1" is the default.
#' @param Random Character vector for random terms.
#' @param Moisture Numeric; the standard moisture content to adjust weights (default is 15.5).
#' @param Bushel Numeric; standard weight of a bushel (default is 56).
#' @param Plot.Length Numeric; the length of each plot in feet (default is 22.5).
#' @param Row.Spacing Numeric; spacing of rows in inches (default is 30).
#' @param Plot.Rows Numeric; the number of rows in each plot (default is 2).
#' @param nest.terms List of character vectors specifying terms to be nested.
#' @param nest.name Character vector specifying new name(s) for nesting term(s).
#' @return A data frame formatted for StageWise, with calculated yield values.
#'
#' @importFrom data.table fread set setDT
#' @importFrom dplyr mutate setdiff %>%
#' @importFrom janitor clean_names
#' @export

harvest.master <- function(Location, Design.File, Combine.File, Fixed = c("id_1"), Random, Moisture = 15.5,
                           Bushel = 56, Plot.Length = 22.5, Row.Spacing = 30, Plot.Rows = 2,
                           nest.terms = NULL, nest.name = NULL) {

  # Define expected columns from Harvest Master system
  H2.columns <- c("Date/Time", "Range", "Row", "Id 1", "Weight (lb)", "Moisture (%)",
                  "Test Weight (lb/bu)", "Quick Note", "Harvest Sequence")

  # Load Harvest Master data, adding the environment (location) column
  df <- data.table::fread(Combine.File) %>%
    dplyr::mutate(env = Location)

  # Identify and add any missing columns to match Harvest Master format
  missing.columns <- dplyr::setdiff(H2.columns, colnames(df))
  if (length(missing.columns) > 0) {
    df[, (missing.columns) := NA]
  }

  # Ensure the selected columns are present and ordered as expected
  H2.columns <- append(H2.columns, "env")
  df <- df[, ..H2.columns]

  # Define target column data types and apply them
  H2.classes <- c("character", "integer", "integer", "character", "numeric",
                  "numeric", "numeric", "character", "integer")
  data.table::setDT(df)  # Ensure df is a data.table object for efficient operations
  for (i in seq_along(H2.classes)) {
    data.table::set(df, j = i, value = as(df[[i]], H2.classes[i]))
  }

  # Update model terms using the design file and specified Fixed and Random effects
  df <- model.terms_update(as.data.frame(df), Design.File, Fixed, Random)

  # Reconvert to data.table for further processing
  data.table::setDT(df)

  # Expand data to ensure spatial completeness if necessary
  if (any(is.na(df$row)) || any(is.na(df$range))) {
    df <- ExpandData(df)
  }

  # Calculate adjusted weight based on moisture content, then compute yield
  df[, Adj_Weight := (weight_lb - (weight_lb * (moisture_percent / 100))) / ((100 - Moisture) / 100)]
  df[, Yield := (Adj_Weight / Bushel) / (1 / (43560 / (Plot.Length * (Row.Spacing / 12) * Plot.Rows)))]

  # Apply nesting terms if provided
  df <- nesting(df, nest.terms = nest.terms, nest.name = nest.name)

  # Clean up column names for final output
  return(janitor::clean_names(df))
}
