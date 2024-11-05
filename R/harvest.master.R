#' harvest.master
#'
#' This function will take the raw output file directly off of a H2 Harvest Master system with the Mirus Software and apply proper formatting so that the output is ready for StageWise. The user must provide a Design File that has columns "plot", "range", "row" that ties the Harvest Master output to that of the respective experimental design features for said plot such as Block. User must specify the location (environment).
#' @param Location Location name for the data.
#' @param Design.File Filename for the Design.File.
#' @param Combine.File Filename for the Combine.File.
#' @param Fixed Character for Fixed Terms, "id_1" is the default.
#' @param Random Character for Random Terms.
#' @param Moisture Numeric for the standard moisture content.
#' @param Bushel Numeric for the standard weight of a bushel.
#' @param Plot.Length Numeric for the length of the plot in feet.
#' @param Row.Spacing Numeric for the spacing of rows in inches.
#' @param Plot.Rows Numeric for the number of rows in the plot.
#' @param nest.terms List of character vectors to be nested.
#' @param nest.name Character vector or new nesting term name(s).
#' @return Dataframe that is properly formated and has had the yield calculated.
#' @export
harvest.master <- function(Location, Design.File, Combine.File, Fixed = c("id_1"), Random, Moisture = 15.5, Bushel = 56, Plot.Length = 22.5,
                           Row.Spacing = 30, Plot.Rows = 2, nest.terms = NULL, nest.name = NULL) {
  H2.columns <- c("Date/Time", "Range", "Row", "Id 1", "Weight (lb)", "Moisture (%)", "Test Weight (lb/bu)", "Quick Note", "Harvest Sequence")
  df <- fread(Combine.File) %>%
    mutate(env = Location)

  # Add missing columns, if any, using `data.table` for speed
  missing.columns <- setdiff(H2.columns, colnames(df))
  if (length(missing.columns) > 0) {
    df[, (missing.columns) := NA]
  }

  # Select the required columns
  H2.columns <- append(H2.columns, "env")
  df <- df[, ..H2.columns]

  # Ensure column classes match
  H2.classes <- c("character", "integer", "integer", "character", "numeric", "numeric", "numeric", "character", "integer")
  setDT(df)  # Ensure df is a data.table
  for (i in seq_along(H2.classes)) {
    set(df, j = i, value = as(df[[i]], H2.classes[i]))
  }

  # Update model terms (convert to data.frame for dplyr compatibility)
  df <- model.terms_update(as.data.frame(df), Design.File, Fixed, Random)

  # Convert back to data.table to allow the use of `:=`
  setDT(df)

  # Expand data if necessary
  if (any(is.na(df$row)) || any(is.na(df$range))) {
    df <- ExpandData(df)
  }

  # Yield calculation (vectorized)
  df[, Adj_Weight := (weight_lb - (weight_lb * (moisture_percent / 100))) / ((100 - Moisture) / 100)]
  df[, Yield := (Adj_Weight / Bushel) / (1 / (43560 / (Plot.Length * (Row.Spacing / 12) * Plot.Rows)))]

  # Handle nesting terms
  df <- nesting(df, nest.terms = nest.terms, nest.name = nest.name)

  # Output the cleaned data
  return(clean_names(df))
}
