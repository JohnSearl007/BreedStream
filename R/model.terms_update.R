#' model.terms_update
#'
#' Updates the model term columns in the input data frame by joining it with a design file that contains experimental design
#' features and covariates. This is necessary as the raw phenotypic data file from Harvest Master may lack these
#' experimental design features (e.g., range, row, etc.) which are required for model building.
#'
#' The design file must include columns: "env", "plot", "range", and "row".
#'
#' @param df Data frame containing the initial phenotypic information from the Harvest Master system.
#' @param Design.File Filename (or path) of the design file containing the design features and covariates.
#' @param Fixed Character vector specifying the fixed effects for the model.
#' @param Random Character vector specifying the random effects for the model.
#' @return A data frame with updated model terms, including design features and covariates for StageWise compatibility.
#'
#' @importFrom data.table fread setDT
#' @importFrom dplyr %>% intersect rename
#' @importFrom janitor clean_names
#' @export

model.terms_update <- function(df, Design.File, Fixed, Random) {

  # Ensure `df` is a data.frame for compatibility with `dplyr` functions
  df <- as.data.frame(df)

  # Generate the columns required for the model terms
  model.columns <- model.terms(Fixed, Random)

  # Rename plot identifier column to standardize with design file columns
  df <- janitor::clean_names(df) %>%
    dplyr::rename(plot = id_1)

  # Read in the design file and clean column names for consistency
  design <- data.table::fread(Design.File) %>%
    janitor::clean_names()

  # Convert `design` columns to match `df` classes where applicable
  convert_class <- function(df_col, design_col) {
    target_class <- class(df_col)
    design_col <- switch(target_class,
                         "character" = as.character(design_col),
                         "numeric" = as.numeric(design_col),
                         "integer" = as.integer(design_col),
                         "factor" = as.factor(design_col),
                         design_col)
    return(design_col)
  }

  columns_to_match <- dplyr::intersect(names(df), names(design))
  for (col in columns_to_match) {
    design[[col]] <- convert_class(df[[col]], design[[col]])
  }

  # Convert to data.table and perform join using 'on'
  data.table::setDT(df)
  data.table::setDT(design)
  df <- df[design, on = .(env, plot, range, row), nomatch = 0]

  return(janitor::clean_names(df))
}
