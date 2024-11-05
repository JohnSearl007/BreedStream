#' model.terms_update
#'
#' The input file with the phenotypic information from the Harvest Master does not contain all the experimental design features and desired covariates for the model. The Model Term columns must be updated with the corresponding values, they cannot just remain NA, this is done with the Design File. The Design File must contain columns for; env, plot, range and row.
#' @param Fixed Character of the Fixed Effects for your model.
#' @param Random Character of the Random Effects for your model.
#' @return Dataframe that helps handle StageWise model syntax.
#' @export
model.terms_update <- function(df, Design.File, Fixed, Random) {
  # Ensure df is a data.frame for dplyr compatibility
  df <- as.data.frame(df)

  model.columns <- model.terms(Fixed, Random)
  df <- clean_names(df) %>%
    rename(plot = id_1)

  design <- fread(Design.File) %>%
    clean_names()

  # Use a map function for more concise class conversion
  convert_class <- function(df_col, design_col) {
    target_class <- class(df_col)
    design_col <- switch(target_class,
                         "character" = as.character(design_col),
                         "numeric" = as.numeric(design_col),
                         "factor" = as.factor(design_col),
                         "integer" = as.integer(design_col),
                         design_col)
    return(design_col)
  }

  # Vectorized conversion and joining, avoids loop overhead
  columns_to_match <- intersect(names(df), names(design))

  # Convert data.table to data.frame for safe column subsetting
  design <- as.data.frame(design)

  # Apply class conversion for matching columns
  design[columns_to_match] <- Map(convert_class, df[columns_to_match], design[columns_to_match])

  # Perform the join with dplyr::left_join()
  df <- left_join(df, design, by = c("env", "plot", "range", "row")) %>%
    select(-matches("\\.x$")) %>%
    rename_with(~ gsub("\\.y$", "", .x))

  return(df)
}
