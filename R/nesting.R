#' nesting
#'
#' This function creates composite "nested" terms by merging specified columns in `nest.terms` and assigning them a new
#' name specified by `nest.name`. This allows StageWise to handle nested terms in the format required for modeling.
#' For example, if Subblock 5 in Block 2 is nested, it can be represented as "Subblock_5_Block_2".
#'
#' @param df Data frame containing the original data.
#' @param nest.terms List of character vectors, where each vector contains column names to be nested together.
#' @param nest.name Character vector of new names for the nested terms.
#' @return A data frame with nested terms created and original columns removed.
#'
#' @importFrom dplyr %>% all_of select pull
#' @importFrom tidyr unite
#' @export

nesting <- function(df, nest.terms, nest.name) {

  # Check if nest.name is provided
  if (!is.null(nest.name)) {
    for (j in seq_along(nest.name)) {

      # Get the name for the new nested column and the columns to be nested
      name <- nest.name[[j]]
      nest <- unlist(nest.terms[[j]])

      # Create the nested term by concatenating specified columns with "_"
      df[[name]] <- df %>%
        tidyr::unite(temp_column, all_of(nest), sep = "_") %>%
        dplyr::pull(temp_column)
    }

    # Remove the original columns that were used to create nested terms
    df <- df %>% dplyr::select(-dplyr::all_of(unlist(nest.terms)))
  }

  return(df)
}
