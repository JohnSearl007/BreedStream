#' nesting
#'
#' In order to handle nesting terms in StageWise you must create a new "composite" term. (Ex. Subblock 5 in Block 2 becomes one term "Subblock_5_Block_2") Therefore I need to handle performing the merging and removal of the old terms. This needs to be handled within the harvest.master function. nest.terms is a list of character vectors and nest.name is just a character vector.
#' @param df Origin data.
#' @param nest.terms List of character vectors to be nested.
#' @param nest.name Character vector or new nesting term name(s).
#' @return Dataframe that helps handle StageWise model syntax.
#' @export
nesting <- function(df, nest.terms, nest.name) {
  if (!is.null(nest.name)) {
    for (j in seq_along(nest.name)) {
      name <- nest.name[[j]]
      nest <- unlist(nest.terms[[j]])

      # Use unite() directly in a chain to avoid redundant selection
      df[[name]] <- df %>%
        unite(temp_column, all_of(nest), sep = "_") %>%
        pull(temp_column)
    }
    # Use select(-nest) once, instead of inside the loop
    df <- df %>% select(-all_of(unlist(nest.terms)))
  }
  return(df)
}
