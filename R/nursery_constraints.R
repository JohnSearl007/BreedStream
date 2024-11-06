#' nursery_constraints
#'
#' This function allocates nursery rows to each biparental cross based on their percent contribution (`value`) to the next cycle.
#' Allocation stops once the specified number of `nursery.rows` has been reached.
#' This function is used for either the DH induction stage or the ear-to-row generational advancement stage, where generation of bi-parental crosses is not limiting.
#'
#' @param df Data frame from OMA with the proposed matings. Must include a `value` column indicating percent contribution of each cross.
#' @param nursery.rows Integer specifying the total nursery rows available for allocation.
#' @return A data frame indicating each cross and the number of nursery rows allocated to each.
#'
#' @export

nursery_constraints <- function(df, nursery.rows) {

  # Calculate the initial allocation of rows for each cross based on its contribution proportion
  df$nursery <- ceiling(df$value * nursery.rows)

  # Initialize cumulative sum of nursery rows allocated
  cumulative_rows <- base::cumsum(df$nursery)

  # Find the row index where cumulative sum first exceeds nursery.rows, truncate there
  cutoff_index <- base::min(which(cumulative_rows >= nursery.rows))

  # If the cumulative allocation meets or exceeds nursery.rows, retain rows up to cutoff
  df <- df[1:cutoff_index, ]

  return(df)
}
