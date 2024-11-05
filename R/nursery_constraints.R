#' nursery_constraints
#'
#' oma and conversely oma_reduced returns the percent contribution for each cross to the next cycle. In a breeding program you are constrained by available resources in regards to the number of available rows. This returns the allocation of rows to each biparental cross at either the DH induction stage or the ear to row generational advancement stage. It is assumed that the generation of the bi-parental crosses is not a limiting factor.
#' @param df Dataframe from OMA with the proposed matings.
#' @return Crosses and how many nursery rows should be allocated to each.
#' @export
nursery_constraints = function(df, nursery.rows) {
  df = df
  df$nursery = ceiling(df$value*nursery.rows)

  for (i in 1:nrow(df)) {
    temp = df[1:i,]
    if (sum(temp$nursery) < nursery.rows) {
      i = i + 1
    } else {
      df = temp
      return(df)
    }
  }

}
