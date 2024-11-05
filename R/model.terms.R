#' model.terms
#'
#' StageWise requires the model terms be passed to the various functions in specific formatting. This includes the input data files. This function allows the user to enter the model in a more straight forward manner and have the formatting handled for you.
#' @param Fixed Character of the Fixed Effects for your model.
#' @param Random Character of the Random Effects for your model.
#' @return Dataframe that helps handle StageWise model syntax.
#' @export
model.terms <- function(Fixed, Random) {
  terms <- c(Fixed, Random)
  fixed_flags <- c(rep(TRUE, length(Fixed)), rep(FALSE, length(Random)))
  factor_flags <- c(rep(FALSE, length(Fixed)), rep(TRUE, length(Random)))

  effects <- data.frame(name = terms, fixed = fixed_flags, factor = factor_flags)
  return(effects)
}
