#' model.terms
#'
#' Prepares the model terms required for StageWise functions, ensuring that the formatting aligns with the specific syntax
#' expectations of StageWise. This function allows users to specify fixed and random effects in a straightforward way,
#' and it handles the necessary formatting automatically.
#'
#' @param Fixed Character vector of the fixed effects for the model (e.g., c("env", "block")).
#' @param Random Character vector of the random effects for the model (e.g., c("genotype")).
#' @return A data frame with columns for the term names, their fixed/random designation, and factor status,
#'         formatted for StageWise compatibility.
#'
#' @export

model.terms <- function(Fixed, Random) {

  # Combine Fixed and Random terms into a single vector of terms
  terms <- c(Fixed, Random)

  # Generate flags indicating which terms are fixed (TRUE for fixed, FALSE for random)
  fixed_flags <- c(rep(TRUE, length(Fixed)), rep(FALSE, length(Random)))

  # Generate flags indicating factor status (FALSE for fixed, TRUE for random)
  # Assuming that Random effects are typically treated as factors
  factor_flags <- c(rep(FALSE, length(Fixed)), rep(TRUE, length(Random)))

  # Create the effects data frame with columns for name, fixed status, and factor status
  effects <- data.frame(name = terms, fixed = fixed_flags, factor = factor_flags)

  return(effects)
}
