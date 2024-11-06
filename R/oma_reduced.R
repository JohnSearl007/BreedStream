#' oma_reduced
#'
#' This function performs an optimized mating allocation (OMA) with reduced computational demands by selecting a top subset of parental candidates via optimal contribution selection (OCS) before running OMA. This is useful when a high number of parent combinations would otherwise exceed memory constraints.
#'
#' @param dF Numeric vector of length 2 indicating inbreeding rate constraints for OMA. Defaults to c(0.01, 0.01).
#' @param geno.file File path for genotype data to be used in OMA.
#' @param kinship.file File path for kinship matrix to be used in OMA.
#' @param ploidy Integer indicating ploidy level of the species. Defaults to 2.
#' @param selection.intensity Numeric value representing the proportion of top parents to select based on OCS. Defaults to 0.10 (top 10%).
#' @param solver Character string specifying the solver COMA should utilize, e.g., "ECOS".
#' @param dF.adapt List specifying adaptive parameters for `dF` if the initial value does not yield a solution (e.g., list(step = 0.005, max = 1)).
#' @param base Character indicating if the base population consists of random matings ("RM") or inbred lines.
#' @return A data frame containing the OMA output with optimal mating combinations.
#'
#' @importFrom COMA ocs oma
#' @export

oma_reduced <- function(dF = c(0.01, 0.01), geno.file, kinship.file, ploidy = 2, selection.intensity = 0.10,
                        solver = "ECOS", dF.adapt = list(step = 0.005, max = 1), base = "RM") {

  # Input validation
  stopifnot(length(dF) == 2L)  # Ensure dF has two elements
  stopifnot(dF[1] <= dF[2])    # Ensure lower bound of dF is not greater than the upper bound

  # Stage 1: Load data and calculate initial OCS allocation for parent selection
  ans1 <- read_data_optimized(geno.file = geno.file,
                              kinship.file = kinship.file,
                              ploidy = ploidy,
                              matings = "all",
                              standardize = TRUE)

  # Calculate maximum number of parents based on selection intensity
  max.parent <- ceiling(length(ans1$parents$id) * selection.intensity)

  # Stage 2: Run OCS to select the top parents based on contribution constraints
  ans2 <- COMA::ocs(parents = data.frame(ans1$parents, min = 0, max = 1 / max.parent),
                    ploidy = ploidy,
                    K = ans1$K,
                    dF = dF[2],
                    dF.adapt = dF.adapt,
                    solver = solver,
                    base = base)

  # Check for solution in OCS results
  if (nrow(ans2$oc) == 0) {
    stop("No solution possible in OCS.")
  }

  # Select the top parents based on their OCS values
  sel1 <- ans2$oc$id[order(ans2$oc$value, decreasing = TRUE)]
  if (length(sel1) > max.parent) {
    sel1 <- sel1[1:max.parent]
  }

  # Stage 3: Reload data, now using only the selected parents for reduced mating combinations
  ans1 <- read_data_optimized(geno.file = geno.file,
                              kinship.file = kinship.file,
                              ploidy = ploidy,
                              matings = sel1,
                              standardize = TRUE)

  # Stage 4: Run OMA with the selected parents and mating combinations
  ans3 <- COMA::oma(parents = data.frame(id = sel1, min = 0, max = 1),
                    matings = data.frame(ans1$matings, min = 0, max = 1),
                    ploidy = ploidy,
                    K = ans1$K,
                    dF = dF,
                    dF.adapt = dF.adapt,
                    solver = solver,
                    base = base)

  # Check for solution in OMA results
  if (nrow(ans3$om) == 0) {
    stop("No solution possible in OMA.")
  }

  # Return the results from OMA
  return(ans3)
}
