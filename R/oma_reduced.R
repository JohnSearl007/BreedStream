#' oma_reduced
#'
#' When the number of parents results in too many pairwise mating combinations to evaluate (due to RAM), this function will first perform OCS and identify the top parental candidates. Then using the OCS output OMA will be conducted in order to find the best matings.
#' @param df Inbreeding rate constraints.
#' @param geno.file Genotype information for COMA.
#' @param kinship.file Kinship matrix for COMA.
#' @param ploidy Ploidy of species.
#' @param selection.intensity The top percent of parents to consider using for large sets of germplasm.
#' @param solver The solver COMA should utilize.
#' @param dF.adapt If the specified dF does not find a solution how should COMA proceed.
#' @param base Are the parents in question inbreds and COMA should use random mating.
#' @return OMA output.
#' @export
oma_reduced = function(dF = c(0.01,0.01), geno.file, kinship.file, ploidy = 2, selection.intensity = 0.10, solver="ECOS", dF.adapt = list(step = 0.005, max = 1), base = "RM") {

  stopifnot(length(dF) == 2L)
  stopifnot(dF[1] <= dF[2])

  ans1 = read_data_optimized(geno.file = geno.file,
                             kinship.file = kinship.file,
                             ploidy = ploidy,
                             matings = "all",
                             standardize = TRUE)

  max.parent = ceiling(as.numeric(length(ans1$parents$id)) * selection.intensity)

  ans2 = COMA::ocs(parents = data.frame(ans1$parents, min = 0, max = 1/max.parent),
                   ploidy = ploidy,
                   K = ans1$K,
                   dF = dF[2],
                   dF.adapt = dF.adapt,
                   solver = solver,
                   base = base)

  if (nrow(ans2$oc) == 0) {
    stop("No solution possible.")
  }

  sel1 = ans2$oc$id[order(ans2$oc$value, decreasing = T)]

  if (length(sel1) > max.parent)
    sel1 = sel1[1:max.parent]

  ans1 = read_data_optimized(geno.file = geno.file,
                             kinship.file = kinship.file,
                             ploidy = ploidy,
                             matings = sel1,
                             standardize = TRUE)

  ans3 = COMA::oma(parents = data.frame(id = sel1, min = 0, max = 1),
                   matings = data.frame(ans1$matings, min = 0, max = 1),
                   ploidy = ploidy,
                   K = ans1$K,
                   dF = dF,
                   dF.adapt = dF.adapt,
                   solver = solver,
                   base = base)

  if (nrow(ans3$om) == 0) {
    stop("No solution possible.")
  }

  return(ans3)
}
