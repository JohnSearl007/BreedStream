#' H.Matrix_Optimized
#'
#' This function tests the proportion of blending to perform in generating the H Matrix so as to minimize the AIC. (If the model does not converge it is more likely an issue of the range of weight values being tested for blending than that of not having the max.iter high enough. If the model does not converge in less than 50 iterations there is generally a bigger issue at hand.)
#' @param pheno.data Phenotype data for StageWise.
#' @param geno.data Genotype data for StageWise.
#' @param pedigree.data Pedigree data for StageWise.
#' @param trait Character of trait(s) to pass StageWise.
#' @param Fixed Fixed effects to pass StageWise.
#' @param Random Random effects to pass StageWise.
#' @param blend.lower The lower limit of the A matrix contribution for blending.
#' @param blend.upper The upper limit of the A matrix contribution for blending.
#' @param ploidy Ploidy of the species.
#' @param max.iter The maxium number of iterations for ASReml.
#' @param map Whether the genotype data includes map coordinate infromation in the first three columns.
#' @param dominance Should dominance effects be calculated with StageWise.
#' @param workspace Memory allocation arguement for ASReml.
#' @param pworkspace Memory allocation arguement for ASReml.
#' @param min.minor.allele How many individuals must have the minor allele for StageWise.
#' @param fix.eff.marker Marker effect arguement for StageWise.
#' @param covariates Covariate arguement for StageWise.
#' @param mask Should StageWise mask some of the data.
#' @param method Arguement for StageWise.
#' @param numCores How many cores should be used for parrele processing.
#' @param weight.vector.length How many different blending weights should be tested.
#' @param non.add How should non additive genetic effects be handled.
#' @param max_global_size Memory allocation arguement.
#' @return StageWise prep object that has the optimal blending of the H Matrix.
#' @export
H.Matrix_Optimized = function(pheno.data, geno.data, pedigree.data, trait, Fixed, Random, blend.lower = 1e-5,
                              blend.upper = 0.05, ploidy = 2, max.iter = 100, map = TRUE, dominance = FALSE,
                              workspace = "12Gb", pworkspace = "12Gb", min.minor.allele = 5, fix.eff.marker = NULL,
                              covariates = NULL, mask = NULL, method = NULL, numCores = 2, weight.vector.length = 2,
                              non.add = "none", max_global_size = 10e+09) {

  # Set the maximum allowed size for globals
  options(future.globals.maxSize = max_global_size)

  # Set up parallel backend
  plan(multisession, workers = numCores)

  # Extract effects from the model terms
  effects = model.terms(Fixed = Fixed, Random = Random)

  # Set ASReml R options
  asreml.options(ai.sing = TRUE, workspace = workspace, pworkspace = pworkspace)

  # Stage 1: Build the initial model
  model <- Stage1(filename = pheno.data, traits = trait, effects = effects, solver = "asreml",
                  workspace = c(workspace, pworkspace))

  # Create a sequence of blending weights to evaluate
  w.vec <- seq(blend.lower, blend.upper, length.out = weight.vector.length)

  # Function to calculate AIC based on blending weight
  optimize_blend <- function(w) {
    # Generate the blended genomic relationship matrix (H matrix)
    geno <- read_geno(geno.data, ploidy = ploidy, map = map, min.minor.allele = min.minor.allele,
                      ped = pedigree.data, w = w, dominance = dominance)

    # Run Stage 2 for the given H matrix
    ans <- Stage2(data = model$blues, vcov = model$vcov, geno = geno, non.add = non.add,
                  workspace = c(workspace, pworkspace), max.iter = max.iter, pairwise = TRUE,
                  fix.eff.marker = fix.eff.marker, covariates = covariates)

    # Return the AIC value for optimization
    return(ans$aic)
  }

  # Use future_lapply to calculate AICs in parallel
  aic_values <- future_lapply(w.vec, optimize_blend)

  # Find the optimal blending weight with minimum AIC
  optimal_index <- which.min(aic_values)
  optimal_weight <- w.vec[optimal_index]
  print(paste0("Optimal blending weight: ", optimal_weight))
  write.table(optimal_weight, file = "optimal_weight.txt", row.names = F, col.names = F, quote = F)

  # Now generate the H matrix using the optimal blending weight
  genoH <- read_geno(geno.data, ploidy = ploidy, map = map, min.minor.allele = min.minor.allele,
                     ped = pedigree.data, w = optimal_weight, dominance = dominance)

  # Run Stage 2 again with the optimal H matrix
  final_ans <- Stage2(data = model$blues, vcov = model$vcov, geno = genoH, non.add = non.add,
                      workspace = c(workspace, pworkspace), max.iter = max.iter, pairwise = TRUE,
                      fix.eff.marker = fix.eff.marker, covariates = covariates)

  # Prepare for BLUP estimation
  prep <- blup_prep(data = model$blues, vcov = model$vcov, geno = genoH, vars = final_ans$vars,
                    mask = mask, method = method)

  # Return the preparation result
  return(prep)
}
