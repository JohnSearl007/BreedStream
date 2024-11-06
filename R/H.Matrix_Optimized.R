#' H.Matrix_Optimized
#'
#' Optimizes the blending weight for generating the H Matrix to minimize the AIC, by testing various proportions
#' of blending. This function is designed to help select an optimal blending weight for a mixed model in ASReml-R
#' based on AIC. Convergence issues are generally related to the blending weight range, not iteration limits.
#'
#' @param pheno.data Data frame containing phenotype data for StageWise.
#' @param geno.data Data frame containing genotype data for StageWise.
#' @param pedigree.data Data frame containing pedigree data for StageWise.
#' @param trait Character vector of trait(s) to be passed to StageWise.
#' @param Fixed Formula or list of fixed effects to be passed to StageWise.
#' @param Random Formula or list of random effects to be passed to StageWise.
#' @param blend.lower Lower limit of the A matrix contribution for blending (default is 1e-5).
#' @param blend.upper Upper limit of the A matrix contribution for blending (default is 0.05).
#' @param ploidy Integer specifying the ploidy level of the species (default is 2).
#' @param max.iter Maximum number of iterations for ASReml (default is 100).
#' @param map Logical; TRUE if genotype data includes map coordinates in the first three columns (default is TRUE).
#' @param dominance Logical; TRUE if dominance effects should be calculated with StageWise (default is FALSE).
#' @param workspace Memory allocation for ASReml (default is "12Gb").
#' @param pworkspace Memory allocation for ASReml's parallel workspace (default is "12Gb").
#' @param min.minor.allele Minimum number of individuals with the minor allele for StageWise (default is 5).
#' @param fix.eff.marker Marker effect argument for StageWise.
#' @param covariates Covariate argument for StageWise.
#' @param mask Logical; if TRUE, StageWise masks some of the data.
#' @param method Method argument for StageWise.
#' @param numCores Number of cores for parallel processing (default is 2).
#' @param weight.vector.length Number of blending weights to test (default is 2).
#' @param non.add Specifies handling of non-additive genetic effects (default is "none").
#' @param max_global_size Memory allocation limit for global options (default is 10e+09).
#' @return A StageWise prep object containing the optimal H Matrix with minimized AIC.
#'
#' @importFrom asreml asreml.options
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom StageWise read_geno Stage1 Stage2 blup_prep
#' @importFrom utils write.table
#' @export

H.Matrix_Optimized <- function(pheno.data, geno.data, pedigree.data, trait, Fixed, Random, blend.lower = 1e-5,
                               blend.upper = 0.05, ploidy = 2, max.iter = 100, map = TRUE, dominance = FALSE,
                               workspace = "12Gb", pworkspace = "12Gb", min.minor.allele = 5, fix.eff.marker = NULL,
                               covariates = NULL, mask = NULL, method = NULL, numCores = 2, weight.vector.length = 2,
                               non.add = "none", max_global_size = 10e+09) {

  # Set the maximum allowed size for global variables
  options(future.globals.maxSize = max_global_size)

  # Initialize parallel processing with the specified number of cores
  future::plan(future::multisession, workers = numCores)

  # Extract fixed and random effects from the model terms
  effects <- model.terms(Fixed = Fixed, Random = Random)

  # Set ASReml options for memory allocation
  asreml::asreml.options(ai.sing = TRUE, workspace = workspace, pworkspace = pworkspace)

  # Stage 1: Build the initial model with specified phenotype data and effects
  model <- StageWise::Stage1(filename = pheno.data, traits = trait, effects = effects, solver = "asreml",
                             workspace = c(workspace, pworkspace))

  # Generate a sequence of blending weights within the specified range
  w.vec <- seq(blend.lower, blend.upper, length.out = weight.vector.length)

  # Function to calculate AIC for each blending weight
  optimize_blend <- function(w) {
    # Generate the H matrix using a blended genomic relationship matrix with weight w
    geno <- StageWise::read_geno(geno.data, ploidy = ploidy, map = map, min.minor.allele = min.minor.allele,
                                 ped = pedigree.data, w = w, dominance = dominance)

    # Run Stage 2 to obtain AIC for the given H matrix
    ans <- StageWise::Stage2(data = model$blues, vcov = model$vcov, geno = geno, non.add = non.add,
                             workspace = c(workspace, pworkspace), max.iter = max.iter, pairwise = TRUE,
                             fix.eff.marker = fix.eff.marker, covariates = covariates)

    return(ans$aic)  # Return the AIC value
  }

  # Use parallel processing to evaluate AIC values for each blending weight
  aic_values <- future.apply::future_lapply(w.vec, optimize_blend)

  # Identify the optimal blending weight with the lowest AIC
  optimal_index <- which.min(aic_values)
  optimal_weight <- w.vec[optimal_index]
  print(paste0("Optimal blending weight: ", optimal_weight))
  utils::write.table(optimal_weight, file = "optimal_weight.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Generate the H matrix using the optimal blending weight
  genoH <- StageWise::read_geno(geno.data, ploidy = ploidy, map = map, min.minor.allele = min.minor.allele,
                                ped = pedigree.data, w = optimal_weight, dominance = dominance)

  # Run Stage 2 with the optimal H matrix
  final_ans <- StageWise::Stage2(data = model$blues, vcov = model$vcov, geno = genoH, non.add = non.add,
                                 workspace = c(workspace, pworkspace), max.iter = max.iter, pairwise = TRUE,
                                 fix.eff.marker = fix.eff.marker, covariates = covariates)

  # Prepare for BLUP estimation based on the final model
  prep <- StageWise::blup_prep(data = model$blues, vcov = model$vcov, geno = genoH, vars = final_ans$vars,
                               mask = mask, method = method)

  return(prep)  # Return the final prep object
}
