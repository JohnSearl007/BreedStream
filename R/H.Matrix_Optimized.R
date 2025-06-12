#' H.Matrix_Optimized
#'
#' Optimizes the blending weight for generating the H Matrix to minimize the AIC, by testing various proportions
#' of blending. This function is designed to help select an optimal blending weight for a mixed model in StageWise
#' based on AIC. Convergence issues are generally related to the blending weight range, not iteration limits.
#'
#' @param pheno.data Data frame containing phenotype data for StageWise.
#' @param geno.data Data frame containing genotype data for StageWise.
#' @param pedigree.data Data frame containing pedigree data for StageWise.
#' @param trait Character vector of trait(s) to be passed to StageWise.
#' @param Fixed Formula or list of fixed effects to be passed to StageWise.
#' @param Random Formula or list of random effects to be passed to StageWise.
#' @param blend.lower Lower limit of the A matrix contribution for blending (default is 1e-5).
#' @param blend.upper Upper limit of the A matrix contribution for blending (default is 0.20).
#' @param ploidy Integer specifying the ploidy level of the species (default is 2).
#' @param max.iter Maximum number of iterations for ASReml (default is 100).
#' @param map Logical; TRUE if genotype data includes map coordinates in the first three columns (default is TRUE).
#' @param dominance Logical; TRUE if dominance effects should be calculated with StageWise (default is FALSE).
#' @param workspace Memory allocation for ASReml (default is "12Gb").
#' @param pworkspace Memory allocation for ASReml's parallel workspace (default is "12Gb").
#' @param min.minor.allele Minimum number of individuals with the minor allele for StageWise (default is 5).
#' @param fix.eff.marker Marker effect argument for StageWise (default is NULL).
#' @param covariates Covariate argument for StageWise (default is NULL).
#' @param mask Logical; if TRUE, StageWise masks some of the data (default is NULL).
#' @param method Method argument for StageWise (default is NULL).
#' @param non.add Specifies handling of non-additive genetic effects (default is "none").
#' @param tol Tolerance level for optimization convergence (default is 1e-4).
#' @return A StageWise prep object containing the optimal H Matrix with minimized AIC and BLUPs prepared.
#'
#' @importFrom asreml asreml.options
#' @importFrom StageWise read_geno Stage1 Stage2 blup_prep
#' @importFrom utils write.table
#' @importFrom stats optimize
#' @export

H.Matrix_Optimized <- function(pheno.data, geno.data, pedigree.data, trait, Fixed, Random,
         blend.lower = 1e-5, blend.upper = 0.20, ploidy = 2, max.iter = 100,
         map = TRUE, dominance = FALSE, workspace = "12Gb", pworkspace = "12Gb",
         min.minor.allele = 5, fix.eff.marker = NULL, covariates = NULL,
         mask = NULL, method = NULL, non.add = "none", tol = 1e-4) {

  # Check for valid blending weight range
  if (blend.lower > blend.upper) {
    stop("Error: blend.lower must be less than or equal to blend.upper")
  }

  # Extract effects and set ASReml options
  effects <- model.terms(Fixed = Fixed, Random = Random)
  asreml::asreml.options(ai.sing = TRUE, workspace = workspace, pworkspace = pworkspace)

  # Stage 1
  model <- suppressMessages(suppressWarnings(
    StageWise::Stage1(filename = pheno.data, traits = trait, effects = effects,
                      solver = "asreml", workspace = c(workspace, pworkspace))
  ))

  # Define function to optimize (only used if optimization is needed)
  optimize_blend <- function(w) {
    tryCatch({
      geno <- suppressMessages(suppressWarnings(
        StageWise::read_geno(geno.data, ploidy = ploidy, map = map,
                             min.minor.allele = min.minor.allele,
                             ped = pedigree.data, w = w, dominance = dominance)
      ))
      ans <- suppressMessages(suppressWarnings(
        StageWise::Stage2(data = model$blues, vcov = model$vcov, geno = geno,
                          non.add = non.add, workspace = c(workspace, pworkspace),
                          max.iter = max.iter, pairwise = TRUE,
                          fix.eff.marker = fix.eff.marker, covariates = covariates)
      ))
      cat(sprintf("Blending weight (w): %.6f, AIC: %.2f\n", w, ans$aic))
      return(ans$aic)
    }, error = function(e) {
      cat(sprintf("Blending weight (w): %.6f, AIC: Inf (convergence failed)\n", w))
      return(Inf)
    })
  }

  # Determine optimal blending weight
  if (blend.lower < blend.upper) {
    # Optimize blending weight when a range is provided
    opt_result <- optimize(optimize_blend, interval = c(blend.lower, blend.upper),
                           maximum = FALSE, tol = tol)
    optimal_weight <- opt_result$minimum
    print(paste0("Optimal blending weight: ", optimal_weight))
  } else {
    # Use the single user-provided weight when blend.lower equals blend.upper
    optimal_weight <- as.numeric(blend.lower)
    print(paste0("Using user-provided blending weight: ", optimal_weight))
  }

  # Save the optimal weight to file
  utils::write.table(optimal_weight, file = "optimal_weight.txt",
                     row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Generate final H matrix with optimal weight
  genoH <- suppressMessages(suppressWarnings(
    StageWise::read_geno(geno.data, ploidy = ploidy, map = map,
                         min.minor.allele = min.minor.allele,
                         ped = pedigree.data, w = optimal_weight, dominance = dominance)
  ))

  # Run Stage 2 with optimal H
  final_ans <- suppressMessages(suppressWarnings(
    StageWise::Stage2(data = model$blues, vcov = model$vcov, geno = genoH,
                      non.add = non.add, workspace = c(workspace, pworkspace),
                      max.iter = max.iter, pairwise = TRUE,
                      fix.eff.marker = fix.eff.marker, covariates = covariates)
  ))

  # Prepare BLUPs
  prep <- suppressMessages(suppressWarnings(
    StageWise::blup_prep(data = model$blues, vcov = model$vcov, geno = genoH,
                         vars = final_ans$vars, mask = mask, method = method)
  ))

  return(prep)
}

