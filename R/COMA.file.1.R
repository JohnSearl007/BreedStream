#' COMA.file.1
#'
#' Creates input file 1 for COMA, containing marker effects and allele dosage information.
#'
#' @param prep A data frame or list containing preprocessed data for StageWise analysis.
#' @param gain.data A data frame with gain data including index coefficients and trait information.
#' @param geno.data Genotype information used by StageWise.
#' @param pedigree.data Pedigree information for individuals in StageWise.
#' @param optimal.weight A numeric value for the optimal H Matrix blending weight.
#' @param ploidy Integer specifying the ploidy level of the species (default is 2).
#' @param map Logical; if TRUE, genotype data includes map information for StageWise.
#' @param min.minor.allele Minimum count of individuals required for a minor allele.
#' @param dominance Logical; if TRUE, calculates dominance effects in addition to additive effects.
#' @return A data frame representing File type 1 required by COMA, containing marker effects.
#'
#' @importFrom StageWise read_geno blup
#' @importFrom dplyr rename left_join select any_of %>%
#' @importFrom data.table fread
#' @importFrom utils write.csv
#' @export

COMA.file.1 <- function(prep, gain.data, geno.data, pedigree.data, optimal.weight, ploidy = 2, map = TRUE, min.minor.allele = 5, dominance = FALSE) {

  # Extract the index coefficients from the gain.data table
  df <- gain.data
  index.coeff <- df$table$index
  names(index.coeff) <- df$table$trait

  # Read genotype data and adjust based on various parameters
  geno <- StageWise::read_geno(geno.data, ploidy = ploidy, map = map, min.minor.allele = min.minor.allele,
                               ped = pedigree.data, w = optimal.weight, dominance = dominance)

  # Calculate additive effects using the blup function
  add.effects <- StageWise::blup(data = prep, geno, what = "AM", index.coeff = index.coeff) %>%
    dplyr::rename(add = effect)  # Renaming effect column to 'add' for additive effects

  # Optionally calculate dominance effects, if specified
  marker.effects <- if (dominance) {
    dom.effects <- StageWise::blup(data = prep, geno, what = "DM", index.coeff = index.coeff) %>%
      dplyr::rename(dom = effect)  # Renaming effect column to 'dom' for dominance effects

    # Merge additive and dominance effects
    dplyr::left_join(add.effects, dom.effects, by = "marker")
  } else {
    add.effects
  }

  # Remove chromosome and position information if present
  marker.effects <- marker.effects %>% dplyr::select(-dplyr::any_of(c("chrom", "position")))

  # Read in dosage information and adjust for ploidy level
  dosage <- data.table::fread(geno.data) %>%
    dplyr::select(-dplyr::any_of(c("chrom", "position"))) %>%
    { .[, 2:ncol(.)] * (ploidy / max(.[, 2:ncol(.)])) }  # Scaling dosage based on ploidy

  # Combine marker column from additive effects with dosage data
  dosage <- cbind(marker = add.effects$marker, dosage)

  # Merge marker effects and dosage data
  output <- dplyr::left_join(marker.effects, dosage, by = "marker")

  # Write output to a CSV file
  utils::write.csv(output, file = "COMA_file_1.csv", quote = FALSE, row.names = FALSE)

  return(output)  # Return the final combined data frame
}
