#' COMA.file.1
#'
#' I need a function that will create input file 1 for COMA that has the marker effects and allele dosage.
#' @param prep StageWise prep object.
#' @param gain.data StageWise gain object.
#' @param geno.data Genotype information for StageWise.
#' @param pedigree.data Pedigree information for StageWise.
#' @param optimal.weight The optimal H Matrix blending weight.
#' @param ploidy Ploidy of the species.
#' @param map Does the Genotype data have map information for StageWise.
#' @param min.minor.allele How many individuals for a minor allele in StageWise.
#' @param dominance Should dominance effects be calculated by StageWise.
#' @return File type 1 required by COMA.
#' @export
COMA.file.1 = function(prep, gain.data, geno.data, pedigree.data, optimal.weight, ploidy = 2, map = TRUE, min.minor.allele = 5, dominance = FALSE) {

  df = gain.data
  index.coeff <- df$table$index
  names(index.coeff) <- df$table$trait

  geno = read_geno(geno.data, ploidy = ploidy, map = map, min.minor.allele = min.minor.allele, ped = pedigree.data,
                   w = optimal.weight, dominance = dominance)

  add.effects = blup(data = prep, geno, what = "AM", index.coeff = index.coeff) %>%
    rename(add = effect)

  marker.effects = if (dominance) {
    dom.effects = blup(data = prep, geno, what = "DM", index.coeff = index.coeff) %>%
      rename(dom = effect)

    left_join(add.effects, dom.effects)
  } else {
    add.effects
  }

  marker.effects = marker.effects %>% select(!any_of(c("chrom", "position")))

  # Columns of Interest & Add Allele Dosage Information
  dosage = fread(geno.data) %>%
    select(!any_of(c("chrom", "position"))) %>%
    { .[, 2:ncol(.)] * (ploidy / max(.[, 2:ncol(.)])) } # Correcting if not coded correctly

  dosage = cbind(marker = add.effects$marker, dosage)

  output = left_join(marker.effects, dosage)

  write.csv(output, file = "COMA_file_1.csv", quote = F, row.names = F)

  return(output)
}
