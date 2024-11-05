#' COMA.file.2
#'
#' I need a function that will create input file 2 for COMA that has the pedigree kinship matrix. The only modification to Amatrix() is that the NA handling and converting to 0 is handled by this function. The input pedigree.data is made outside the function as naming conventions (text separators) are not consistent, requiring the user to make a dataframe with three columns c("id", "parent1", "parent2").
#' @param pedigree.data Pedigree information for StageWise.
#' @param ploidy Ploidy of the species.
#' @return File type 2 required by COMA.
#' @export
COMA.file.2 = function(pedigree.data, ploidy = 2) {
  pedigree.data[is.na(pedigree.data)] <- 0
  A = Amatrix(data = pedigree.data, ploidy = ploidy)

  write.csv(A, file = "COMA_file_2.csv", quote = F)

  return(A)
}
