#' COMA.file.2
#'
#' Creates input file 2 for COMA, containing the pedigree kinship matrix.
#' Handles missing values in `pedigree.data` by replacing `NA`s with `0`s.
#'
#' @param pedigree.data A data frame with columns "id", "parent1", and "parent2" for pedigree information.
#' @param ploidy Integer specifying the ploidy level of the species (default is 2).
#' @return A kinship matrix representing File type 2 for COMA, saved as "COMA_file_2.csv".
#'
#' @importFrom AGHmatrix Amatrix
#' @importFrom utils write.csv
#' @export

COMA.file.2 <- function(pedigree.data, ploidy = 2) {

  # Replace NA values with 0 to handle missing parent information
  pedigree.data[is.na(pedigree.data)] <- 0

  # Generate the kinship matrix using the Amatrix function with specified ploidy
  A <- AGHmatrix::Amatrix(data = pedigree.data, ploidy = ploidy)

  # Write the resulting matrix to a CSV file
  utils::write.csv(A, file = "COMA_file_2.csv", quote = FALSE)

  return(A)  # Return the kinship matrix
}
