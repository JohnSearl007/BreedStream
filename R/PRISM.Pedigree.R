#' PRISM.Pedigree
#'
#' PRISM is a database system commonly used in the commercial sector for pedigree tracking. This function processes an export from PRISM to create a standardized three-column pedigree dataframe
#' (individual, female parent, male parent) suitable for generating a kinship matrix with COMA.file.2.
#' The PRISM export must contain the columns "Pedigree", "PedId", "Female PedId", and "Male PedId".
#'
#' @param data A CSV file containing the PRISM-exported pedigree history.
#' @return A three-column pedigree dataframe with the columns: `id` (individual), `parent1` (female parent), and `parent2` (male parent).
#'
#' @importFrom data.table fread
#' @importFrom dplyr %>% distinct filter left_join mutate select transmute
#' @export

PRISM.Pedigree <- function(data) {
  # Read in the CSV file
  data <- data.table::fread(data)

  # Initial filtering and renaming of columns for consistency
  data <- data %>%
    dplyr::filter(!grepl(" X ", Pedigree)) %>%  # Remove any hybrids with " X " in the Pedigree column
    dplyr::select(Pedigree, PedId, Female_PedId = `Female PedId`, Male_PedId = `Male PedId`) %>%
    dplyr::mutate(Male_PedId = ifelse(Male_PedId == 0, Female_PedId, Male_PedId)) %>%  # Use Female_PedId if Male_PedId is zero
    dplyr::distinct()  # Remove any duplicate rows

  # Create a lookup table to map PedId values to their corresponding Pedigree names
  lookup <- data %>% dplyr::select(Pedigree, PedId)

  # Join with lookup to get the names of the female and male parents based on their PedId values
  data_updated <- data %>%
    dplyr::left_join(lookup, by = c("Female_PedId" = "PedId"), suffix = c("", "_female")) %>%
    dplyr::left_join(lookup, by = c("Male_PedId" = "PedId"), suffix = c("", "_male")) %>%
    dplyr::transmute(id = Pedigree,                # Individual ID
                     parent1 = Pedigree_female,    # Female parent
                     parent2 = Pedigree_male) %>%  # Male parent
    dplyr::distinct()  # Ensure uniqueness

  return(data_updated)
}
