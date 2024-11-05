#' PRISM.Pedigree
#'
#' PRISM is the standard database system in the commercial sector (or a variant based off of it). PRISM is capable of exporting Pedigree History from a defined History Start. This function will allow you to create a complete three column pedigree dataframe from the PRISM export to use with COMA.file.2 for generating a pedigree kinship matrix. The PRISM export must contain the following columns; "Pedigree", "PedId", "Female PedId", and "Male PedId".
#' @param data CSV of PRISM exported pedigree history.
#' @return Three column pedigree dataframe properly formated.
#' @export
PRISM.Pedigree = function(data) {
  # Read the data
  data <- fread(data)

  # Clean and prepare data
  data <- data %>%
    filter(!grepl(" X ", Pedigree)) %>%
    select(Pedigree, PedId, Female_PedId = `Female PedId`, Male_PedId = `Male PedId`) %>%
    mutate(Male_PedId = ifelse(Male_PedId == 0, Female_PedId, Male_PedId)) %>%
    distinct()

  # Create lookup for Pedigree to PedId
  lookup <- data %>% select(Pedigree, PedId)

  # Join and reshape
  data_updated <- data %>%
    left_join(lookup, by = c("Female_PedId" = "PedId"), suffix = c("", "_female")) %>%
    left_join(lookup, by = c("Male_PedId" = "PedId"), suffix = c("", "_male")) %>%
    transmute(id = Pedigree, parent1 = Pedigree_female, parent2 = Pedigree_male) %>%
    distinct()

  return(data_updated)
}
