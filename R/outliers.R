#' outliers
#'
#' Interactive detection of outliers based on histograms. Writes a new dataframe with erroneous data containing plots set to NA. Multiple traits can be evaluated at once using trait = c().
#' @param data Raw input file to be cleaned.
#' @param trait Character of trait(s) to be cleaned in the input file.
#' @return Dataframe that is properly formated and has had the yield calculated.
#' @export
outliers = function(data, trait) {
  traits = trait
  df = data

  for (trait in trait) {
    hist(as.numeric(unlist(df[, ..trait])), main = paste0("Histogram of: ",trait), xlab = trait)
    exit = "N"

    while (!(exit == "Y")) {
      upper.bound = as.numeric(readline("Enter the upper bound: "))
      lower.bound = as.numeric(readline("Enter the lower bound: "))

      temp = subset(df, unlist(df[, ..trait]) >= lower.bound & unlist(df[, ..trait]) <= upper.bound)
      hist(as.numeric(unlist(temp[, ..trait])), main = paste0("Histogram of: ",trait), xlab = trait)

      exit = as.character(readline("Would you like to exit (Y/N)?: "))
    }

    upper.cutoff = as.numeric(readline("Enter your final upper cutoff: "))
    lower.cutoff = as.numeric(readline("Enter your final lower cutoff: "))

    df = subset(df, unlist(df[, ..trait]) >= lower.cutoff & unlist(df[, ..trait]) <= upper.cutoff)
  }

  df2 = data %>% select(!any_of(traits))
  clean = left_join(df2, df)
  clean[,"Quick Note"] = NA

  return(clean)
}
