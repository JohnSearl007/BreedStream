#' outliers
#'
#' Detects outliers interactively by displaying histograms of specified traits, allowing the user to set upper and lower bounds iteratively.
#' The function updates the data by removing outliers and sets the values of flagged outlier data to NA.
#' Multiple traits can be evaluated simultaneously by passing them as a character vector in `trait`.
#'
#' @param data Dataframe containing raw input data to be cleaned.
#' @param trait Character vector of trait(s) to be cleaned in the input file.
#' @return Dataframe with outliers set to NA for the specified traits.
#'
#' @importFrom dplyr %>% all_of left_join select
#' @export

outliers <- function(data, trait) {
  df <- data  # Copy of input data to avoid modification of original data

  for (trait_name in trait) {
    # Display histogram for the current trait
    hist(as.numeric(df[[trait_name]]), main = paste("Histogram of:", trait_name), xlab = trait_name)
    exit <- "N"

    # Interactive loop to iteratively adjust bounds
    while (toupper(exit) != "Y") {
      # Obtain bounds from user input
      upper_bound <- as.numeric(readline(prompt = "Enter the upper bound: "))
      lower_bound <- as.numeric(readline(prompt = "Enter the lower bound: "))

      # Subset data to values within the specified bounds
      temp <- subset(df, df[[trait_name]] >= lower_bound & df[[trait_name]] <= upper_bound)

      # Display histogram for the filtered data
      hist(as.numeric(temp[[trait_name]]), main = paste("Filtered Histogram of:", trait_name), xlab = trait_name)

      # Prompt user to exit or adjust bounds
      exit <- readline(prompt = "Would you like to exit (Y/N)?: ")
    }

    # Final bounds for the current trait
    upper_cutoff <- as.numeric(readline(prompt = "Enter your final upper cutoff: "))
    lower_cutoff <- as.numeric(readline(prompt = "Enter your final lower cutoff: "))

    # Set outliers to NA in the original dataframe
    df[[trait_name]][df[[trait_name]] < lower_cutoff | df[[trait_name]] > upper_cutoff] <- NA
  }

  # Merge cleaned traits back with original data
  clean <- data %>%
    dplyr::select(-dplyr::all_of(trait))  # Remove original columns for traits being cleaned
  clean <- dplyr::left_join(clean, df, by = names(clean))  # Join cleaned columns back

  # Add "Quick Note" column to flag cleaned entries if desired
  clean$`Quick Note` <- NA

  return(clean)
}
