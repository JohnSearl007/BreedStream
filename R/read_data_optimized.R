#' read_data_optimized
#'
#' This function optimizes the `read_data()` functionality from COMA for large datasets by processing mating merit calculations in chunks,
#' reducing memory usage by saving intermediary results to disk. The output is compatible with COMA's kinship matrix and genetic merit calculations.
#'
#' @param geno.file Character. Path to the genotype file for COMA.
#' @param kinship.file Character. Path to the kinship matrix for COMA.
#' @param ploidy Numeric. Ploidy level of the species.
#' @param sex Dataframe or NULL. Contains 'id' and 'female' columns for COMA (optional).
#' @param matings Character or Dataframe. Specifies mating pairs; can be "none" or "all" for default options.
#' @param standardize Logical. If TRUE, standardizes the genetic merit.
#' @param n.core Integer. Number of cores for parallel processing.
#' @param partition Logical. Argument passed to COMA.
#' @param chunk_size Integer. Number of mating pairs to process per core in each chunk.
#' @param result_file Character. Filename for output.
#' @return List containing the kinship matrix `K`, parental data `parents`, and mating data `matings`.
#'
#' @importFrom data.table data.table fread fwrite
#' @importFrom parallel clusterExport makeCluster parSapply stopCluster
#' @importFrom utils read.csv
#' @export

read_data_optimized <- function(geno.file, kinship.file, ploidy, sex=NULL,
                                matings="none", standardize=FALSE,
                                n.core=4, partition=FALSE, chunk_size=100000, result_file="matings_results.csv") {

  # Helper function to calculate mid-parent heterosis (MPH)
  MPH <- function(parents, ploidy, geno) {
    Xi <- geno[,parents[1]]
    Xj <- geno[,parents[2]]
    ploidy/4/(ploidy-1)*((Xi-Xj)^2 + 2/ploidy*Xi*Xj - (Xi+Xj))
  }

  # Function to process mating chunk in parallel
  process_chunk_parallel <- function(matings_chunk, ploidy, geno, effects, n.core) {
    mate.list <- split(as.matrix(matings_chunk[, 1:2]), f=1:nrow(matings_chunk))
    cl <- parallel::makeCluster(n.core)
    parallel::clusterExport(cl, varlist=c("MPH", "ploidy", "geno", "effects"), envir=environment())
    ans <- parallel::parSapply(cl, mate.list, MPH, ploidy=ploidy, geno=geno)
    parallel::stopCluster(cl)
    return(ans)
  }

  # Load genotype data and kinship matrix
  data <- utils::read.csv(file = geno.file, check.names=FALSE, row.names=1)
  dominance <- (colnames(data)[2] == "dom")
  geno.start <- 2 + as.integer(dominance)
  effects <- as.matrix(data[,1:(geno.start-1), drop=FALSE])
  geno <- as.matrix(data[,geno.start:ncol(data)])
  p <- apply(geno, 1, mean, na.rm=TRUE) / ploidy
  rownames(geno) <- rownames(data)
  id <- colnames(geno)

  # Set up coefficients for additive and dominance effects
  Pmat <- kronecker(matrix(p, nrow=1, ncol=nrow(geno)), matrix(1, nrow=ncol(geno), ncol=1))
  coeff <- t(geno) - ploidy * Pmat
  dimnames(coeff) <- list(id, rownames(geno))
  coeff[is.na(coeff)] <- 0

  K <- as.matrix(utils::read.csv(kinship.file, row.names=1, check.names=FALSE))
  K <- K[id, id]

  parents <- data.frame(id=id, add=as.numeric(coeff %*% effects[,1,drop=FALSE]))

  if (dominance) {
    coeff.D <- -2*choose(ploidy,2)*Pmat^2 + 2*(ploidy-1)*Pmat*t(geno) - t(geno)*(t(geno)-1)
    coeff.D[is.na(coeff.D)] <- 0
    gamma <- (ploidy/2 - 1)/(ploidy - 1)
    parents$dom <- as.numeric(coeff.D %*% effects[,2,drop=FALSE])
    parents$merit <- parents$add + gamma*parents$dom
  } else {
    parents$merit <- parents$add
  }

  sd.merit <- sd(parents$add)
  if (standardize)
    parents$merit <- parents$merit/sd.merit

  if (!is.null(sex)) {
    colnames(sex) <- c("id", "female")
    parents <- merge(parents, sex)
  }

  # Handle 'matings' parameter
  if (is.character(matings) && length(matings) == 1 && matings == "none") {
    return(list(K=K, parents=parents[, setdiff(colnames(parents), c("add", "dom"))]))
  } else if (is.character(matings) && length(matings) == 1 && matings == "all") {
    matings <- parents$id
  }

  if ("female" %in% colnames(parents)) {
    females <- intersect(parents$id[parents$female], matings)
    males <- intersect(parents$id[!parents$female], matings)
    matings <- data.frame(expand.grid(female=females, male=males, stringsAsFactors=FALSE))
  } else {
    id2 <- intersect(matings, parents$id)
    matings <- data.frame(expand.grid(female=id2, male=id2, stringsAsFactors=FALSE))
    matings <- matings[matings$female >= matings$male,]
  }

  matings$female <- as.character(matings$female)
  matings$male <- as.character(matings$male)

  # Initialize the results file
  empty_df <- data.table::data.table(female=character(), male=character(), merit=numeric())
  if (dominance) empty_df[, `:=` (MPA=numeric(), MPD=numeric(), MPH=numeric())]
  data.table::fwrite(empty_df, file=result_file)

  # Process mating chunks and save results iteratively
  for (i in seq(1, nrow(matings), by=chunk_size)) {
    mating_chunk <- matings[i:min(i + chunk_size - 1, nrow(matings)), ]
    mating_chunk$merit <- (parents$add[match(mating_chunk$female, parents$id)] +
                             parents$add[match(mating_chunk$male, parents$id)]) / 2
    if (dominance) {
      mating_chunk$MPA <- mating_chunk$merit
      mating_chunk$MPD <- (parents$dom[match(mating_chunk$female, parents$id)] +
                             parents$dom[match(mating_chunk$male, parents$id)]) / 2

      ans <- process_chunk_parallel(mating_chunk, ploidy, geno, effects, n.core)
      mating_chunk$MPH <- as.numeric(crossprod(ans, effects[,2]))

      if (standardize) {
        mating_chunk$MPA <- mating_chunk$MPA / sd.merit
        mating_chunk$MPD <- mating_chunk$MPD / sd.merit
        mating_chunk$MPH <- mating_chunk$MPH / sd.merit
      }
      mating_chunk$merit <- mating_chunk$MPA + mating_chunk$MPD + mating_chunk$MPH
    } else if (standardize) {
      mating_chunk$merit <- mating_chunk$merit / sd(parents$add)
    }

    data.table::fwrite(mating_chunk, file=result_file, append=TRUE)
  }

  matings <- data.table::fread(result_file, select = c("female", "male", "merit"))
  colnames(matings) <- c("parent1", "parent2", "merit")

  return(list(K=K, parents=parents[, setdiff(colnames(parents), c("add", "dom"))], matings=matings))
}
