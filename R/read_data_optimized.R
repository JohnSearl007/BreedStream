#' read_data_optimized
#'
#' The read_data() function from COMA, presents limitations for large datasets due to the way that it utilizes RAM. The read_data_optimized() function features the same functionality and output as read_data() but breaks the calculation of merit for the desired matings into chunks and writes the intermediary output to disk for better memory allocation.
#' @param geno.file Genotype information for COMA.
#' @param kinship.file Kinship matrix for COMA.
#' @param ploidy Ploidy of species.
#' @param sex Sex arguement for COMA.
#' @param matings Matings arguement for COMA.
#' @param standardize Should COMA standardize.
#' @param n.core The number of cores for parelle processing.
#' @param partition Partition arguement for COMA.
#' @param chunk_size How many matings pairwise combinations can be handled per core.
#' @param result_file Filename for the output.
#' @return Input data object for COMA.
#' @export
read_data_optimized <- function(geno.file, kinship.file, ploidy, sex=NULL,
                                matings="none", standardize=FALSE,
                                n.core=4, partition=FALSE, chunk_size=100000, result_file="matings_results.csv") {

  # Helper function to calculate MPH (mid-parent heterosis)
  MPH <- function(parents, ploidy, geno) {
    Xi <- geno[,parents[1]]
    Xj <- geno[,parents[2]]
    ploidy/4/(ploidy-1)*((Xi-Xj)^2 + 2/ploidy*Xi*Xj - (Xi+Xj))
  }

  # Function to process each chunk in parallel
  process_chunk_parallel <- function(matings_chunk, ploidy, geno, effects, n.core) {
    mate.list <- split(as.matrix(matings_chunk[, 1:2]), f=1:nrow(matings_chunk))

    cl <- makeCluster(n.core)

    # Export necessary variables to worker nodes
    clusterExport(cl, varlist=c("MPH", "ploidy", "geno", "effects"), envir=environment())

    # Run the parallel calculation
    ans <- parSapply(cl, mate.list, MPH, ploidy=ploidy, geno=geno)

    stopCluster(cl)

    return(ans)
  }

  # Load geno data and kinship matrix using read.csv due to input format constraints
  data <- read.csv(file = geno.file, check.names=FALSE, row.names=1)
  dominance <- (colnames(data)[2] == "dom")
  geno.start <- 2 + as.integer(dominance)

  effects <- as.matrix(data[,1:(geno.start-1), drop=FALSE])
  geno <- as.matrix(data[,geno.start:ncol(data)])
  p <- apply(geno, 1, mean, na.rm=TRUE) / ploidy

  rownames(geno) <- rownames(data)
  id <- colnames(geno)

  Pmat <- kronecker(matrix(p, nrow=1, ncol=nrow(geno)), matrix(1, nrow=ncol(geno), ncol=1))
  coeff <- t(geno) - ploidy * Pmat
  dimnames(coeff) <- list(id, rownames(geno))
  coeff[is.na(coeff)] <- 0

  # Load kinship matrix
  K <- as.matrix(read.csv(kinship.file, row.names=1, check.names=FALSE))
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

  # Initialize result file with an empty data frame
  empty_df <- data.table(female=character(), male=character(), merit=numeric())
  if (dominance) {
    empty_df[, `:=` (MPA=numeric(), MPD=numeric(), MPH=numeric())]
  }

  # Write the empty data frame to initialize the result file
  fwrite(empty_df, file=result_file)

  # Process matings in chunks
  for (i in seq(1, nrow(matings), by=chunk_size)) {
    mating_chunk <- matings[i:min(i + chunk_size - 1, nrow(matings)), ]
    mating_chunk$merit <- (parents$add[match(mating_chunk$female, parents$id)] +
                             parents$add[match(mating_chunk$male, parents$id)])/2

    if (dominance) {
      mating_chunk$MPA <- mating_chunk$merit
      mating_chunk$MPD <- (parents$dom[match(mating_chunk$female, parents$id)] +
                             parents$dom[match(mating_chunk$male, parents$id)])/2

      # Parallel merit calculation
      ans <- process_chunk_parallel(mating_chunk, ploidy, geno, effects, n.core)
      mating_chunk$MPH <- as.numeric(crossprod(ans, effects[,2]))

      if (standardize) {
        mating_chunk$MPA <- mating_chunk$MPA / sd.merit
        mating_chunk$MPD <- mating_chunk$MPD / sd.merit
        mating_chunk$MPH <- mating_chunk$MPH / sd.merit
      }
      mating_chunk$merit <- mating_chunk$MPA + mating_chunk$MPD + mating_chunk$MPH
    } else if (standardize) {
      mating_chunk$merit <- mating_chunk$merit / sd.merit
    }

    # Append chunk results to the file
    fwrite(mating_chunk, file=result_file, append=TRUE)
  }

  # Read and Format Matings
  matings = fread(result_file, select = c("female", "male", "merit"))
  colnames(matings) = c("parent1", "parent2", "merit")

  # Return parents, K matrix, and matings
  return(list(K=K, parents=parents[, setdiff(colnames(parents), c("add", "dom"))], matings=matings))
}
