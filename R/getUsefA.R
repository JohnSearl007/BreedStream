#' Predict Usefulness for Crosses with Additive Trait Effects
#'
#' Based on the SimpleMating::getUsefA() function.
#' Predicts the usefulness (mean performance and genetic variance) for a set of crosses
#' based on additive marker effects, accounting for linkage via a genetic map or
#' linkage disequilibrium matrix. Supports Doubled-Haploids (DH) and Recombinant
#' Inbred Lines (RIL) mating systems.
#'
#' @param MatePlan A data frame with at least two columns specifying the parent pairs for crosses (Parent1, Parent2).
#' @param Markers A numeric matrix of marker genotypes for candidate parents, coded as 0, 2, or NA (missing).
#' @param addEff A numeric vector of additive marker effects, with names matching marker names in \code{Markers}.
#' @param K A numeric relationship matrix with row and column names matching genotypes in \code{MatePlan}.
#' @param Map.In A data frame with at least three columns: chromosome identifier, position (e.g., cM), and marker ID.
#' @param linkDes An optional numeric linkage disequilibrium matrix; used if \code{Map.In} is NULL.
#' @param propSel A numeric value between 0 and 1 specifying the proportion of individuals selected (default: 0.05).
#' @param Type A character string specifying the mating system: "DH" (Doubled-Haploids) or "RIL" (Recombinant Inbred Lines).
#' @param Generation A positive integer indicating the generation for DH or RIL populations (default: 1).
#'
#' @return A list of two data frames:
#' \itemize{
#'   \item{[1]}{A data frame with columns: \code{Parent1}, \code{Parent2}, \code{Cross.ID}, \code{Mean} (predicted mean), \code{Variance} (genetic variance), \code{sd} (standard deviation), \code{Usefulness} (mean + selection intensity * sd), and \code{K} (kinship).}
#'   \item{[2]}{A data frame with columns: \code{Parent1}, \code{Parent2}, \code{Y} (usefulness), and \code{K} (kinship), sorted by \code{Y} in descending order, excluding NA values.}
#' }
#'
#' @details
#' The function computes the mean performance as the average breeding value of the parents,
#' imputes missing marker data with column means, and estimates genetic variance using
#' recombination matrices (via \code{Map.In}) or linkage disequilibrium (\code{linkDes}).
#' Variance calculations incorporate segregating markers and account for mating type and generation.
#' Usefulness is calculated as the mean plus the selection intensity (based on \code{propSel})
#' times the standard deviation.
#' @export

getUsefA <- function(MatePlan, Markers, addEff, K, Map.In, linkDes = NULL,
                     propSel = 0.05, Type = 'DH', Generation = 1) {

  # Input validation
  if (!inherits(MatePlan, "data.frame")) {
    stop("Argument 'MatePlan' must be a data.frame.")
  }
  gnames <- unique(c(MatePlan[, 1], MatePlan[, 2]))
  if (!all(gnames %in% rownames(Markers))) {
    stop("Some individuals in 'MatePlan' are missing in 'Markers'.")
  }
  if (!is.numeric(propSel) || propSel <= 0 || propSel >= 1) {
    stop("Argument 'propSel' must be a numeric value between 0 and 1.")
  }
  if (!is.matrix(Markers)) {
    stop("Argument 'Markers' must be a matrix.")
  }
  if (!all(Markers %in% c(0, 2, NA), na.rm = TRUE)) {
    stop("Markers matrix must contain only 0, 2, or NA.")
  }
  if (is.null(linkDes) && is.null(Map.In)) {
    stop("Provide either 'linkDes' or 'Map.In'.")
  }
  if (!is.null(Map.In) && ncol(Map.In) < 3) {
    stop("'Map.In' must have at least 3 columns: chromosome, position, marker ID.")
  }
  if (length(addEff) != ncol(Markers)) {
    stop("'addEff' length must match number of columns in 'Markers'.")
  }

  # Handle addEff names
  if (is.null(names(addEff))) {
    warning("'addEff' is unnamed; assuming positional alignment with 'Markers' columns.")
    names(addEff) <- colnames(Markers)
  }

  # Prepare MatePlan
  MatePlan$idcross <- paste0(MatePlan[, 1], "_", MatePlan[, 2])
  colnames(MatePlan) <- c("Parent1", "Parent2", "Cross.ID")

  # Align markers with Map.In
  if (!is.null(Map.In)) {
    map_markers <- Map.In[, 3]
    marker_names <- colnames(Markers)
    common_markers <- intersect(map_markers, marker_names)
    if (length(common_markers) == 0) {
      if (nrow(Map.In) == ncol(Markers)) {
        warning("No common marker IDs between 'Map.In' and 'Markers'; assuming positional alignment.")
        Map.In[, 3] <- marker_names
        common_markers <- marker_names
      } else {
        stop("No common markers and lengths differ between 'Map.In' and 'Markers'.")
      }
    } else {
      # Subset to common markers
      Markers <- Markers[, common_markers, drop = FALSE]
      addEff <- addEff[common_markers]
      Map.In <- Map.In[Map.In[, 3] %in% common_markers, ]
      Map.In <- Map.In[match(common_markers, Map.In[, 3]), ]
    }
  }

  # Impute missing marker values and compute means
  Markers <- apply(Markers, 2, function(wna) {
    sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm = TRUE), ina))
  })
  est.bredv <- Markers %*% addEff
  MatePlan$Mean <- apply(MatePlan, 1, function(tmp) {
    mean_val <- (est.bredv[tmp[1], 1] + est.bredv[tmp[2], 1]) / 2
    if (!is.finite(mean_val)) mean_val <- 0
    round(mean_val, 5)
  })

  # Process based on linkage information
  if (is.null(linkDes)) {
    # Using genetic map
    Markers_names <- colnames(Markers)

    # Split by chromosome, ensuring consistent ordering
    chroms <- unique(Map.In[, 1])
    Map.Chr <- lapply(chroms, function(ch) Map.In[Map.In[, 1] == ch, ])
    names(Map.Chr) <- as.character(chroms)
    Map.Pos <- lapply(chroms, function(ch) Map.In[Map.In[, 1] == ch, 3])
    names(Map.Pos) <- as.character(chroms)
    Map.Eff <- lapply(chroms, function(ch) addEff[Map.In[Map.In[, 1] == ch, 3]])
    names(Map.Eff) <- as.character(chroms)

    # Compute recombination matrices
    rMat <- lapply(Map.Chr, theta)
    names(rMat) <- as.character(chroms)

    # Compute MCov based on Type and Generation
    if (Type == "DH") {
      if (Generation == 1) {
        MCov <- lapply(rMat, function(ctheta) 1 - (2 * ctheta))
      } else if (Generation >= 10) {
        MCov <- lapply(rMat, function(ctheta) (1 - (2 * ctheta)) / (1 + (2 * ctheta)))
      } else {
        RHS_MCov <- lapply(rMat, function(ctheta) {
          popInfo <- 0.5 * (1 - (2 * ctheta))
          Reduce(`+`, lapply(seq(Generation), function(k) popInfo^k))
        })
        LHS_MCov <- lapply(rMat, function(ctheta) (0.5 * (1 - (2 * ctheta)))^Generation)
        MCov <- Map("+", RHS_MCov, LHS_MCov)
      }
    } else if (Type == "RIL") {
      if (Generation >= 10) {
        MCov <- lapply(rMat, function(ctheta) (1 - (2 * ctheta)) / (1 + (2 * ctheta)))
      } else {
        MCov <- lapply(rMat, function(ctheta) {
          popInfo <- 0.5 * (1 - (2 * ctheta))
          Reduce(`+`, lapply(seq(Generation), function(k) popInfo^k))
        })
      }
    }

    # Ensure MCov is aligned
    MCov <- MCov[as.character(chroms)]

    # Validate MCov dimensions against Map.Pos
    for (chr in names(MCov)) {
      if (is.null(MCov[[chr]]) || is.null(Map.Pos[[chr]]) ||
          any(dim(MCov[[chr]]) != length(Map.Pos[[chr]]))) {
        stop(sprintf("MCov dimensions for chromosome %s (%dx%d) do not match marker count in Map.Pos (%d)",
                     chr, nrow(MCov[[chr]]), ncol(MCov[[chr]]), length(Map.Pos[[chr]])))
      }
    }

    Markers <- Markers - 1

    # Function to calculate information content
    calc.info <- function(Markers) {
      diff <- Markers[1, , drop = FALSE] - Markers[2, , drop = FALSE]
      crossprod(diff) / 4
    }

    # Variance prediction per cross
    crospredPar <- function(Ncross) {
      cross_variance <- vector("list", nrow(Ncross))
      for (i in seq_along(cross_variance)) {
        Matepair <- as.character(Ncross[i, ])
        Total_SNP <- Markers[Matepair, , drop = FALSE]
        SNPseg <- which(!colMeans(Total_SNP) %in% c(1, -1))
        SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_names[SNPseg])

        Pair.Var <- 0
        for (chr in names(Map.Pos)) {
          seg_markers <- SNPseg.Chr[[chr]]
          if (length(seg_markers) > 0) {
            pos <- which(Map.Pos[[chr]] %in% seg_markers)
            if (any(pos > length(Map.Pos[[chr]]))) {
              stop(sprintf("Invalid pos indices for chromosome %s in cross %d: %s",
                           chr, i, paste(pos, collapse = ",")))
            }
            parGen_chr <- Markers[Matepair, seg_markers, drop = FALSE]
            D_chr <- calc.info(parGen_chr)
            MCov_chr <- MCov[[chr]][pos, pos, drop = FALSE]
            VarCov_chr <- D_chr * MCov_chr
            EffA_chr <- Map.Eff[[chr]][pos]
            var_contrib <- crossprod(EffA_chr, VarCov_chr %*% EffA_chr)
            if (!is.finite(var_contrib)) {
              warning(sprintf("Non-finite variance contribution for cross %d, chromosome %s", i, chr))
              var_contrib <- 0
            }
            Pair.Var <- Pair.Var + var_contrib
          }
        }
        cross_variance[[i]] <- data.frame(Parent1 = Matepair[1], Parent2 = Matepair[2],
                                          Variance = abs(Pair.Var),
                                          stringsAsFactors = FALSE, row.names = NULL)
      }
      do.call("rbind", cross_variance)
    }

    # Compute variances
    cros2cores <- list(`1` = MatePlan[, c(1:2)])
    tmp_var <- lapply(cros2cores, crospredPar)
    MateVar <- do.call("rbind", tmp_var)
    MateVar$Cross.ID <- paste0(MateVar$Parent1, "_", MateVar$Parent2)
    MatePlan <- merge(MatePlan, MateVar[, c("Cross.ID", "Variance")], by = "Cross.ID", all.x = TRUE)

    # Fill NA variances with 0
    MatePlan$Variance[is.na(MatePlan$Variance)] <- 0
  } else {
    # Using linkage disequilibrium matrix
    Markers_names <- colnames(Markers)
    map_markers <- Map.In[, 2]
    common_markers <- intersect(map_markers, Markers_names)
    if (length(common_markers) == 0) {
      if (nrow(Map.In) == ncol(Markers)) {
        warning("No common marker IDs; assuming positional alignment for 'linkDes'.")
        Map.In[, 2] <- Markers_names
        common_markers <- Markers_names
      } else {
        stop("No common markers and lengths differ for 'linkDes'.")
      }
    }
    Markers <- Markers[, common_markers, drop = FALSE]
    addEff <- addEff[common_markers]
    Map.In <- Map.In[Map.In[, 2] %in% common_markers, ]
    Map.In <- Map.In[match(common_markers, Map.In[, 2]), ]
    rownames(linkDes) <- colnames(linkDes) <- common_markers

    Map.Pos <- lapply(unique(Map.In[, 1]), function(ch) Map.In[Map.In[, 1] == ch, 2])
    names(Map.Pos) <- as.character(unique(Map.In[, 1]))
    Map.Eff <- lapply(unique(Map.In[, 1]), function(ch) addEff[Map.In[Map.In[, 1] == ch, 2]])
    names(Map.Eff) <- as.character(unique(Map.In[, 1]))
    block_sizes <- table(Map.In[, 1])
    rMat <- list()
    for (i in seq_along(block_sizes)) {
      iMat <- sum(block_sizes[1:(i-1)]) + 1
      jMat <- sum(block_sizes[1:i])
      rMat[[i]] <- linkDes[iMat:jMat, iMat:jMat]
    }
    names(rMat) <- names(block_sizes)

    # Compute MCov
    if (Type == "DH") {
      if (Generation == 1) {
        MCov <- lapply(rMat, function(ctheta) 1 - (2 * ctheta))
      } else if (Generation >= 10) {
        MCov <- lapply(rMat, function(ctheta) (1 - (2 * ctheta)) / (1 + (2 * ctheta)))
      } else {
        RHS_MCov <- lapply(rMat, function(ctheta) {
          popInfo <- 0.5 * (1 - (2 * ctheta))
          Reduce(`+`, lapply(seq(Generation), function(k) popInfo^k))
        })
        LHS_MCov <- lapply(rMat, function(ctheta) (0.5 * (1 - (2 * ctheta)))^Generation)
        MCov <- Map("+", RHS_MCov, LHS_MCov)
      }
    } else if (Type == "RIL") {
      if (Generation >= 10) {
        MCov <- lapply(rMat, function(ctheta) (1 - (2 * ctheta)) / (1 + (2 * ctheta)))
      } else {
        MCov <- lapply(rMat, function(ctheta) {
          popInfo <- 0.5 * (1 - (2 * ctheta))
          Reduce(`+`, lapply(seq(Generation), function(k) popInfo^k))
        })
      }
    }

    MCov <- MCov[as.character(names(block_sizes))]

    # Validate MCov dimensions
    for (chr in names(MCov)) {
      if (is.null(MCov[[chr]]) || is.null(Map.Pos[[chr]]) ||
          any(dim(MCov[[chr]]) != length(Map.Pos[[chr]]))) {
        stop(sprintf("MCov dimensions for chromosome %s (%dx%d) do not match marker count in Map.Pos (%d)",
                     chr, nrow(MCov[[chr]]), ncol(MCov[[chr]]), length(Map.Pos[[chr]])))
      }
    }

    Markers <- Markers - 1

    calc.info <- function(Markers) {
      diff <- Markers[1, , drop = FALSE] - Markers[2, , drop = FALSE]
      crossprod(diff) / 4
    }

    crospredPar <- function(Ncross) {
      cross_variance <- vector("list", nrow(Ncross))
      for (i in seq_along(cross_variance)) {
        Matepair <- as.character(Ncross[i, ])
        Total_SNP <- Markers[Matepair, , drop = FALSE]
        SNPseg <- which(!colMeans(Total_SNP) %in% c(1, -1))
        SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_names[SNPseg])

        Pair.Var <- 0
        for (chr in names(Map.Pos)) {
          seg_markers <- SNPseg.Chr[[chr]]
          if (length(seg_markers) > 0) {
            pos <- which(Map.Pos[[chr]] %in% seg_markers)
            if (any(pos > length(Map.Pos[[chr]]))) {
              stop(sprintf("Invalid pos indices for chromosome %s in cross %d: %s",
                           chr, i, paste(pos, collapse = ",")))
            }
            parGen_chr <- Markers[Matepair, seg_markers, drop = FALSE]
            D_chr <- calc.info(parGen_chr)
            MCov_chr <- MCov[[chr]][pos, pos, drop = FALSE]
            VarCov_chr <- D_chr * MCov_chr
            EffA_chr <- Map.Eff[[chr]][pos]
            var_contrib <- crossprod(EffA_chr, VarCov_chr %*% EffA_chr)
            if (!is.finite(var_contrib)) {
              warning(sprintf("Non-finite variance contribution for cross %d, chromosome %s", i, chr))
              var_contrib <- 0
            }
            Pair.Var <- Pair.Var + var_contrib
          }
        }
        cross_variance[[i]] <- data.frame(Parent1 = Matepair[1], Parent2 = Matepair[2],
                                          Variance = abs(Pair.Var),
                                          stringsAsFactors = FALSE, row.names = NULL)
      }
      do.call("rbind", cross_variance)
    }

    cros2cores <- list(`1` = MatePlan[, c(1:2)])
    tmp_var <- lapply(cros2cores, crospredPar)
    MateVar <- do.call("rbind", tmp_var)
    MateVar$Cross.ID <- paste0(MateVar$Parent1, "_", MateVar$Parent2)
    MatePlan <- merge(MatePlan, MateVar[, c("Cross.ID", "Variance")], by = "Cross.ID", all.x = TRUE)

    # Fill NA variances with 0
    MatePlan$Variance[is.na(MatePlan$Variance)] <- 0
  }

  # Compute usefulness
  MatePlan$sd <- sqrt(MatePlan$Variance)
  selin <- dnorm(qnorm(1 - propSel)) / propSel
  calcuf <- function(x) {
    mean <- as.numeric(x[["Mean"]])
    variance <- as.numeric(x[["Variance"]])
    std <- selin * sqrt(variance)
    if (!is.finite(mean) || !is.finite(std)) {
      return(NA)
    }
    round(mean + std, 5)
  }
  MatePlan$Usefulness <- apply(MatePlan, 1, calcuf)
  MatePlan <- MatePlan[order(MatePlan$Usefulness, decreasing = TRUE), ]
  rownames(MatePlan) <- NULL

  # Prepare output with kinship
  melted_rel <- meltKUsef(K)
  par_K <- data.frame(Cross.ID = paste0(melted_rel$Parent2, "_", melted_rel$Parent1),
                      K = melted_rel$K)
  df <- merge(MatePlan, par_K, by = "Cross.ID", all.x = TRUE)
  crosses2opt <- list(MatePlan)
  crosses2opt[[2]] <- data.frame(Parent1 = df$Parent1, Parent2 = df$Parent2,
                                 Y = df$Usefulness, K = df$K)
  crosses2opt[[2]] <- crosses2opt[[2]][order(crosses2opt[[2]]$Y, decreasing = TRUE), ]
  crosses2opt[[2]] <- crosses2opt[[2]][!is.na(crosses2opt[[2]]$Y), ]
  rownames(crosses2opt[[2]]) <- NULL

  cat("Usefulness predicted for ", nrow(MatePlan), " crosses.\n")

  return(crosses2opt)
}

# Helper functions
theta <- function(map) {
  distMat <- as.matrix(dist(map[, 2], upper = TRUE, diag = TRUE, method = "manhattan"))
  0.5 * (1 - exp(-2 * (distMat / 100)))
}

meltKUsef <- function(X) {
  namesK <- rownames(X)
  X <- cbind(which(!is.na(X), arr.ind = TRUE), na.omit(as.vector(X)))
  X <- as.data.frame(X)
  X[, 1] <- namesK[X[, 1]]
  X[, 2] <- namesK[X[, 2]]
  colnames(X) <- c("Parent1", "Parent2", "K")
  rownames(X) <- NULL
  X
}
