#' oma_optimized
#'
#' Optimizes COMA::oma for use with large datasets. The key difference is the use of sparse matrices as opposed to dense matrices.
#' This allows for more efficient RAM allocation enabling the use of larger datasets.
#' Additionally, R has a 32-bit integer limit of 2^31-1 so if ((# parents) * (# matings)) >= 2^31-1,
#' then matings is filtered to maximize the number of matings evaluated while remaining below the 32-bit limit (issue resides in CVXR solver).
#'
#' @param dF Targeted inbreeding rate.
#' @param parents Data frame containing the parents.
#' @param matings Data frame containing specified matings.
#' @param ploidy Numeric value for ploidy.
#' @param K Kinship matrix.
#' @param tol Tolerance threshold.
#' @param dF.adapt Allows for inbreeding rate adjustment.
#' @param solver Solver method for CVXR.
#' @param base Population type.
#' @param sparse_threshold Number of crosses to begin using sparse matrices, user option allows trade off for speed vs RAM efficiency.
#' @return A list containing response data.frame, optimal contributions data.frame, and optimal allocations data.frame.
#'
#' @importFrom CVXR Variable Maximize sum_entries quad_form Problem solve
#' @importFrom Matrix sparse.model.matrix t diag
#' @export

oma_optimized <- function(dF, parents, matings, ploidy, K, tol=1e-6,
                          dF.adapt=NULL, solver="ECOS", base="current", sparse_threshold=1000) {

  requireNamespace("CVXR", quietly=TRUE)
  requireNamespace("Matrix", quietly=TRUE)

  # Initial input summary
  cat("Initial input sizes:\n")
  cat("  Parents:", nrow(parents), "rows\n")
  cat("  Matings:", nrow(matings), "rows\n")
  cat("  Kinship matrix:", dim(K), "\n")
  cat("  Sparse matrix used:", nrow(matings) > sparse_threshold, "\n")

  # Input validation
  if (!is.numeric(dF) || length(dF) > 2L) stop("dF must be numeric with length 1 or 2")
  if (length(dF) == 2L && dF[1] > dF[2]) stop("dF[1] must be <= dF[2]")
  dF <- rep(dF, length.out=2)
  if (!solver %in% c("SCS", "ECOS")) stop("Solver must be 'SCS' or 'ECOS'")
  if (!is.numeric(sparse_threshold) || sparse_threshold <= 0) stop("sparse_threshold must be positive")

  parents <- parents[parents$max > 0, , drop=FALSE]
  matings <- matings[matings$max > 0, , drop=FALSE]
  if (nrow(parents) == 0 || nrow(matings) == 0) stop("No valid parents or matings after filtering")

  if (!all(c("id", "min", "max") == colnames(parents)[1:3])) stop("parents must have columns 'id', 'min', 'max'")
  parents$id <- as.character(parents$id)
  if (any(parents$max > 1 | parents$min < 0)) stop("parents$max must be <= 1, min >= 0")

  if (!all(c("merit", "min", "max") == colnames(matings)[3:5])) stop("matings must have columns 'merit', 'min', 'max'")
  if (any(matings$max > 1 | matings$min < 0)) stop("matings$max must be <= 1, min >= 0")
  colnames(matings) <- c("female", "male", "merit", "min", "max")

  # Filter matings to ensure n * p * (1 + n/p) < 2^31 - 1 (CVXR constraint matrix aprox. size)
  max_elements <- as.integer(0.97 * 2^31 - 1)  # 32-bit limit (2,147,483,647) minus a "fudge-factor" as I can only approximate
  p <- as.numeric(nrow(matings))
  parent_ids <- unique(c(matings$female, matings$male))
  n <- as.numeric(length(parent_ids))
  was_filtered <- FALSE

  if (n * p * (1 + n/p) >= max_elements) {
    # Solve for p: n * p * (1 + n/p) < max_elements
    # This simplifies to: n * p + n * n < max_elements
    p_new <- floor((max_elements - n * n) / n)
    if (p_new < 1) stop("Number of parents (n) too large to satisfy n * p * (1 + n/p) < 2^31-1 even with one mating")

    # Sort matings by merit and keep top p_new
    matings <- matings[order(matings$merit, decreasing = TRUE)[1:p_new], , drop = FALSE]
    p <- nrow(matings)
    parent_ids <- unique(c(matings$female, matings$male))
    n <- length(parent_ids)
    was_filtered <- TRUE
  }

  # Summarize retained data
  cat("Retained data sizes", ifelse(was_filtered, " (after filtering):", ":"), "\n")
  cat("  Parents:", n, "rows\n")
  cat("  Matings:", p, "rows\n")
  cat("  Kinship matrix:", n, "x", n, "\n")
  cat("  Sparse matrix used:", p > sparse_threshold, "\n")

  ploidy <- as.integer(ploidy)
  if (ploidy < 2) stop("ploidy must be >= 2")
  Fi_full <- matrix((ploidy * Matrix::diag(K) - 1) / (ploidy - 1), nrow=1)
  Ft0 <- mean(as.numeric(Fi_full))
  if (toupper(base) == "RM") {
    Ft0 <- (ploidy/2 * mean(K) + (ploidy/2 - 1) * Ft0) / (ploidy - 1)
  }

  parent.id <- sort(parent_ids)
  ix <- match(parent.id, parents$id)
  if (anyNA(ix)) stop("parents data.frame missing individuals in matings")
  parents <- parents[ix, , drop=FALSE]
  K <- K[parent.id, parent.id]

  use_sparse <- p > sparse_threshold
  K <- if (use_sparse) as(K, "dgCMatrix") else as.matrix(K)

  matings$female <- factor(matings$female, levels=parent.id)
  matings$male <- factor(matings$male, levels=parent.id)
  M <- if (use_sparse) {
    Lf <- Matrix::sparse.model.matrix(~female - 1, matings)
    Lm <- Matrix::sparse.model.matrix(~male - 1, matings)
    Matrix::t(Lf + Lm) / 2
  } else {
    Lf <- model.matrix(~female - 1, matings)
    Lm <- model.matrix(~male - 1, matings)
    t(Lf + Lm) / 2
  }
  rownames(M) <- parent.id

  Kvec <- K[cbind(match(matings$female, parent.id), match(matings$male, parent.id))]
  Kvec <- matrix(Kvec, nrow=1)
  Fi <- matrix((ploidy * Matrix::diag(K) - 1) / (ploidy - 1), nrow=1)

  oc <- data.frame(id=parents$id, value=numeric(n))
  om <- data.frame(matings[, 1:2], value=numeric(p))

  y <- CVXR::Variable(n)
  x <- CVXR::Variable(p)
  objective <- CVXR::Maximize(matrix(matings$merit, nrow=1) %*% x)

  con.list <- list(
    y <= parents$max,
    y >= parents$min,
    y == M %*% x,
    x <= matings$max,
    x >= matings$min,
    CVXR::sum_entries(x) == 1
  )

  constraints <- NULL
  if ("female" %in% colnames(parents)) {
    if (anyDuplicated(parents$id) || !is.logical(parents$female))
      stop("parents$id must be unique and female column must be logical")
    if (sum(parents$female) == 0 || sum(!parents$female) == 0)
      stop("At least one female and one male parent required")
    parents$female <- as.integer(parents$female)
    con.list[[7]] <- CVXR::sum_entries(parents$female * y) == 0.5
    sexed <- TRUE
    if (ncol(parents) > 4) constraints <- parents[, 5:ncol(parents), drop=FALSE]
  } else {
    sexed <- FALSE
    if (ncol(parents) > 3) constraints <- parents[, 4:ncol(parents), drop=FALSE]
  }

  if (!is.null(constraints)) {
    tmp <- colnames(constraints)
    signs <- substr(tmp, 1, 2)
    stopifnot(signs %in% c("lt", "gt", "eq"))
    values <- as.numeric(substr(tmp, 3, nchar(tmp)))
    lt <- which(signs == "lt")
    gt <- which(signs == "gt")
    eq <- which(signs == "eq")
    if (length(gt) > 0) {
      values[gt] <- -values[gt]
      constraints[, gt] <- -constraints[, gt]
      signs[gt] <- "lt"
      lt <- which(signs == "lt")
    }
    if (length(lt) > 0) {
      con.list <- c(con.list,
                    list(t(as.matrix(constraints[, lt, drop=FALSE])) %*% y <= values[lt]))
    }
    if (length(eq) > 0) {
      con.list <- c(con.list,
                    list(t(as.matrix(constraints[, eq, drop=FALSE])) %*% y == values[eq]))
    }
  }

  dF2 <- dF[2]
  done <- FALSE
  while (!done) {
    cat("Attempting solve with dF2 =", dF2, "\n")

    Ft1.max <- dF2 + (1 - dF2) * Ft0
    Ft1.min <- dF[1] + (1 - dF[1]) * Ft0
    Ft2 <- dF2 * (2 - dF2) + (1 - dF2)^2 * Ft0

    con2 <- c(con.list,
              list((ploidy/2) * Kvec %*% x + (ploidy/2 - 1) * Fi %*% y <= (ploidy - 1) * Ft1.max,
                   (ploidy/2) * Kvec %*% x + (ploidy/2 - 1) * Fi %*% y >= (ploidy - 1) * Ft1.min),
              list((ploidy/2) * (ploidy - 1) * CVXR::quad_form(y, K) +
                     (ploidy/2 - 1) * ((ploidy/2 - 1) * Fi %*% y + ploidy/2 * Kvec %*% x) <=
                     (ploidy - 1)^2 * Ft2))

    problem <- CVXR::Problem(objective, con2)

    cat("Starting CVXR solver with n =", n, ", p =", p, "\n")
    result <- CVXR::solve(problem, solver=solver, verbose=TRUE)
    cat("Solver finished with status:", result$status, "\n")

    if (result$status %in% c("optimal", "optimal_inaccurate")) {
      done <- TRUE
    } else if (!is.null(dF.adapt) && dF2 < dF.adapt$max) {
      dF2 <- dF2 + dF.adapt$step
      cat("Adjusting dF2 to:", dF2, "\n")
    } else {
      done <- TRUE
    }
  }

  if (result$status %in% c("optimal", "optimal_inaccurate")) {
    if (result$status == "optimal_inaccurate") warning("Optimal inaccurate solution")

    x.opt <- as.numeric(result$getValue(x))
    x.opt[x.opt < tol] <- 0
    x.opt <- x.opt / sum(x.opt)
    y.opt <- as.numeric(M %*% x.opt)

    Ft1 <- as.numeric((ploidy/2) * Kvec %*% x.opt + (ploidy/2 - 1) * Fi %*% y.opt) / (ploidy - 1)
    dF1 <- (Ft1 - Ft0) / (1 - Ft0)
    Ft2 <- ((ploidy/2) * (ploidy - 1) * crossprod(y.opt, as.numeric(K %*% y.opt)) +
              (ploidy/2 - 1) * ((ploidy/2 - 1) * Fi %*% y.opt + ploidy/2 * Kvec %*% x.opt)) / (ploidy - 1)^2
    dF2 <- 1 - sqrt((1 - Ft2) / (1 - Ft0))

    oc$value <- y.opt
    om$value <- x.opt

    om <- om[om$value >= tol, , drop=FALSE]
    oc <- oc[match(oc$id, c(om$female, om$male), nomatch=0) > 0, , drop=FALSE]

    if (!sexed) colnames(om)[1:2] <- c("parent1", "parent2")

    return(list(response=data.frame(dF1=round(dF1, 4), dF2=round(dF2, 4),
                                    merit=result$value, n.parent=nrow(oc), n.mate=nrow(om)),
                oc=oc, om=om))
  } else {
    return(list(response=data.frame(dF1=NA_real_, dF2=NA_real_, merit=NA_real_,
                                    n.parent=0L, n.mate=0L),
                oc=oc[0, ], om=om[0, ]))
  }
}
