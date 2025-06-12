#' planCross
#'
#' Optimizes SimpleMating::planCross to handle large crossing designs more efficiently. The primary difference is the vectorization of building crosses and not iteratively building the crossing schedule.
#'
#' @param TargetPop Parents for crossing.
#' @param MateDesign Mating design choice.
#' @param TargetPop2 Parents to cross with TargetPop.
#' @param Indiv2keep Filter for parents to use.
#' @return A data.frame object with both parents for a cross.
#'
#' @export

planCross <- function(TargetPop, MateDesign = "half", TargetPop2 = NULL, Indiv2keep = NULL) {
  if (!MateDesign %in% c("full_p", "full", "half", "half_p", "maxAvoid", "circularPlan")) {
    stop("Please, choose a valid mate design plan: full_p, full, half_p, half, maxAvoid, circularPlan")
  }

  # Early filtering
  if (!is.null(Indiv2keep)) TargetPop <- TargetPop[TargetPop %in% Indiv2keep]

  # Two-population case (North Carolina 2 design)
  if (!is.null(TargetPop2)) {
    parent_comb <- as.matrix(expand.grid(TargetPop, TargetPop2))
  } else {
    group1 <- group2 <- TargetPop
    n <- length(group1)

    if (MateDesign == "half") {
      parent_comb <- t(combn(n, 2))
      parent_comb <- cbind(group1[parent_comb[, 1]], group2[parent_comb[, 2]])
    } else if (MateDesign == "half_p") {
      idx <- expand.grid(1:n, 1:n)
      parent_comb <- idx[idx[, 1] <= idx[, 2], ]
      parent_comb <- cbind(group1[parent_comb[, 1]], group2[parent_comb[, 2]])
    } else if (MateDesign == "full") {
      parent_comb <- t(combn(n, 2))
      parent_comb <- rbind(parent_comb, parent_comb[, c(2, 1)])
      parent_comb <- cbind(group1[parent_comb[, 1]], group2[parent_comb[, 2]])
    } else if (MateDesign == "full_p") {
      parent_comb <- as.matrix(expand.grid(1:n, 1:n))
      parent_comb <- cbind(group1[parent_comb[, 1]], group2[parent_comb[, 2]])
    } else if (MateDesign == "maxAvoid") {
      if (n %% 2 != 0) TargetPop <- TargetPop[-n]
      Plan <- matrix(1:length(TargetPop), ncol = 2, byrow = TRUE)
      orderTmp <- TargetPop[c(seq(1, length(TargetPop), by = 2), seq(2, length(TargetPop), by = 2))]
      parent_comb <- cbind(orderTmp[Plan[, 1]], orderTmp[Plan[, 2]])
    } else if (MateDesign == "circularPlan") {
      Plan <- cbind(1:n, (1:n %% n) + 1)
      parent_comb <- cbind(TargetPop[Plan[, 1]], TargetPop[Plan[, 2]])
    }
  }

  MatePlan <- data.frame(Parent1 = parent_comb[, 1], Parent2 = parent_comb[, 2])
  cat("Number of crosses generated:", nrow(MatePlan), "\n")
  return(MatePlan)
}
