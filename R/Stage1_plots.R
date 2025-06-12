#' Stage1_plots
#'
#' Generate diagnostic plots for the Stage 1 model fitted to a single trait using the SPATS solver.
#' The function creates a boxplot of residuals, a QQ plot, a heritability plot, and a spatial field map.
#'
#' @param filename Character string specifying the path to the phenotype data file for StageWise.
#' @param trait Character string specifying a single trait to be analyzed in StageWise.
#' @param Fixed Formula or list of fixed effects to be passed to StageWise.
#' @param Random Formula or list of random effects to be passed to StageWise.
#' @return A list containing four ggplot objects:
#' \itemize{
#'   \item \code{boxplot}: Boxplot of model residuals.
#'   \item \code{qqplot}: QQ plot of model residuals.
#'   \item \code{heritability}: Boxplot of broad-sense heritability (H^2) by environment.
#'   \item \code{spatial}: Spatial field map of residuals and trends.
#' }
#'
#' @importFrom ggplot2 ggplot aes stat_boxplot xlab ylab ylim
#' @importFrom StageWise Stage1
#' @export

Stage1_plots <- function(filename, trait, Fixed, Random) {
  effects <- model.terms(Fixed = Fixed, Random = Random)
  model <- StageWise::Stage1(filename = filename, traits = trait, effects = effects,
                             solver = "spats", spline = c("row", "range"))
  plots <- list(
    boxplot = model$resid$boxplot,
    qqplot = model$resid$qqplot,
    heritability = ggplot2::ggplot(data = model$fit, aes(x = env, y = H2)) +
      ggplot2::stat_boxplot(outlier.color = "red") +
      ggplot2::xlab("Location") +
      ggplot2::ylab(expression(paste("Broad-sense ", H^2, " (plot basis)"))) +
      ggplot2::ylim(0, 1),
    spatial = model$resid$spatial
  )
  return(plots)
}
