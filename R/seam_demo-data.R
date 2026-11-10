#' Simulated expression and trait data for two groups and two modules
#'
#' A reproducible example with 48 samples per group, 24 variables, two
#' modules (`blue` and `red`), and two continuous traits. Different subsets
#' of variables dominate the module signals in groups A and B; the first
#' trait has an opposite relationship to the blue module across groups.
#' Data are simulated and have no biological interpretation.
#'
#' @format A list with four components:
#' \describe{
#'   \item{expression_data}{Named list of two 48-by-24 numeric matrices.}
#'   \item{trait_data}{Named list of two 48-by-2 numeric data frames.}
#'   \item{module_assignment}{Character vector of 24 module labels.}
#'   \item{labels}{Character vector `c("A", "B")`.}
#' }
#' @source Simulated with a fixed random seed; see `data-raw/seam_demo.R`.
#' @keywords datasets
"seam_demo"
