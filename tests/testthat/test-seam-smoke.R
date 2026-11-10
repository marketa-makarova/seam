# one small fitted dataset shared by all tests
set.seed(42)

seam_modules = rep(c("blue", "red"), each = 8)
seam_genes = paste0("gene", seq_len(16))
seam_samples = paste0("s", seq_len(30))
seam_expression = lapply(seq_len(2), function(i) {
  x = matrix(rnorm(30 * 16), 30, 16,
             dimnames = list(seam_samples, seam_genes))
  x[, 1:8] = x[, 1:8] + outer(rnorm(30), rep(1, 8))
  x[, 9:16] = x[, 9:16] + outer(rnorm(30), rep(1, 8))
  x
})
seam_labels = c("A", "B")
seam_traits = lapply(seam_expression, function(x) {
  data.frame(trait1 = x[, 1] + rnorm(nrow(x), sd = 0.2),
             trait2 = x[, 9] + rnorm(nrow(x), sd = 0.2),
             row.names = rownames(x))
})
seam_fit = spca_all_modules(seam_expression, seam_modules,
                            lambda = 0.00001, bootstrap = FALSE)
seam_mtc = compute_module_trait_correlations(
  seam_traits, seam_fit, seam_modules, mode = "spca")
# Ensure the plot-selection threshold selects a module/trait pair.
seam_plot_mtc = seam_mtc
seam_plot_mtc$mtc_spca = lapply(seam_mtc$mtc_spca, function(m) {
  m[1, 1] = 0.95
  m
})

# function which checks whether the result is a list of ggplots
seam_is_ggplot_list = function(x) {
  is.list(x) && length(x) > 0 &&
    all(vapply(x, inherits, logical(1), "ggplot"))
}

test_that("spca_all_modules fits and rejects bad input", {
  expect_length(seam_fit$vs, 2)
  expect_equal(seam_fit$module_colors, c("blue", "red"))
  expect_error(spca_all_modules(seam_expression, seam_modules, lambda = -1), "non-negative lambda")
})

test_that("spca_all_modules_grid fits and rejects bad lambdas", {
  fit = spca_all_modules_grid(seam_expression, seam_modules, lambdas = c(0.00001, 0.00001), bootstrap = FALSE)
  expect_length(fit$vs, 2)
  expect_error(spca_all_modules_grid(seam_expression, seam_modules, lambdas = -1, bootstrap = FALSE), "one finite non-negative lambda")
})

test_that("compute_activity returns variable statuses and checks labels", {
  x = compute_activity(seam_fit, labels = seam_labels, return_information = FALSE)
  expect_s3_class(x, "data.frame")
  expect_true(all(c("variable_id", "module", "activity", "w_A", "w_B") %in% names(x)))
  expect_error(compute_activity(seam_fit, labels = "A"), "unique labels")
})

test_that("compute_module_trait_correlations returns matrices and checks mode", {
  x = compute_module_trait_correlations(seam_traits, seam_fit, seam_modules, mode = "spca")
  expect_equal(dim(x$mtc_spca[[1]]), c(2, 2))
  expect_error(compute_module_trait_correlations(seam_traits, seam_fit, seam_modules, mode = "invalid"), "Mode must")
})

test_that("compute_variable_trait_correlations returns columns and checks groups", {
  x = compute_variable_trait_correlations(seam_traits, seam_fit, labels = seam_labels)
  expect_length(x, 2)
  expect_true(all(c("trait1_cor", "trait1_w_cor") %in% names(x[[1]])))
  expect_error(compute_variable_trait_correlations(seam_traits[1], seam_fit, labels = seam_labels), "one trait dataset")
})

test_that("insight_correlations_lollipop returns plots and checks type", {
  x = insight_correlations_lollipop(seam_plot_mtc, seam_fit, labels = seam_labels, interest_threshold = 0.1)
  expect_true(seam_is_ggplot_list(x))
  expect_error(insight_correlations_lollipop(seam_plot_mtc, seam_fit, labels = seam_labels, type = "invalid"), "Type must")
})

test_that("insight_correlations_reconstruct returns plots and checks top n", {
  x = insight_correlations_reconstruct(seam_plot_mtc, seam_fit, labels = seam_labels, top_n_lower = 1, top_n_upper = 4, interest_threshold = 0.1)
  expect_true(seam_is_ggplot_list(x))
  expect_error(insight_correlations_reconstruct(seam_plot_mtc, seam_fit, labels = seam_labels, top_n_lower = 5, top_n_upper = 2), "positive integers")
})

test_that("insight_correlations_sign_changes returns plots and checks threshold", {
  x = insight_correlations_sign_changes(seam_plot_mtc, seam_fit, labels = seam_labels, interest_threshold = 0.1)
  expect_true(seam_is_ggplot_list(x))
  expect_error(insight_correlations_sign_changes(seam_plot_mtc, seam_fit, labels = seam_labels, interest_threshold = 2), "between 0 and 1")
})

test_that("insight_weights_lollipop returns plots and checks group count", {
  x = insight_weights_lollipop(seam_fit, labels = seam_labels, n = 5)
  expect_true(seam_is_ggplot_list(x))
  expect_error(insight_weights_lollipop(seam_fit, labels = "A"), "exactly two groups")
})

test_that("plot_activity_bar returns a plot and checks module names", {
  x = plot_activity_bar(seam_fit, labels = seam_labels, return_information = FALSE)
  expect_s3_class(x, "ggplot")
  expect_error(plot_activity_bar(seam_fit, labels = seam_labels, modules = "unknown"), "Unknown module")
})

test_that("plot_activity_bubble returns a plot and checks scale range", {
  x = plot_activity_bubble(seam_fit, labels = seam_labels, return_information = FALSE)
  expect_s3_class(x, "ggplot")
  expect_error(plot_activity_bubble(seam_fit, labels = seam_labels, scale_range = 2), "scale_range")
})

test_that("plot_explained_variance returns a plot and checks modules", {
  x = plot_explained_variance(seam_fit, labels = seam_labels)
  expect_s3_class(x, "ggplot")
  expect_error(plot_explained_variance(seam_fit, labels = seam_labels, modules = "unknown"), "Unknown module")
})

test_that("plot_sparsity returns a plot and checks modules", {
  x = plot_sparsity(seam_fit, labels = seam_labels)
  expect_s3_class(x, "ggplot")
  expect_error(plot_sparsity(seam_fit, labels = seam_labels, modules = "unknown"), "Unknown module")
})

test_that("plot_module_heatmaps returns plots and checks modules", {
  x = plot_module_heatmaps(seam_fit, labels = seam_labels, return_information = FALSE)
  expect_true(seam_is_ggplot_list(x))
  expect_error(plot_module_heatmaps(seam_fit, labels = seam_labels, modules = "unknown"), "Unknown module")
})

test_that("plot_module_trait_correlations returns plots and checks mode", {
  x = plot_module_trait_correlations(seam_mtc, labels = seam_labels, mode = "spca")
  expect_true(seam_is_ggplot_list(x))
  expect_error(plot_module_trait_correlations(seam_mtc, labels = seam_labels, mode = "invalid"), "Mode must")
})

test_that("plot_lambdas returns plots and checks lambda grid", {
  x = plot_lambdas(seam_expression, seam_modules, lambda_grid = c(0.00001, 0.00002), labels = seam_labels)
  expect_s3_class(x$explained_variance, "ggplot")
  expect_s3_class(x$sparsity, "ggplot")
  expect_error(plot_lambdas(seam_expression, seam_modules, lambda_grid = -1, labels = seam_labels), "finite non-negative")
})
