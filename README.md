# seam

**seam** is an R package for exploring and comparing sparse principal component analysis (SPCA) across groups. It supports module-wise analysis, comparison of variable loadings and sparsity, associations with sample traits, and visual summaries of group-specific patterns.

The package uses `elasticnet::spca()` to estimate sparse components. When a module assignment is supplied, SPCA is fitted separately within each module and group. Optional bootstrap resampling can be used to examine how consistently variables receive nonzero loadings.

## Installation

```r
install.packages("seam")
```

## Quick start

The example below uses simulated data with two groups and two predefined modules. Rows represent samples; columns represent the same variables in each group.

```r
library(seam)

set.seed(42)
modules = rep(c("blue", "red"), each = 8)
genes = paste0("gene", seq_len(16))

expression_data = lapply(seq_len(2), function(i) {
  x = matrix(rnorm(30 * 16), nrow = 30, ncol = 16,
             dimnames = list(NULL, genes))
  x[, 1:8] = x[, 1:8] + rnorm(30)
  x[, 9:16] = x[, 9:16] + rnorm(30)
  x
})

fit = spca_all_modules(
  expression_data = expression_data,
  module_assignment = modules,
  lambda = 1e-5,
  bootstrap = FALSE
)

plot_sparsity(fit, labels = c("Group A", "Group B"))
plot_explained_variance(fit, labels = c("Group A", "Group B"))
```

The `lambda` argument controls sparsity and is passed to the `para` argument of `elasticnet::spca()`. It is distinct from that function's own `lambda` argument.

## Explore differences between groups

The fitted object stores component scores, variable loadings, explained variance, and module information. The following functions provide complementary visualizations:

```r
plot_activity_bar(fit, labels = c("Group A", "Group B"),
                  return_information = FALSE)

insight_weights_lollipop(fit, labels = c("Group A", "Group B"), n = 5)
```

You can also explore the effect of the sparsity penalty:

```r
plot_lambdas(
  expression_data = expression_data,
  module_assignment = modules,
  lambda_grid = c(1e-6, 1e-5, 1e-4),
  labels = c("Group A", "Group B"),
  bootstrap = FALSE
)
```

## Assess selection stability

With `bootstrap = TRUE`, samples are resampled with replacement within each group, and SPCA is refitted. The bootstrap loading matrices can be used to calculate the fraction of resamples in which each variable is selected:

```r
set.seed(42)
fit_boot = spca_all_modules(
  expression_data = expression_data,
  module_assignment = modules,
  lambda = 1e-5,
  bootstrap = TRUE,
  n = 100
)

# Selection frequencies for the first module in Group A
selection_frequency = colMeans(fit_boot$bootstrap$vs[[1]][[1]] != 0)
print(selection_frequency)
```

Selection frequencies describe the stability of the fitted procedure under resampling. They are not p-values, do not establish statistical significance, and do not by themselves demonstrate differences in population-level variable support between groups.

## Trait associations

When sample-level traits are available, `seam` can calculate and visualize module–trait and variable–trait associations using `compute_module_trait_correlations()`, `compute_variable_trait_correlations()`, and their associated plotting functions. See the individual help pages for input requirements and examples.

## Documentation

Use `help(package = "seam")` to browse available functions and `?spca_all_modules` to get started. The package includes examples for its exported functions.

## Citation

If you use `seam` in your research, please cite the package. The underlying SPCA method is implemented through the **elasticnet** R package.
