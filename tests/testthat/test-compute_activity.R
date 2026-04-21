library(testthat)

test_that("compute_activity assigns correct activity labels and column names", {
  
  spca_data <- list(
    PCs = list(matrix(1), matrix(1)), # two groups, one module
    vs = list(
      list(M1 = c(var1 = 1, var2 = 0)), # weights for module 1
      list(M1 = c(var1 = 0, var2 = 1)) # weights for module 2
    ),
    module_matrices = list(
      list(M1 = matrix( # module 1
        1:4,
        nrow = 2,
        dimnames = list(c("var1", "var2"), c("var1", "var2"))
      )),
      list(M1 = matrix( # module 2
        2:5,
        nrow = 2,
        dimnames = list(c("var1", "var2"), c("var1", "var2"))
      ))
    )
  )
  
  result <- compute_activity(
    spca_data = spca_data,
    return_information = FALSE,
    labels = c("A", "B")
  )
  
  expect_s3_class(result, "data.frame")
  
  expect_identical(
    colnames(result),
    c("variable_id", "module", "w_A", "w_B", "activity")
  )

  expect_setequal(result$variable_id, c("var1", "var2"))
  
  expect_equal(
    result$activity[result$variable_id == "var1"],
    "active_A"
  )
  
  expect_equal(
    result$activity[result$variable_id == "var2"],
    "active_B"
  )
})

test_that("compute_activity validates inputs", {
  
  expect_error( # invalid spca data
    compute_activity(
      spca_data = "not_a_list",
      labels = c("A")
    ),
    "Invalid spca_data argument"
  )
  
  expect_error( # wrong length labels
    compute_activity(
      spca_data = list(PCs = list(matrix(1)), vs = list(), module_matrices = list()),
      labels = c("A", "B")
    ),
    "Length of labels does not coincide"
  )
  
  expect_error( # forbidden symbol in labels
    compute_activity(
      spca_data = list(PCs = list(matrix(1)), vs = list(), module_matrices = list()),
      labels = c("A-B")
    ),
    "Do not use - in labels"
  )
})

