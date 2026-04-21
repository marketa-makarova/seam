test_that("compute_eigengenes returns correctly shaped eigengene matrices", {
  
  module_matrices <- list( # one group
    list( # two modules
      M1 = matrix( # module 1
        c(1, 2,
          3, 4),
        nrow = 2,
        dimnames = list(c("s1", "s2"), c("v1", "v2"))
      ),
      M2 = matrix( # module 2
        c(2, 1,
          4, 3),
        nrow = 2,
        dimnames = list(c("s1", "s2"), c("v3", "v4"))
      )
    )
  )
  
  result <- compute_eigengenes(
    module_matrices = module_matrices,
    module_assignment =  c(1, 2),
    module_colors = c("blue", "red")
  )
  
  expect_type(result, "list")
  expect_length(result, 1)
  expect_true(is.matrix(result[[1]]))
  
  expect_equal(
    dim(result[[1]]),
    c(2, 2)  # samples / modules
  )
  
  expect_true(all(is.finite(result[[1]])))
  expect_true(is.numeric(result[[1]]))
  
})
