test_that("format_parameters output is the right size", {
  params <- list(
    a = 1,
    b = 2,
    c = 3
  )
  result <- data.frame(EIR_All = rep(4, 5 * 365))
  output <- format_parameters(params, seq(365), 5, result)
  expect_equal(output, c(4, 1, 2, 3, seq(365)))
})
