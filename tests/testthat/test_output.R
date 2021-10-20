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

test_that("format_output summarises correctly", {
  result <- data.frame(
    timestep = rep(seq(10 * 365), 5),
    repetition = rep(seq(5), each=10 * 365),
    n_detect_730_3650 = rep(50, 5 * 10 * 365),
    n_730_3650 = rep(100, 5 * 10 * 365),
    n_inc_clinical_0_36500 = rep(100, 5 * 10 * 365),
    n_0_36500 = rep(1000, 5 * 10 * 365),
    EIR_All = rep(50, 5 * 10 * 365)
  )
  actual <- format_outputs(result, 5)
  expected <- matrix(
    c(
      rep(.5, 5),
      rep(0, 5),
      rep(.1, 5),
      rep(0, 5),
      rep(50, 5),
      rep(0, 5)
    ),
    nrow = 5,
    ncol = 6
  )
  expect_equal(actual, expected)
})
