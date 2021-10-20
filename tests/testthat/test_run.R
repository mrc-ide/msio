test_that('params_from_sample creates valid params', {
  actual <- params_from_sample(
    list(
      average_age = 20 * 365,
      init_EIR = 20,
      Q0 = .5,
      phi_indoors = .5,
      phi_bednets = .5,
      sigma_squared = 1,
      du = 30,
      ct = .5,
      cd = .5,
      gamma1 = 1,
      cu = .5,
      kb = 1,
      ub = 5,
      uc = 5,
      ud = 5,
      kc = 1,
      b0 = .5,
      b1_prop = .5,
      ib0 = 5,
      ic0 = 5
    )
  )
  expected <- list(
    average_age = 20 * 365,
    Q0 = .5,
    phi_indoors = .5,
    phi_bednets = .5,
    sigma_squared = 1,
    du = 30,
    ct = .5,
    cd = .5,
    gamma1 = 1,
    cu = .5,
    kb = 1,
    ub = 5,
    uc = 5,
    ud = 5,
    kc = 1,
    b0 = .5,
    ib0 = 5,
    ic0 = 5,
    b1 = .25
  )
  expect_equal(actual, expected)
})
