
sample_params <- function(n) {
  r <- lhs::randomLHS(n, 26)
  
  mosquito_proportions <- matrix(runif(n * 3), nrow = n, ncol = 3)
  mosquito_proportions <- mosquito_proportions / rowSums(mosquito_proportions)

  list(
    basic = sample_basic_params(r, n),
    mosquito_params = sample_mosquito_params(r, n),
    mosquito_proportions = mosquito_proportions
  )
}

sample_basic_params <- function(r, n) {
  db <- 10
  dc <- 30
  a0 <- 2920
  rho <- 0.85
  params <- data.frame(
    average_age = qunif(r[,1], min=20 * 365, max=40 * 365),
    init_EIR = qunif(r[,2], min=0, max=1000),
    sigma_squared = qunif(r[,3], min=1, max=3),
    du = qunif(r[,4], min=30, max=100),
    kb = qunif(r[,5], min=0.01, max=10),
    ub = qunif(r[,6], min=1, max=1000),
    uc = qunif(r[,7], min=1, max=1000), # not documented
    ud = qunif(r[,8], min=1, max=1000), # not documented
    kc = qunif(r[,9], min=0.01, max=10),
    b0 = qunif(r[,10], min=0.01, max=0.99),
    ct = qunif(r[,11], min=0, max=1),
    cd = qunif(r[,12], min=0, max=1),
    #ca = runif(n, 0, 1), #TODO: how is this changed?
    cu = qunif(r[,13], min=0, max=1)
  )

  daily_EIR <- params$init_EIR / 365

  params$ib0 <- daily_EIR * db / (daily_EIR * params$ub + 1) * vapply(
    qunif(r[,14], min=1, max=30) * 365,
    function(a) immunity_scale(a, db, a0, rho),
    numeric(1)
  )

  params$ic0 <- daily_EIR * dc * vapply(
    qunif(r[,15], min=1, max=30) * 365,
    function(a) immunity_scale(a, dc, a0, rho),
    numeric(1)
  )

  params$b1 <- qunif(r[,16], min=0, max=1) * params$b0

  params
}

sample_mosquito_params <- function(r, n) {
  species <- c('gamb', 'arab', 'fun')
  lapply(
    seq(n),
    function(i) {
      lapply(
        seq_along(species),
        function(j) {
          list(
            species = species[[j]],
            mum = 0.1253333,
            blood_meal_rates = 1 / 3,
            foraging_time = 0.69,
            Q0 = qunif(r[i, 17 + j], min=0, max=1),
            phi_indoors = qunif(r[i, 20 + j], min=0, max=1),
            phi_bednets = qunif(r[i, 23 + j], min=0, max=1)
          )
        }
      )
    }
  )
}

sample_df <- function(df, n) {
  df[sample(nrow(df), n, replace = TRUE), ]
}

sample_intervention <- function(df, n) {
  df <- df[,!(names(df) %in% c('Name', 'Gaul_Code'))]
  df <- df[!duplicated(df),]
  df <- df[sample(nrow(df), n, replace = TRUE),]
  df
}

immunity_scale <- function(a, d, a0, rho) {
  (1 + rho / (d / a0 - 1) * exp(-a /a0) - (1 + rho/(d / a0 + 1)) * exp(-a / d))
}
