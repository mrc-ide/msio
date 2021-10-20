format_results <- function(
  params,
  seasonality, 
  mosquito_params,
  mosquito_proportions,
  nets,
  spraying,
  treatment,
  warmup,
  result
  ) {
  rainfall <- get_rainfall(seasonality)
  list(
    parameters = format_parameters(
      params,
      rainfall,
      mosquito_params,
      mosquito_proportions,
      warmup,
      result
    ),
    timed_parameters = format_timed(nets, spraying, treatment),
    outputs = format_outputs(result, warmup)
  )
}

get_spec <- function() {
  list(
    parameters = c(
      'baseline_EIR',
      'average_age',
      'sigma_squared',
      'du',
      'kb',
      'ub',
      'uc',
      'ud',
      'kc',
      'b0',
      'b1',
      'ib0',
      'ic0',
      'ct',
      'cd',
      'cu',
      'Q0_gamb',
      'phi_indoors_gamb',
      'phi_bednets_gamb',
      'Q0_fun',
      'phi_indoors_fun',
      'phi_bednets_fun',
      'Q0_arab',
      'phi_indoors_arab',
      'phi_bednets_arab',
      'prop_gamb',
      'prop_fun',
      'prop_arab',
      seq(365)
    ),
    timed_parameters = c('nets', 'spraing', 'treatment'),
    outputs = as.character(seq(365))
  )
}

get_EIR <- function(result) {
  colSums(rbind(result$EIR_gamb, result$EIR_fun, result$EIR_arab))
}

format_parameters <- function(params, rainfall, m_params, m_prop, warmup, result) {
  year <- 365
  species_vector <- unlist(lapply(
    m_params,
    function(p) {
      c(
        p$Q0,
        p$phi_indoors,
        p$phi_bednets
      )
    }
  ))
  row <- params
  c(
    mean(get_EIR(result)[seq((warmup - 1) * year, warmup * year)]),
    row$average_age,
    row$sigma_squared,
    row$du,
    row$kb,
    row$ub,
    row$uc,
    row$ud,
    row$kc,
    row$b0,
    row$b1,
    row$ib0,
    row$ic0,
    row$ct,
    row$cd,
    row$cu,
    species_vector,
    m_prop,
    rainfall
  )
}

format_timed <- function(nets, spraying, treatment) {
  matrix(
    c(
      as.numeric(nets),
      as.numeric(spraying),
      as.numeric(treatment)
    ),
    ncol = 3,
    nrow = length(nets)
  )
}

format_outputs <- function(result, warmup) {
  year <- 365
  result <- result[seq(warmup * year + 1, nrow(result)),]
  series <- result$n_detect_730_3650 / result$n_730_3650
  t(matrix(series, nrow=year))
}

get_rainfall <- function(seas_row) {
  vapply(
    1:365,
    function(t) malariasimulation:::rainfall(
      t,
      g0 = seas_row$seasonal_a0,
      g = c(seas_row$seasonal_a1, seas_row$seasonal_a2, seas_row$seasonal_a3),
      h = c(seas_row$seasonal_b1, seas_row$seasonal_b2, seas_row$seasonal_b3)
    ),
    numeric(1)
  )
}

