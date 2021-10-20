format_results <- function(
  params,
  seasonality, 
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
      warmup,
      result
    ),
    timed_parameters = format_timed(nets, spraying, treatment),
    outputs = format_outputs(result, warmup)
  )
}

get_spec <- function(params) {
  n <- names(params)
  n[[1]] <- 'baseline_EIR'
  list(
    parameters = c(
      n,
      seq(365)
    ),
    timed_parameters = c('nets', 'spraing', 'treatment'),
    outputs = as.character(seq(365))
  )
}

get_EIR <- function(result) result$EIR_All

format_parameters <- function(params, rainfall, warmup, result) {
  year <- 365
  row <- params
  row$init_EIR <- NULL
  as.numeric(c(
    mean(get_EIR(result)[seq((warmup - 1) * year, warmup * year)]),
    row,
    rainfall
  ))
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
