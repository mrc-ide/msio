format_results <- function(
  params,
  seasonality, 
  species_proportions, 
  average_age, 
  interventions,
  nets,
  spraying,
  treatment,
  result,
  outputs,
  seasonality_output
  ) {
  if (seasonality_output == 'daily') {
    rainfall <- get_rainfall(seasonality)
  } else if (seasonality_output == 'coefficients') {
    rainfall <- seasonality
  } else {
    stop('unknown seasonality output')
  }
  if (is.null(result)) {
    return(list())
  }
  list(
    parameters = format_parameters(
      params,
      average_age,
      species_proportions,
      rainfall,
      result$pre
    ),
    timed_parameters = format_timed(
      interventions,
      nets,
      spraying,
      treatment
    ),
    outputs = format_outputs(result$post, outputs)
  )
}

get_spec <- function(params, interventions) {
  n <- names(params)
  n[n == 'init_EIR'] <- 'baseline_EIR'
  list(
    parameters = c(
      n,
      c('arab_prop', 'fun_prop', 'gamb_prop'),
      'average_age',
      seq(365)
    ),
    timed_parameters = interventions,
    outputs = as.character(seq(365))
  )
}

get_EIR <- function(result) {
  (result$EIR_arab + result$EIR_fun + result$EIR_gamb) * 365 / 10000
}

format_parameters <- function(
  params,
  average_age,
  species_proportions,
  rainfall,
  pre
  ) {
  row <- params
  row$init_EIR <- estimate_baseline(pre)
  as.numeric(c(
    row,
    average_age,
    species_proportions,
    rainfall
  ))
}

estimate_baseline <- function(pre) {
  mean(get_EIR(tail(pre, 365)))
}

format_timed <- function(interventions, nets, spraying, treatment) {
  data <- NULL
  if ('nets' %in% interventions) {
    data <- cbind(data, as.numeric(nets))
  }
  if ('spraying' %in% interventions) {
    data <- cbind(data, as.numeric(spraying))
  }
  if ('treatment' %in% interventions) {
    data <- cbind(data, as.numeric(treatment))
  }
  data
}

format_outputs <- function(post, outputs) {
  year <- 365
  row_year <- floor((post$timestep - 1) / year)
  data <- list()
  if ('prev' %in% outputs) {
    prev <- post$n_detect_730_3650 / post$n_730_3650
    data[[length(data) + 1]] <- summarise_daily(prev, row_year, post$timestep)
  }
  if ('prev_6_59' %in% outputs) {
    prev <- post$n_detect_180_1799 / post$n_180_1799
    data[[length(data) + 1]] <- summarise_daily(prev, row_year, post$timestep)
  }
  if ('inc' %in% outputs) {
    inc <- post$n_inc_clinical_0_36500 / post$n_0_36500
    inc[is.na(inc)] <- 0
    data[[length(data) + 1]] <- summarise_daily(inc, row_year, post$timestep)
  }
  if ('eir' %in% outputs) {
    eir <- get_EIR(post)
    data[[length(data) + 1]] <- summarise_daily(eir, row_year, post$timestep)
  }
  do.call(cbind, data)
}

summarise_daily <- function(metric, year, timestep) {
  avg_prof <- aggregate(
    metric,
    by = list(year = year),
    FUN = as.numeric
  )
  avg_prof[order(avg_prof$year), 'x']
}

get_rainfall <- function(seas_row) {
  vapply(
    1:365,
    function(t) malariasimulation:::rainfall(
      t,
      g0 = seas_row$seasonal_a0,
      g = c(seas_row$seasonal_a1, seas_row$seasonal_a2, seas_row$seasonal_a3),
      h = c(seas_row$seasonal_b1, seas_row$seasonal_b2, seas_row$seasonal_b3),
      floor = 0.001
    ),
    numeric(1)
  )
}
