all_outputs <- c('prev', 'inc', 'eir')

format_results <- function(
  params,
  seasonality, 
  interventions,
  nets,
  spraying,
  treatment,
  warmup,
  result,
  outputs,
  aggregation
  ) {
  rainfall <- get_rainfall(seasonality)
  list(
    parameters = format_parameters(
      params,
      rainfall,
      warmup,
      result
    ),
    timed_parameters = format_timed(interventions, nets, spraying, treatment),
    outputs = format_outputs(result, warmup, outputs, aggregation)
  )
}

get_spec <- function(params, interventions) {
  n <- names(params)
  n[n == 'init_EIR'] <- 'baseline_EIR'
  list(
    parameters = c(
      n,
      seq(365)
    ),
    timed_parameters = interventions,
    outputs = as.character(seq(365))
  )
}

get_EIR <- function(result) result$EIR_All * 365 / 10000

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

format_timed <- function(interventions, nets, spraying, treatment) {
  data <- NULL
  if ('nets' %in% interventions) {
    data <- c(data, as.numeric(nets))
  }
  if ('spraying' %in% interventions) {
    data <- c(data, as.numeric(spraying))
  }
  if ('treatment' %in% interventions) {
    data <- c(data, as.numeric(treatment))
  }
  matrix(
    data,
    ncol = length(interventions),
    nrow = length(nets)
  )
}

format_outputs <- function(result, warmup, outputs, aggregation) {
  year <- 365
  result <- result[result$timestep > (warmup * year), ]
  row_year <- floor((result$timestep - 1) / year)
  if (aggregation == 'yearly') {
    summarise <- function(x) summarise_yearly(x, row_year, result$repetition)
  } else {
    summarise <- function(x) summarise_daily(x, row_year, result$timestep)
  }
  data <- list()
  if ('prev' %in% outputs) {
    prev <- result$n_detect_730_3650 / result$n_730_3650
    data[[length(data) + 1]] <- summarise(prev)
  }
  if ('inc' %in% outputs) {
    inc <- result$n_inc_clinical_0_36500 / result$n_0_36500
    inc[is.na(inc)] <- 0
    data[[length(data) + 1]] <- summarise(inc)
  }
  if ('eir' %in% outputs) {
    eir <- get_EIR(result)
    data[[length(data) + 1]] <- summarise(eir)
  }
  do.call(cbind, data)
}

#' @importFrom stats aggregate sd
summarise_yearly <- function(metric, year, repetition) {
  avg_year <- aggregate(
    metric,
    by = list(year = year, r = repetition),
    FUN = mean
  )
  avg_rep <- aggregate(
    avg_year$x,
    by = list(year = avg_year$year),
    FUN = mean
  )
  sd_rep <- aggregate(
    avg_year$x,
    by = list(year = avg_year$year),
    FUN = sd
  )
  matrix(c(avg_rep$x, sd_rep$x), ncol = 2, nrow = length(avg_rep$x))
}

#' @importFrom stats aggregate sd
summarise_daily <- function(metric, year, timestep) {
  avg_rep <- aggregate(
    metric,
    by = list(year = year, r = timestep),
    FUN = mean
  )
  avg_prof <- aggregate(
    avg_rep$x,
    by = list(year = avg_rep$year),
    FUN = as.numeric
  )
  avg_prof$x
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
