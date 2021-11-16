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
    outputs = format_outputs(result, warmup, outputs, aggregation),
    notes = list(warmup_eirs = warmup_eirs(result, warmup))
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
  row <- params
  row$init_EIR <- estimate_baseline(result, warmup)
  as.numeric(c(
    row,
    rainfall
  ))
}

estimate_baseline <- function(result, warmup) {
  year <- 365
  period <- (
    result$timestep >= (warmup - 1) * year
  ) & (result$timestep < warmup * year)
  period_df <- result[period,]
  mean(
    aggregate(
      get_EIR(period_df),
      by=list(rep = period_df$repetition),
      FUN=mean
    )$x
  )
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
  result <- result[result$timestep > warmup * year, ]
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
  avg_rep <- avg_rep[order(avg_rep$year), 'x']
  avg_sd <- avg_sd[order(avg_sd$year), 'x']
  matrix(c(avg_rep, sd_rep), ncol = 2, nrow = length(avg_rep))
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
  avg_prof[order(avg_prof$year), 'x']
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

warmup_eirs <- function(results, warmup) {
  year <- 365
  results <- results[results$timestep < warmup * year,]
  result_year <- floor(results$timestep / year)
  per_run <- aggregate(
    get_EIR(results),
    by=list(year = result_year, rep = results$repetition),
    FUN=mean
  )
  per_year <- aggregate(
    per_run$x,
    by=list(year = per_run$year),
    FUN=mean
  )
  as.numeric(per_year[order(per_year$year),'x'])
}
