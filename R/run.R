#' @title Run simulations from data
#' @param node an arbitrary number, useful for tracking multiple executions
#' @param warmup number of years to warm up for
#' @param paramset a list of parameters from sample.R
#' @param seed random seed
#' @param datadir directory for global data to sample from
#' @param interventions character vector of interventions to include
#' @param batch_size number of runs per batch
#' @param n_batches number of batches
#' @param outdir directory to save outputs
#' @param outputs character vector of outputs to include
#' @param aggregation type of aggregation for outputs, either 'daily' or 'yearly'
#' @importFrom utils read.csv
#' @export
run_simulations_from_data <- function(
  node = 1,
  warmup = 5,
  reps = 10,
  paramset = basic_params,
  seed = 42,
  datadir = NULL,
  interventions = 'nets',
  batch_size = 1,
  n_batches = 1,
  outdir = '.',
  outputs = 'prev',
  aggregation = 'daily'
  ) {
  stopifnot(length(interventions) > 0)
  n <- n_batches * batch_size
  set.seed(seed)
  r <- lhs::randomLHS(n, length(paramset) + 3)
  species_proportions <- synthetic_species(r[,c(1, 2)])
  demography <- synthetic_demography(r[,3])
  params <- sample_params(n, paramset, r[,seq_along(paramset) + 3])
  if (is.null(datadir)) {
    datadir <- system.file('default', package='msio')
  }
  seasonality <- sample_df(
    read.csv(file.path(datadir, 'seasonality.csv')),
    n
  )
  nets <- sample_intervention(read.csv(file.path(datadir, 'nets.csv')), n)
  spraying <- sample_intervention(read.csv(file.path(datadir, 'spraying.csv')), n)
  treatment <- sample_intervention(read.csv(file.path(datadir, 'treatment.csv')), n)

  run_simulations(
    node,
    params,
    seasonality,
    species_proportions,
    demography,
    interventions,
    nets,
    spraying,
    treatment,
    warmup,
    reps,
    batch_size,
    n_batches,
    outdir,
    outputs,
    aggregation
  )
}

#' @title Run synthetic simulations
#' @param node an arbitrary number, useful for tracking multiple executions
#' @param n_years number of years to simulate
#' @param warmup number of years to warm up for
#' @param paramset a list of parameters from sample.R
#' @param seed random seed
#' @param interventions vector of interventions to include
#' @param batch_size number of runs per batch
#' @param n_batches number of batches
#' @param outdir directory to save outputs
#' @param outputs character vector of outputs to include
#' @param aggregation type of aggregation for outputs, either 'daily' or 'yearly'
#' @export
run_synthetic_simulations <- function(
  node = 1,
  n_years = 5,
  warmup = 5,
  reps = 10,
  paramset = basic_params,
  seed = 42,
  sites = NULL,
  interventions = 'nets',
  batch_size = 1,
  n_batches = 1,
  outdir = '.',
  outputs = 'prev',
  aggregation = 'daily'
  ) {
  n <- n_batches * batch_size
  set.seed(seed)
  if (is.null(sites)) {
    r <- lhs::randomLHS(n, length(paramset) + 10)
    seasonality <- synthetic_seasonality(r[,seq(7)])
    species_proportions <- synthetic_species(r[,c(8, 9)])
    demography <- synthetic_demography(r[,10])
    params <- sample_params(n, paramset, r[seq_along(paramset) + 9])
  } else {
    r <- lhs::randomLHS(n, length(paramset))
    params <- sample_params(n, paramset, r)
    site_df <- read.csv(file.path(sites))
    site_df <- sample_df(site_df, n)
    seasonality <- site_df[,c(
      'seasonal_a0',
      'seasonal_a1',
      'seasonal_a2',
      'seasonal_a3',
      'seasonal_b1',
      'seasonal_b2',
      'seasonal_b3'
    )]
    species_proportions <- site_df[,c(
      'arab_prop',
      'fun_prop',
      'gamb_prop'
    )]
    demography <- site_df$average_age
  }

  nets <- synthetic_nets(n, n_years)
  spraying <- synthetic_spraying(n, n_years)
  treatment <- synthetic_tx(n, n_years)

  run_simulations(
    node,
    params,
    seasonality,
    species_proportions,
    demography,
    interventions,
    nets,
    spraying,
    treatment,
    warmup,
    reps,
    batch_size,
    n_batches,
    outdir,
    outputs,
    aggregation
  )
}

run_simulations <- function(
  node,
  params,
  seasonality,
  species_proportions,
  demography,
  interventions,
  nets,
  spraying,
  treatment,
  warmup,
  reps,
  batch_size,
  n_batches,
  outdir,
  outputs,
  aggregation
  ) {
  print(paste0('beginning node ', node))
  n <- n_batches * batch_size
  batches <- split(
    seq(n),
    (seq(n)-1) %/% batch_size
  )
  for (batch_i in seq_along(batches)) {
    outpath <- file.path(
      outdir,
      paste0('realisation_', node, '_batch_', batch_i, '.json')
    )
    if (!file.exists(outpath)) {
      start_time <- Sys.time()
      print(paste0('node ', node, ' batch ', batch_i, ' starting'))
      # do the work
      results <- lapply(
        batches[[batch_i]],
        function(i) {
          run_row(
            params[i,,drop=FALSE],
            seasonality[i,], 
            species_proportions[i,], 
            demography[[i]], 
            interventions,
            nets[i,],
            spraying[i,],
            treatment[i,],
            warmup,
            reps
          )
        }
      )

      output <- lapply(
        seq_along(results),
        function(i) {
          format_results(
            params[i,,drop=FALSE],
            seasonality[i,], 
            species_proportions[i,], 
            demography[[i]], 
            interventions,
            nets[i,],
            spraying[i,],
            treatment[i,],
            warmup,
            results[[i]],
            outputs,
            aggregation
          )
        }
      )

      jsonlite::write_json(output, outpath, auto_unbox=TRUE, pretty=TRUE)
      print(paste0('node ', node, ' batch ', batch_i, ' completed'))
      print(Sys.time())
      print(Sys.time() - start_time)
    }
  }
}

run_row <- function(
  params,
  seasonality, 
  species_proportions,
  demography,
  interventions,
  nets,
  spraying,
  treatment,
  warmup,
  reps
  ) {
  year <- 365
  row <- params_from_sample(params)
  
  seas_row <- seasonality
  parameters <- malariasimulation::get_parameters(
    c(
      list(
        human_population = 10000,
        average_age = demography,
        individual_mosquitoes = FALSE,
        model_seasonality = TRUE,
        g0 = seas_row$seasonal_a0,
        g = c(seas_row$seasonal_a1, seas_row$seasonal_a2, seas_row$seasonal_a3),
        h = c(seas_row$seasonal_b1, seas_row$seasonal_b2, seas_row$seasonal_b3),
        prevalence_rendering_min_ages = c(2 * year),
        prevalence_rendering_max_ages = c(10 * year),
        clinical_incidence_rendering_min_ages = 0,
        clinical_incidence_rendering_max_ages = 100 * year
      ),
      row
    )
  )

  parameters <- malariasimulation::set_drugs(
    parameters,
    list(
      malariasimulation::DHA_PQP_params,
      malariasimulation::AL_params,
      malariasimulation::SP_AQ_params
    )
  )

  parameters <- malariasimulation::set_species(
    parameters,
    species = list(
      malariasimulation::arab_params,
      malariasimulation::fun_params,
      malariasimulation::gamb_params
    ),
    proportions = species_proportions / sum(species_proportions)
  )

  parameters <- malariasimulation::set_equilibrium(parameters, params$init_EIR)
  
  period <- length(nets)
  one_round_timesteps <- seq(0, period - 1) * year

  if ('nets' %in% interventions) {
    # bednets
    parameters <- malariasimulation::set_bednets(
      parameters,
      timesteps = one_round_timesteps + warmup * year,
      coverages = as.numeric(nets),
      retention = 5 * year,
      dn0 = matrix(.533, nrow=period, ncol=3),
      rn = matrix(.56, nrow=period, ncol=3),
      rnm = matrix(.24, nrow=period, ncol=3),
      gamman = rep(2.64 * year, period)
    )
  }

  if ('spraying' %in% interventions) {
    # spraying
    parameters <- malariasimulation::set_spraying(
      parameters,
      timesteps = one_round_timesteps + warmup * year,
      coverages = as.numeric(spraying),
      ls_theta = matrix(2.025, nrow=period, ncol=1),
      ls_gamma = matrix(-0.009, nrow=period, ncol=1),
      ks_theta = matrix(-2.222, nrow=period, ncol=1),
      ks_gamma = matrix(0.008, nrow=period, ncol=1),
      ms_theta = matrix(-1.232, nrow=period, ncol=1),
      ms_gamma = matrix(-0.009, nrow=period, ncol=1)
    )
  }

  if ('treatment' %in% interventions) {
    # tx
    parameters <- malariasimulation::set_clinical_treatment(
      parameters,
      drug = 2,
      timesteps = one_round_timesteps + warmup * year,
      coverages = as.numeric(treatment)
    )
  }

  malariasimulation::run_simulation_with_repetitions(
    (period + warmup) * year,
    reps,
    parameters
  )[c(
    'timestep',
    'repetition',
    'n_inc_clinical_0_36500',
    'n_0_36500',
    'n_detect_730_3650',
    'n_730_3650',
    'EIR_arab',
    'EIR_fun',
    'EIR_gamb'
  )]
}

params_from_sample <- function(params) {
  params$init_EIR <- NULL
  if (!is.null(params$b1_prop)) {
    params$b1 <- params$b0 * params$b1_prop
    params$b1_prop <- NULL
  }
  if (length(params) == 1) {
    return(NULL)
  }
  params
}
