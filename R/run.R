#' @title Run simulations from data
#' @param node an arbitrary number, useful for tracking multiple executions
#' @param warmup number of years to warm up for
#' @param seed random seed
#' @param datadir directory for global data to sample from
#' @param batch_size number of runs per batch
#' @param n_batches number of batches
#' @param outdir directory to save outputs
#' @export
run_simulations_from_data <- function(
  node,
  warmup,
  seed,
  datadir,
  batch_size,
  n_batches,
  outdir
  ) {
  n <- n_batches * batch_size
  set.seed(seed)
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
    seasonality,
    nets,
    spraying,
    treatment,
    warmup,
    batch_size,
    n_batches,
    outdir
  )
}

run_synthetic_simulations <- function(
  node,
  warmup,
  seed,
  n_years,
  batch_size,
  n_batches,
  outdir
  ) {
  set.seed(seed)
  seasonality <- synthetic_seasonality(n)
  nets <- synthetic_nets(n)
  spraying <- synthetic_spraying(n)
  treatment <- synthetic_tx(n)
  run_simulations(
    node,
    seasonality,
    nets,
    spraying,
    treatment,
    warmup,
    batch_size,
    n_batches,
    outdir
  )
}

run_simulations <- function(
  node,
  seasonality,
  nets,
  spraying,
  treatment,
  warmup,
  batch_size,
  n_batches,
  outdir
  ) {
  print(paste0('beginning node ', node))
  n <- n_batches * batch_size
  params <- sample_params(n)
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
            params$basic[i,],
            seasonality[i,], 
            params$mosquito_params[[i]],
            params$mosquito_proportions[i,],
            nets[i,],
            spraying[i,],
            treatment[i,],
            warmup
          )
        }
      )

      output <- lapply(
        seq_along(results),
        function(i) {
          format_results(
            params$basic[i,],
            seasonality[i,], 
            params$mosquito_params[[i]],
            params$mosquito_proportions[i,],
            nets[i,],
            spraying[i,],
            treatment[i,],
            warmup,
            results[[i]]
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
  mosquito_params,
  mosquito_proportions,
  nets,
  spraying,
  treatment,
  warmup
  ) {
  year <- 365
  row <- params
  seas_row <- seasonality
  parameters <- malariasimulation::get_parameters(list(
    human_population = 10000,
    individual_mosquitoes = FALSE,
    model_seasonality = TRUE,
    g0 = seas_row$seasonal_a0,
    g = c(seas_row$seasonal_a1, seas_row$seasonal_a2, seas_row$seasonal_a3),
    h = c(seas_row$seasonal_b1, seas_row$seasonal_b2, seas_row$seasonal_b3),
    average_age = row$average_age,
    sigma_squared = row$sigma_squared,
    du = row$du,
    kb = row$kb,
    ub = row$ub,
    uc = row$uc,
    ud = row$ud,
    kc = row$kc,
    b0 = row$b0,
    b1 = row$b1,
    ib0 = row$ib0,
    ic0 = row$ic0,
    ct = row$ct,
    cd = row$cd,
    #ca = row$ca, #TODO: as above
    cu = row$cu,
    prevalence_rendering_min_ages = c(2 * year),
    prevalence_rendering_max_ages = c(10 * year),
    clinical_incidence_rendering_min_ages = 0,
    clinical_incidence_rendering_max_ages = 100 * year
  ))

  parameters <- malariasimulation::set_species(
    parameters,
    mosquito_params,
    proportions = mosquito_proportions
  )

  parameters <- malariasimulation::set_drugs(
    parameters,
    list(
      malariasimulation::DHA_PQP_params,
      malariasimulation::AL_params,
      malariasimulation::SP_AQ_params
    )
  )

  parameters <- malariasimulation::set_equilibrium(parameters, row$init_EIR)
  
  period <- length(nets)
  one_round_timesteps <- seq(0, period - 1) * year

  # bednets
  parameters <- malariasimulation::set_bednets(
    parameters,
    timesteps = one_round_timesteps + (warmup * year),
    coverages = as.numeric(nets),
    retention = 5 * year,
    dn0 = matrix(.533, nrow=period, ncol=3),
    rn = matrix(.56, nrow=period, ncol=3),
    rnm = matrix(.24, nrow=period, ncol=3),
    gamman = rep(2.64 * year, period)
  )

  # spraying
  parameters <- malariasimulation::set_spraying(
    parameters,
    timesteps = one_round_timesteps + (warmup * year),
    coverages = as.numeric(spraying),
    ls_theta = matrix(2.025, nrow=period, ncol=3),
    ls_gamma = matrix(-0.009, nrow=period, ncol=3),
    ks_theta = matrix(-2.222, nrow=period, ncol=3),
    ks_gamma = matrix(0.008, nrow=period, ncol=3),
    ms_theta = matrix(-1.232, nrow=period, ncol=3),
    ms_gamma = matrix(-0.009, nrow=period, ncol=3)
  )

  # tx
  parameters <- malariasimulation::set_clinical_treatment(
    parameters,
    drug = 2,
    timesteps = one_round_timesteps + warmup * year,
    coverages = as.numeric(treatment)
  )

  malariasimulation::run_simulation(
    (period + warmup) * year,
    parameters = parameters
  )[c(
    'n_inc_clinical_0_36500',
    'n_0_36500',
    'n_detect_730_3650',
    'n_730_3650',
    'EIR_gamb',
    'EIR_arab',
    'EIR_fun'
  )]
}
