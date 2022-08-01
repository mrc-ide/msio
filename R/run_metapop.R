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
run_synthetic_metapop_simulations <- function(
  node = 1,
  n_years = 5,
  warmup = 5,
  paramset = basic_params,
  seed = 42,
  sites = NULL,
  interventions = 'nets',
  batch_size = 1,
  n_batches = 1,
  outdir = '.',
  outputs = c('prev', 'eir', 'infectivity'),
  aggregation = 'daily',
  synthetic_intervention_method='lhs',
  human_population = 1e5,
  n_pop = 2
  ) {
  n <- n_batches * batch_size
  set.seed(seed)
  samples <- lapply(
    seq(n_pop),
    function(.) {
      create_samples(
        n_years,
        paramset,
        sites,
        interventions,
        synthetic_intervention_method,
        n
      )
    }
  )

  r <- lhs::randomLHS(n, n_pop ^ 2)
  mixing_matrices <- lapply(
    seq(nrow(r)),
    function(i) {
      m <- matrix(r[i,], nrow = n_pop, ncol = n_pop)
      m / rowSums(m)
    }
  )

  run_metapop_simulations(
    node,
    samples,
    interventions,
    warmup,
    batch_size,
    n_batches,
    outdir,
    outputs,
    aggregation,
    human_population,
    mixing_matrices
  )
}

run_metapop_simulations <- function(
  node,
  samples,
  interventions,
  warmup,
  batch_size,
  n_batches,
  outdir,
  outputs,
  aggregation,
  human_population,
  mixing_matrices
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
          run_metapop_row(
            lapply(seq_along(samples), function(pop) samples[[pop]][[i]]),
            interventions,
            warmup,
            human_population,
            mixing_matrices[[i]]
          )
        }
      )

      output <- NULL

      for (i in seq_along(results)) {
        for (pop in seq_along(samples)) {
          output <- c(
            output,
            list(
              format_metapop_results(
                samples[[pop]][[i]],
                interventions,
                warmup,
                results[[i]][[pop]],
                outputs,
                aggregation,
                mixing_matrices[[i]],
                pop
              )
            )
          )
        }
      }

      jsonlite::write_json(output, outpath, auto_unbox=TRUE, pretty=TRUE)
      print(paste0('node ', node, ' batch ', batch_i, ' completed'))
      print(Sys.time())
      print(Sys.time() - start_time)
    }
  }
}

run_metapop_row <- function(
  samples,
  interventions,
  warmup,
  human_population,
  mixing
  ) {
  year <- 365
  month <- 30
  period <- length(samples[[1]]$nets)
  
  parameters <- lapply(
    samples,
    function(s) params_from_sample(s, warmup, human_population, interventions)
  )

  lapply(
    malariasimulation::run_metapop_simulation(
      (period + warmup) * year,
      parameters,
      NULL,
      mixing
    ),
    function(df) {
      df$repetition <- 1
      df[c(
        'timestep',
        'repetition',
        'n_inc_clinical_0_36500',
        'n_0_36500',
        'n_detect_730_3650',
        'n_730_3650',
        'n_detect_180_1799',
        'n_180_1799',
        'EIR_arab',
        'EIR_fun',
        'EIR_gamb',
        'infectivity'
      )]
    }
  )
}
