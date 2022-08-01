params_from_sample <- function(sample, warmup, human_population, interventions) {
  year <- 365
  month <- 30

  row <- static_params_from_sample(sample$params)
  seas_row <- sample$seasonality
  parameters <- malariasimulation::get_parameters(
    c(
      list(
        human_population = human_population,
        average_age = sample$params$average_age,
        individual_mosquitoes = FALSE,
        model_seasonality = TRUE,
        g0 = seas_row$seasonal_a0,
        g = c(seas_row$seasonal_a1, seas_row$seasonal_a2, seas_row$seasonal_a3),
        h = c(seas_row$seasonal_b1, seas_row$seasonal_b2, seas_row$seasonal_b3),
        prevalence_rendering_min_ages = c(6 * month, 2 * year),
        prevalence_rendering_max_ages = c(60 * month - 1, 10 * year),
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
    proportions = sample$species_proportions / sum(sample$species_proportions)
  )

  parameters <- malariasimulation::set_equilibrium(
    parameters,
    sample$params$init_EIR
  )

  period <- length(sample$nets)
  one_round_timesteps <- seq(0, period - 1) * year

  if ('nets' %in% interventions) {
    # bednets
    parameters <- malariasimulation::set_bednets(
      parameters,
      timesteps = one_round_timesteps + warmup * year,
      coverages = as.numeric(sample$nets),
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
      coverages = as.numeric(sample$spraying),
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
      coverages = as.numeric(sample$treatment)
    )
  }

  parameters
}

static_params_from_sample <- function(params) {
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

create_samples <- function(
  n_years,
  paramset = basic_params,
  sites = NULL,
  interventions = 'nets',
  synthetic_intervention_method='lhs',
  n
  ) {
  if (is.null(sites)) {
    r <- lhs::randomLHS(n, length(paramset) + 11)
    seasonality <- synthetic_seasonality(r[,seq(7), drop=FALSE])
    species_proportions <- synthetic_species(r[,c(8, 9, 10), drop=FALSE])
    params <- sample_params(n, paramset, r[seq_along(paramset) + 11,drop=FALSE])
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
  }

  if (synthetic_intervention_method == 'lhs') {
    r <- lhs::randomLHS(n * n_years, length(interventions))
    i <- 1
    if ('nets' %in% interventions) {
      nets <- synthetic_nets_lhs(n, n_years, r[, i])
      i <- i + 1
    }
    if ('spraying' %in% interventions) {
      spraying <- synthetic_spraying_lhs(n, n_years, r[, i])
      i <- i + 1
    }
    if ('treatment' %in% interventions) {
      treatment <- synthetic_tx_lhs(n, n_years, r[, i])
      i <- i + 1
    }
  } else if (synthetic_intervention_method == 'bgw') {
    nets <- synthetic_nets(n, n_years)
    spraying <- synthetic_spraying(n, n_years)
    treatment <- synthetic_tx(n, n_years)
  } else if (synthetic_intervention_method == 'historic') {
    datadir <- system.file('default', package='msio')
    nets <- sample_intervention(read.csv(file.path(datadir, 'nets.csv')), n)
    spraying <- sample_intervention(
      read.csv(file.path(datadir, 'spraying.csv')),
      n
    )
    treatment <- sample_intervention(
      read.csv(file.path(datadir, 'treatment.csv')),
      n
    )
  }

  lapply(
    seq(n),
    function(i) {
      x <- list(
        params = params[i,,drop=FALSE],
        seasonality = seasonality[i,],
        species_proportions = species_proportions[i,]
      )
      if ('nets' %in% interventions) {
        x$nets = nets[i,]
      }
      if ('spraying' %in% interventions) {
        x$spraying = spraying[i,]
      }
      if ('treatment' %in% interventions) {
        x$treatment = treatment[i,]
      }
      x
    }
  )
}

