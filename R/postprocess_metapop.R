all_outputs <- c('prev', 'inc', 'eir', 'infectivity')

format_metapop_results <- function(
  sample_i,
  interventions,
  warmup,
  result,
  outputs,
  aggregation,
  mixing_matrix,
  mixing_index
  ) {
  formatted <- format_results(
    sample_i,
    interventions,
    warmup,
    result,
    outputs,
    aggregation
  )
  formatted$notes$mixing_matrix <- mixing_matrix
  formatted$notes$mixing_index <- mixing_index
  formatted
}
