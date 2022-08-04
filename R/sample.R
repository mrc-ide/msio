#' @title basic model parameters
#' @description parameters for simple global runs
#' @export
basic_params <- list(
  init_EIR = list(min=0, max=500)
)

#' @title extended parameters for inference
#' @description parameters for inference on:
#'
#' * heterogeneity
#' * infectivity towards mosquitoes
#' * duration of subpatent infection
#' @export
extended_params <- c(
  basic_params,
  list(
    sigma_squared = list(min=1, max=3),
    du = list(min=30, max=100),
    ct = list(min=0, max=1),
    cd = list(min=0, max=1),
    gamma1 = list(min=0.01, max=10),
    cu = list(min=0, max=1)
  )
)

#' @title all parameters for inference
#' @description all inference parameters including immunity
#' @export
all_params <- c(
  extended_params,
  list(
    kb = list(min=0.01, max=10),
    ub = list(min=1, max=10),
    uc = list(min=1, max=10),
    ud = list(min=1, max=10),
    kc = list(min=0.01, max=10),
    b0 = list(min=0.01, max=0.99),
    b1_prop = list(min=0, max=1),
    ib0 = list(min=1, max=100),
    ic0 = list(min=1, max=100)
  )
)

all_interventions <- c(
  'nets',
  'spraying',
  'treatment'
)

#' @description sample a dataframe of params for a paramset
#' @importFrom stats qunif
#' @noRd
sample_params <- function(n, paramset, r) {
  r <- lhs::randomLHS(n, length(paramset))
  cols <- lapply(
    seq_along(paramset),
    function(i) {
      qunif(r[,i], min=paramset[[i]]$min, max=paramset[[i]]$max)
    }
  )
  names(cols) <- names(paramset)
  data.frame(cols)
}

sample_df <- function(df, n) {
  df[sample(nrow(df), n, replace = TRUE), ]
}

sample_intervention <- function(df, n) {
  df <- df[,!(names(df) %in% c('Name', 'Gaul_Code'))]
  df <- df[!duplicated(df),]
  df <- df[sample(nrow(df), n, replace = TRUE),]
  df
}
