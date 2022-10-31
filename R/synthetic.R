synthetic_seasonality <- function(r) {
  data.frame(list(
    g0 = qunif(r[,1], min=-10, max=10),
    g1 = qunif(r[,2], min=-10, max=10),
    g2 = qunif(r[,3], min=-10, max=10),
    g3 = qunif(r[,4], min=-10, max=10),
    h1 = qunif(r[,5], min=-10, max=10),
    h2 = qunif(r[,6], min=-10, max=10),
    h3 = qunif(r[,7], min=-10, max=10)
  ))
}

synthetic_nets <- function(n, n_years) {
  start <- matrix(0, nrow=n, ncol=1)
  bounded_gauss_random_walk(start, n_years, .1)
}

synthetic_spraying <- function(n, n_years) {
  start <- matrix(0, nrow=n, ncol=1)
  bounded_gauss_random_walk(start, n_years, .1)
}

synthetic_tx <- function(n, n_years) {
  start <- matrix(0, nrow=n, ncol=1)
  bounded_gauss_random_walk(start, n_years, .1)
}

synthetic_nets_lhs <- function(n, n_years, r) {
  matrix(r, nrow=n, ncol=n_years)
}

synthetic_spraying_lhs <- function(n, n_years, r) {
  matrix(r, nrow=n, ncol=n_years)
}

synthetic_tx_lhs <- function(n, n_years, r) {
  matrix(r, nrow=n, ncol=n_years)
}

bounded_gauss_random_walk <- function(start, n, sigma) {
  walk <- start
  for (i in seq(2, n)) {
    step <- rnorm(nrow(start), 0, sigma)
    position <- walk[,i - 1] + step
    position <- pmin(pmax(position, 0), 1)
    walk <- cbind(walk, position)
  }
  walk
}

synthetic_species <- function(r) {
  r <- r / rowSums(r)
  data.frame(list(
    arab_prop = r[,1],
    fun_prop = r[,2],
    gamb_prop = r[,3]
  ))
}

synthetic_demography <- function(r) {
  qunif(r[,1], min=15*365, max=60*365)
}
