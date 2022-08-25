# msio

Samples malariasimulation inputs and calculates corresponding outputs for surrogate modelling.

# Installation

From the root of this directory, you can install this package with the following
command:

```
install.packages('.', repos=NULL, type='source')
```

# Usage

You can run a quick simulation with the following R function:

```
msio::run_synthetic_simulations(
  paramset=msio::basic_params,
  outputs='prev',
  sites=NULL,
  batch_size=10,
  n_batches=1,
  reps=1,
  synthetic_intervention_method='lhs',
  human_population=100
)
```

For more usage please consult the function documentation in an R terminal, i.e.
`?msio::run_synthetic_simulations`.
