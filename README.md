# COVID-19 Testing

Comparison of different sampling and estimation strategies for COVID-19 prevalence

## Simulations

- src/R/strs_sim.R contains R code to run the stratified random sampling simulations. It calls either src/stan/strs_mrp.stan or src/stan/strs_pool_mrp.stan to fit the models

- src/R/mcrs_sim.R contains R code to run the multistage cluster random sampling simulations. It calls either src/stan/mcrs_mrp.stan or src/stan/mcrs_pool_mrp.stan to fit the models

- src/R/srs_sim.R contains R code to run the simple random sampling simulations. It calls either src/stan/srs_mrp.stan or src/stan/srs_pool_mrp.stan to fit the models