# hierarchicalFusion

Code to implement experiments from [Divide-and-Conquer Monte Carlo Fusion](https://arxiv.org/abs/2110.07265) by Ryan S.Y. Chan, Adam M. Johansen, Murray Pollock and Gareth O. Roberts.

**Note**: package has been renamed to `DCFusion` but the repo is still called `hierarchicalFusion` for now since that is what the current arxiv and submitted version has linked to. This will change when this gets updated.

## Installation

Simply run: `devtools::install_github('rchan26/hierarchicalFusion')`

## Running the experiments

The experiments were ran on [Microsoft Azure](https://azure.microsoft.com/en-gb/) using *Data Science Virtual Machine's (DSVM)* with either 16 core (Section 4) or 64 core machines (Section 5). The code utilises parallel computing (via the base [`parallel`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) package) and by default uses all the cores available on the machine. To change this, modify the `n_cores` variable in the functions which perform the methodology (this is set to `parallel::detectCores()` by default).

* Section 4.1: [varying_rho_replicates.R](https://github.com/rchan26/hierarchicalFusion/blob/main/scripts/bivariate_Gaussian/importance_sampling/varying_rho_replicates.R)
* Section 4.2: [varying_C_experiments_uniG_smc_replicates.R](https://github.com/rchan26/hierarchicalFusion/blob/main/scripts/univariate_Gaussian/importance_sampling/varying_C_experiments_uniG_smc_replicates.R)
* Section 4:3: [separate_modes_smc.R](https://github.com/rchan26/hierarchicalFusion/blob/main/scripts/univariate_Gaussian/importance_sampling/separate_modes_smc.R) and [separate_modes_with_tempering.R](https://github.com/rchan26/hierarchicalFusion/blob/main/scripts/univariate_Gaussian/importance_sampling/separate_modes_with_tempering.R)
* Section 5.1: [logistic_regression/simulated_data/](https://github.com/rchan26/hierarchicalFusion/tree/main/scripts/logistic_regression/simulated_data)
* Section 5.2: [logistic_regression/credit_card/](https://github.com/rchan26/hierarchicalFusion/tree/main/scripts/logistic_regression/credit_card)

## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
