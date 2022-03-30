library(HMCBRR)
library(readxl)

seed <- 2022
nsamples <- 30000
warmup <- 20000
sigma <- 1
prior_means <- rep(0, 5)
prior_variances <- rep(10, 5)
n_cores <- parallel::detectCores()

##### Loading in Data #####

# Features consist of hourly average ambient variables
# - Temperature (T) in the range 1.81°C and 37.11°C,
# - Ambient Pressure (AP) in the range 992.89-1033.30 milibar,
# - Relative Humidity (RH) in the range 25.56% to 100.16%
# - Exhaust Vacuum (V) in teh range 25.36-81.56 cm Hg
# - Net hourly electrical energy output (EP) 420.26-495.76 MW

load_pp_data <- function(file, standardise_variables = TRUE) {
  original_data <- as.data.frame(readxl::read_xlsx(file))
  colnames(original_data) <- c('T', 'V', 'AP', 'RH', 'EP')
  if (standardise_variables) {
    X <- subset(original_data, select = -c(EP))
    variable_means <- rep(NA, ncol(X))
    variable_sds <- rep(NA, ncol(X))
    for (col in 1:ncol(X)) {
      variable_means[col] <- mean(X[,col])
      variable_sds[col] <- sd(X[,col])
      X[,col] <- (X[,col]-variable_means[col])/variable_sds[col]
    }
    design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
    colnames(design_mat)[1] <- 'intercept'
    return(list('data' = cbind('EP' = original_data$EP, X),
                'y' = original_data$EP,
                'X' = design_mat,
                'variable_means' = variable_means,
                'variable_sds' = variable_sds))
  } else {
    X <- subset(original_data, select = -c(EP))
    design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
    colnames(design_mat)[1] <- 'intercept'
    return(list('data' = original_data,
                'y' = original_data$EP,
                'X' = design_mat))
  }
}

power_plant <- load_pp_data('scripts/robust_regression/power_plant/power_plant.xlsx')

##### Sampling from full posterior #####

full_posterior <-  hmc_sample_BRR(noise_error = 'laplace',
                                  y = power_plant$y,
                                  X = power_plant$X,
                                  C = 1,
                                  sigma = sigma,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  iterations = nsamples + warmup,
                                  warmup = warmup,
                                  chains = 1,
                                  seed = seed,
                                  output = T)

##### Sampling from sub-posterior C=2 #####

data_split_2 <- split_data(power_plant$data,
                           y_col_index = 1,
                           X_col_index = 2:5,
                           C = 2,
                           as_dataframe = F)
sub_posteriors_2 <- hmc_base_sampler_BRR(noise_error = 'laplace',
                                         nsamples = nsamples,
                                         data_split = data_split_2,
                                         C = 2,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         warmup = warmup,
                                         seed = seed,
                                         output = T)

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(power_plant$data,
                           y_col_index = 1,
                           X_col_index = 2:5,
                           C = 4,
                           as_dataframe = F)
sub_posteriors_4 <- hmc_base_sampler_BRR(noise_error = 'laplace',
                                         nsamples = nsamples,
                                         data_split = data_split_4,
                                         C = 4,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         warmup = warmup,
                                         seed = seed,
                                         output = T)

##### Sampling from sub-posterior C=8 #####

data_split_8 <- split_data(power_plant$data,
                           y_col_index = 1,
                           X_col_index = 2:5,
                           C = 8,
                           as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_BRR(noise_error = 'laplace',
                                         nsamples = nsamples,
                                         data_split = data_split_8,
                                         C = 8,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         warmup = warmup,
                                         seed = seed,
                                         output = T)

##### Sampling from sub-posterior C=16 #####

data_split_16 <- split_data(power_plant$data,
                            y_col_index = 1,
                            X_col_index = 2:5,
                            C = 16,
                            as_dataframe = F)
sub_posteriors_16 <- hmc_base_sampler_BRR(noise_error = 'laplace',
                                          nsamples = nsamples,
                                          data_split = data_split_16,
                                          C = 16,
                                          sigma = sigma,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          warmup = warmup,
                                          seed = seed,
                                          output = T)
