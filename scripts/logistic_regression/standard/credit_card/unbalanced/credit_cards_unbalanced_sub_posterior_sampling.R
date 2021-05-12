library(hierarchicalFusion)
library(HMCBLR)

seed <- 2016

##### Loading in Data #####

original_data <- read.csv('../credit_cards.csv', header = T)
original_data <- original_data[2:nrow(original_data),]
credit_cards <- original_data[,c(25, 3, 4)]
colnames(credit_cards) <- c('y', 'sex', 'education')

# y (default payment): 1 for yes, 0 for no
credit_cards_full <- data.frame(y = as.numeric(credit_cards$y==1))
# sex: 1 for male, 0 for female
credit_cards_full$sex <- as.numeric(credit_cards$sex==1)
# education: ED1: graduate university
credit_cards_full$ED1 <- as.numeric(credit_cards$education==1)
# education: ED2: undergraduate university
credit_cards_full$ED2 <- as.numeric(credit_cards$education==2)
# education: ED3: highschool 
credit_cards_full$ED3 <- as.numeric(credit_cards$education==3)

# function to check the frequencies of 'active' variables in the dataset
check_activity <- function(data, proportion = T) {
  if (!is.data.frame(data)) {
    stop("check_activity: data is not in a data frame format")
  }
  
  # create new data frame with same column names
  active_df <- data.frame(matrix(nrow = 1, ncol = ncol(data)))
  colnames(active_df) <- colnames(data)
  
  # loop through columns and check how many are 'active'
  # i.e. have a 1 in the instance
  for (j in 1:ncol(data)) {
    if (proportion) {
      active_df[,j] <- sum(data[,j]==1) / nrow(data)
    } else {
      active_df[,j] <- sum(data[,j]==1)
    }
  }
  
  print('proportion that each variable is active in the dataset:')
  for (j in 1:ncol(active_df)) {
    print(paste(colnames(active_df)[j], ':', active_df[1,j]))
  }
}

check_activity(credit_cards_full)

##### Making the data unbalanced #####

length(which(credit_cards_full$y==1))

credit_cards_full <- credit_cards_full[order(credit_cards_full$y), ]

n_swap <- 5000
no_default_sample_indices <- sample(1:23364, n_swap)
default_sample_indices <- sample(23365:30000, n_swap)
for (i in 1:n_swap) {
  transfer_to_default <- credit_cards_full[no_default_sample_indices[i],]
  transfer_to_no_default <- credit_cards_full[default_sample_indices[i],]
  # swap the values
  credit_cards_full[no_default_sample_indices[i],] <- transfer_to_no_default
  credit_cards_full[default_sample_indices[i],] <- transfer_to_default
}

y <- credit_cards_full$y
X <- as.matrix(cbind(rep(1, nrow(credit_cards_full[,2:5])), credit_cards_full[,2:5]))
colnames(X)[1] <- 'intercept'

##### Sampling from full posterior #####

credit_card_data <- list()
credit_card_data$y <- credit_cards_full$y
credit_card_data$X <- as.matrix(cbind(rep(1, nrow(credit_cards_full[,2:5])), credit_cards_full[,2:5]))
colnames(credit_card_data$X)[1] <- 'intercept'

full_posterior <- hmc_sample_BLR(data = credit_card_data, 
                                 C = 1, 
                                 prior_means = rep(0, 5),
                                 prior_variances = rep(1, 5), 
                                 iterations = 40000, 
                                 warmup = 10000, 
                                 chains = 1, 
                                 power = 1,
                                 seed = seed,
                                 output = T)

standard_normal_prior <- function(beta, C) {
  dim <- length(beta)
  return(sum(dnorm(x = beta, mean = 0, sd = sqrt(C), log = TRUE)))
}

full_posterior_mcmc <- as.matrix(MCMCpack::MCMClogit(formula = y ~ sex + ED1 + ED2 + ED3, 
                                                     data = credit_cards_full,
                                                     burnin = 10000, 
                                                     mcmc = 40000, 
                                                     user.prior.density = standard_normal_prior,
                                                     C = 1,
                                                     logfun = TRUE))

glm(formula = y ~ sex + ED1 + ED2 + ED3, 
    data = credit_cards_full, 
    family = 'binomial')

apply(full_posterior, 2, mean)
apply(full_posterior_mcmc, 2, mean)

##### Sampling from sub-posterior C=2 #####

data_split_2 <- split_data(credit_cards_full, y_col_index = 1, X_col_index = 2:5, C = 2, as_dataframe = F)
sub_posteriors_2 <- hmc_base_sampler_BLR(nsamples = 30000,
                                         data_split = data_split_2,
                                         C = 2, 
                                         prior_means = rep(0, 5),
                                         prior_variances = rep(1, 5),
                                         warmup = 10000,
                                         seed = seed,
                                         output = T)

compare_samples_bivariate(sub_posteriors_2, c('red' ,'blue'), c(-4, 4))

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(credit_cards_full, y_col_index = 1, X_col_index = 2:5, C = 4, as_dataframe = F)
sub_posteriors_4 <- hmc_base_sampler_BLR(nsamples = 30000,
                                         data_split = data_split_4,
                                         C = 4, 
                                         prior_means = rep(0, 5),
                                         prior_variances = rep(1, 5),
                                         warmup = 10000,
                                         seed = seed,
                                         output = T)

compare_samples_bivariate(sub_posteriors_4, c(rep('red', 2), rep('blue', 2)), c(-4, 4))

##### Sampling from sub-posterior C=8 #####

data_split_8 <- split_data(credit_cards_full, y_col_index = 1, X_col_index = 2:5, C = 8, as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_BLR(nsamples = 30000,
                                         data_split = data_split_8,
                                         C = 8, 
                                         prior_means = rep(0, 5),
                                         prior_variances = rep(1, 5),
                                         warmup = 10000,
                                         seed = seed,
                                         output = T)

compare_samples_bivariate(sub_posteriors_8, c(rep('red', 2), rep('blue', 2), rep('green', 2), rep('purple', 2)), c(-4, 4))

bandwidths <- rep(0.1, 5)

consensus_mat_8 <- consensus_scott(S = 8, samples_to_combine = sub_posteriors_8, indep = F)
consensus_sca_8 <- consensus_scott(S = 8, samples_to_combine = sub_posteriors_8, indep = T)
integrated_abs_distance(full_posterior, consensus_mat_8$samples, bandwidths)
integrated_abs_distance(full_posterior, consensus_sca_8$samples, bandwidths)

neiswanger_8_true <- neiswanger(S = 8, 
                                samples_to_combine = sub_posteriors_8, 
                                bw = bandwidths, 
                                anneal = TRUE)
neiswanger_8_false <- neiswanger(S = 8, 
                                 samples_to_combine = sub_posteriors_8, 
                                 bw = bandwidths, 
                                 anneal = FALSE)
integrated_abs_distance(full_posterior, neiswanger_8_true$samples, bw = bandwidths)
integrated_abs_distance(full_posterior, neiswanger_8_false$samples, bw = bandwidths)

weierstrass_8_importance <- weierstrass(Samples = sub_posteriors_8, method = 'importance')
weierstrass_8_rejection <- weierstrass(Samples = sub_posteriors_8, method = 'reject')
integrated_abs_distance(full_posterior, weierstrass_8_importance$samples, bandwidths)
integrated_abs_distance(full_posterior, weierstrass_8_rejection$samples, bandwidths)

##### Sampling from sub-posterior C=16 #####

data_split_16 <- split_data(dataframe = credit_cards_full, y_col_index = 1, X_col_index = 2:5, C = 16, as_dataframe = F)
sub_posteriors_16 <- hmc_base_sampler_BLR(nsamples = 30000,
                                          data_split = data_split_16,
                                          C = 16, 
                                          prior_means = rep(0, 5),
                                          prior_variances = rep(1, 5),
                                          warmup = 10000,
                                          seed = seed,
                                          output = T)

consensus_mat_16 <- consensus_scott(S = 16, samples_to_combine = sub_posteriors_16, indep = F)
consensus_sca_16 <- consensus_scott(S = 16, samples_to_combine = sub_posteriors_16, indep = T)
integrated_abs_distance(full_posterior, consensus_mat_16$samples, bandwidths)
integrated_abs_distance(full_posterior, consensus_sca_16$samples, bandwidths)

neiswanger_16_true <- neiswanger(S = 16, 
                                 samples_to_combine = sub_posteriors_16, 
                                 bw = bandwidths, 
                                 anneal = TRUE)
neiswanger_16_false <- neiswanger(S = 16, 
                                  samples_to_combine = sub_posteriors_16, 
                                  bw = bandwidths, 
                                  anneal = FALSE)
integrated_abs_distance(full_posterior, neiswanger_16_true$samples, bw = bandwidths)
integrated_abs_distance(full_posterior, neiswanger_16_false$samples, bw = bandwidths)

weierstrass_16_importance <- weierstrass(Samples = sub_posteriors_16, method = 'importance')
weierstrass_16_rejection <- weierstrass(Samples = sub_posteriors_16, method = 'reject')
integrated_abs_distance(full_posterior, weierstrass_16_importance$samples, bandwidths)
integrated_abs_distance(full_posterior, weierstrass_16_rejection$samples, bandwidths)
##### Sampling from sub-posterior C=32 #####

data_split_32 <- split_data(credit_cards_full, y_col_index = 1, X_col_index = 2:5, C = 32, as_dataframe = F)
sub_posteriors_32 <- hmc_base_sampler_BLR(nsamples = 30000,
                                          data_split = data_split_32,
                                          C = 32, 
                                          prior_means = rep(0, 5),
                                          prior_variances = rep(1, 5),
                                          warmup = 10000,
                                          seed = seed,
                                          output = T)

consensus_mat_32 <- consensus_scott(S = 32, samples_to_combine = sub_posteriors_32, indep = F)
consensus_sca_32 <- consensus_scott(S = 32, samples_to_combine = sub_posteriors_32, indep = T)
integrated_abs_distance(full_posterior, consensus_mat_32$samples, bw = bandwidths)
integrated_abs_distance(full_posterior, consensus_sca_32$samples, bw = bandwidths)

neiswanger_32_true <- neiswanger(S = 32, 
                                 samples_to_combine = sub_posteriors_32, 
                                 bw = bandwidths, 
                                 anneal = TRUE)
neiswanger_32_false <- neiswanger(S = 32, 
                                  samples_to_combine = sub_posteriors_32, 
                                  bw = bandwidths, 
                                  anneal = FALSE)
integrated_abs_distance(full_posterior, neiswanger_32_true$samples, bw = bandwidths)
integrated_abs_distance(full_posterior, neiswanger_32_false$samples, bw = bandwidths)

weierstrass_32_importance <- weierstrass(Samples = sub_posteriors_32,
                                         method = 'importance')
weierstrass_32_rejection <- weierstrass(Samples = sub_posteriors_32,
                                        method = 'reject')
integrated_abs_distance(full_posterior, weierstrass_32_importance$samples, bw = bandwidths)
integrated_abs_distance(full_posterior, weierstrass_32_rejection$samples, bw = bandwidths)
