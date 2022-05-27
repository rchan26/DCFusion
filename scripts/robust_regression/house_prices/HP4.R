library(DCFusion)
library(HMCBRR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples_MCF <- 10000
nsamples_GBF <- 10000
warmup <- 10000
time_choice <- 1
nu <- 5
sigma <- 1
prior_means <- rep(0, 12)
prior_variances <- rep(10, 12)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.2
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

load_house_prices_data <- function(file, seed = NULL) {
  original_data <- read.csv(file)
  house_prices <- data.frame(price = original_data$price/1000,
                             bedrooms = original_data$bedrooms,
                             bathrooms = original_data$bathrooms,
                             sqft_living = original_data$sqft_living,
                             sqft_lot = original_data$sqft_lot,
                             sqft_above = original_data$sqft_above,
                             sqft_basement = original_data$sqft_basement,
                             floors = original_data$floors,
                             condition = as.numeric(original_data$condition > 3),
                             view = as.numeric(original_data$view > 0),
                             avg_grade = as.numeric(7 <= original_data$grade & original_data$grade <= 9),
                             good_grade = as.numeric(original_data$grade > 9))
  house_prices <- house_prices[complete.cases(house_prices),]
  if (!is.null(seed)) {
    set.seed(seed)
    house_prices <- house_prices[sample(1:nrow(house_prices)),]
  }
  means <- list('sqft_living' = mean(house_prices$sqft_living),
                'sqft_lot' = mean(house_prices$sqft_lot),
                'sqft_above' = mean(house_prices$sqft_above),
                'sqft_basement' = mean(house_prices$sqft_basement))
  sds <- list('sqft_living' = sd(house_prices$sqft_living),
              'sqft_lot' = sd(house_prices$sqft_lot),
              'sqft_above' = sd(house_prices$sqft_above),
              'sqft_basement' = sd(house_prices$sqft_basement))
  house_prices$sqft_living <- (house_prices$sqft_living-means$sqft_living) / sds$sqft_living
  house_prices$sqft_lot <- (house_prices$sqft_lot-means$sqft_lot) / sds$sqft_lot
  house_prices$sqft_above <- (house_prices$sqft_above-means$sqft_above) / sds$sqft_above
  house_prices$sqft_basement <- (house_prices$sqft_basement-means$sqft_basement) / sds$sqft_basement
  X <- subset(house_prices, select = -c(price))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind(subset(design_mat, select = -c(intercept)),
                             'price' = house_prices$price),
              'y' = house_prices$price,
              'X' = design_mat,
              'variable_means' = means,
              'variable_sds' = sds))
}

house_prices <- load_house_prices_data('scripts/robust_regression/house_prices/kc_house_data.csv', seed = seed)

##### Sampling from full posterior #####

full_posterior <-  hmc_sample_BRR(noise_error = 'student_t',
                                  y = house_prices$y,
                                  X = house_prices$X,
                                  C = 1,
                                  nu = nu,
                                  sigma = sigma,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  iterations = nsamples_MCF + warmup,
                                  warmup = warmup,
                                  chains = 1,
                                  seed = seed,
                                  output = T)

apply(full_posterior, 2, mean)

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(house_prices$data,
                           y_col_index = 12,
                           X_col_index = 1:11,
                           C = 4,
                           as_dataframe = F)
sub_posteriors_4 <- hmc_base_sampler_BRR(noise_error = 'student_t',
                                         nsamples = nsamples_MCF,
                                         data_split = data_split_4,
                                         C = 4,
                                         nu = nu,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         warmup = warmup,
                                         seed = seed,
                                         output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_4 <- consensus_scott(S = 4, samples_to_combine = sub_posteriors_4, indep = F)
consensus_sca_4 <- consensus_scott(S = 4, samples_to_combine = sub_posteriors_4, indep = T)
neiswanger_true_4 <- neiswanger(S = 4,
                                samples_to_combine = sub_posteriors_4,
                                anneal = TRUE)
neiswanger_false_4 <- neiswanger(S = 4,
                                 samples_to_combine = sub_posteriors_4,
                                 anneal = FALSE)
weierstrass_importance_4 <- weierstrass(Samples = sub_posteriors_4,
                                        method = 'importance')
weierstrass_rejection_4 <- weierstrass(Samples = sub_posteriors_4,
                                       method = 'reject')

for (i in 1:4) {
  print(apply(sub_posteriors_4[[i]], 2, mean))
}
##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_4 <- bal_binary_fusion_SMC_BRR(N_schedule = rep(nsamples_MCF, 2),
                                     m_schedule = rep(2, 2),
                                     time_schedule = rep(time_choice, 2),
                                     base_samples = sub_posteriors_4,
                                     L = 3,
                                     dim = 5,
                                     data_split = data_split_4,
                                     nu = nu,
                                     sigma = sigma,
                                     prior_means = prior_means,
                                     prior_variances = prior_variances,
                                     C = 4,
                                     precondition = TRUE,
                                     resampling_method = 'resid',
                                     ESS_threshold = ESS_threshold,
                                     cv_location = 'hypercube_centre',
                                     diffusion_estimator = 'NB',
                                     seed = seed,
                                     n_cores = n_cores)
NB_hc_4$particles <- resample_particle_y_samples(particle_set = NB_hc_4$particles[[1]],
                                                 multivariate = TRUE,
                                                 resampling_method = 'resid',
                                                 seed = seed)
NB_hc_4$proposed_samples <- NB_hc_4$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_4$particles$y_samples))

##### Generalised Bayesian Fusion #####

##### all at once #####
GBF_4 <- list('reg' = bal_binary_GBF_BRR(N_schedule = nsamples_GBF,
                                         m_schedule = 4,
                                         time_mesh = NULL,
                                         base_samples = sub_posteriors_4,
                                         L = 2,
                                         dim = 5,
                                         data_split = data_split_4,
                                         nu = nu,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         C = 4,
                                         precondition = TRUE,
                                         resampling_method = 'resid',
                                         ESS_threshold = ESS_threshold,
                                         adaptive_mesh = FALSE,
                                         mesh_parameters = list('condition' = 'SH',
                                                                'CESS_0_threshold' = CESS_0_threshold,
                                                                'CESS_j_threshold' = CESS_j_threshold,
                                                                'vanilla' = FALSE),
                                         diffusion_estimator = diffusion_estimator,
                                         seed = seed),
              'adaptive' = bal_binary_GBF_BRR(N_schedule = nsamples_GBF,
                                              m_schedule = 4,
                                              time_mesh = NULL,
                                              base_samples = sub_posteriors_4,
                                              L = 2,
                                              dim = 5,
                                              data_split = data_split_4,
                                              nu = nu,
                                              sigma = sigma,
                                              prior_means = prior_means,
                                              prior_variances = prior_variances,
                                              C = 4,
                                              precondition = TRUE,
                                              resampling_method = 'resid',
                                              ESS_threshold = ESS_threshold,
                                              adaptive_mesh = TRUE,
                                              mesh_parameters = list('condition' = 'SH',
                                                                     'CESS_0_threshold' = CESS_0_threshold,
                                                                     'CESS_j_threshold' = CESS_j_threshold,
                                                                     'vanilla' = FALSE),
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed))

# regular mesh
GBF_4$reg$particles <- resample_particle_y_samples(particle_set = GBF_4$reg$particles[[1]],
                                                   multivariate = TRUE,
                                                   resampling_method = 'resid',
                                                   seed = seed)
print(integrated_abs_distance(full_posterior, GBF_4$reg$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            GBF_4$reg$proposed_samples[[1]],
                                            GBF_4$reg$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))
# adaptive mesh
GBF_4$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_4$adaptive$particles[[1]],
                                                        multivariate = TRUE,
                                                        resampling_method = 'resid',
                                                        seed = seed)
print(integrated_abs_distance(full_posterior, GBF_4$adaptive$particles$y_samples))
compare_samples_bivariate(posteriors = list(full_posterior,
                                            GBF_4$adaptive$proposed_samples[[1]],
                                            GBF_4$adaptive$particles$y_samples),
                          colours = c('black', 'green', 'red'),
                          common_limit = c(-4, 4))

##### bal binary combining two sub-posteriors at a time #####
balanced_C4 <- list('reg' = bal_binary_GBF_BRR(N_schedule = rep(nsamples_DCGBF, 2),
                                               m_schedule = rep(2, 2),
                                               time_mesh = NULL,
                                               base_samples = sub_posteriors_4,
                                               L = 3,
                                               dim = 5,
                                               data_split = data_split_4,
                                               nu = nu,
                                               sigma = sigma,
                                               prior_means = prior_means,
                                               prior_variances = prior_variances,
                                               C = 4,
                                               precondition = TRUE,
                                               resampling_method = 'resid',
                                               ESS_threshold = ESS_threshold,
                                               adaptive_mesh = FALSE,
                                               mesh_parameters = list('condition' = 'SH',
                                                                      'CESS_0_threshold' = CESS_0_threshold,
                                                                      'CESS_j_threshold' = CESS_j_threshold,
                                                                      'vanilla' = FALSE),
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed),
                    'adaptive' = bal_binary_GBF_BRR(N_schedule = rep(nsamples_DCGBF, 2),
                                                    m_schedule = rep(2, 2),
                                                    time_mesh = NULL,
                                                    base_samples = sub_posteriors_4,
                                                    L = 3,
                                                    dim = 5,
                                                    data_split = data_split_4,
                                                    nu = nu,
                                                    sigma = sigma,
                                                    prior_means = prior_means,
                                                    prior_variances = prior_variances,
                                                    C = 4,
                                                    precondition = TRUE,
                                                    resampling_method = 'resid',
                                                    ESS_threshold = ESS_threshold,
                                                    adaptive_mesh = TRUE,
                                                    mesh_parameters = list('condition' = 'SH',
                                                                           'CESS_0_threshold' = CESS_0_threshold,
                                                                           'CESS_j_threshold' = CESS_j_threshold,
                                                                           'vanilla' = FALSE),
                                                    diffusion_estimator = diffusion_estimator,
                                                    seed = seed))

# regular mesh
balanced_C4$reg$particles <- resample_particle_y_samples(particle_set = balanced_C4$reg$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
balanced_C4$reg$proposed_samples <- balanced_C4$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C4$reg$particles$y_samples))
# adaptive mesh
balanced_C4$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C4$adaptive$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
balanced_C4$adaptive$proposed_samples <- balanced_C4$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C4$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, GBF_4$reg$particles$y_samples)
integrated_abs_distance(full_posterior, GBF_4$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C4$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C4$adaptive$particles$y_samples)
# integrated_abs_distance(full_posterior, Poisson_hc_4$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_4$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_4$samples)
integrated_abs_distance(full_posterior, consensus_sca_4$samples)
integrated_abs_distance(full_posterior, neiswanger_true_4$samples)
integrated_abs_distance(full_posterior, neiswanger_false_4$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_4$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_4$samples)

save.image('HP4.RData')
