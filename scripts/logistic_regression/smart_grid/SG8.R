library(DCFusion)
library(HMCBLR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples_MCF <- 10000
nsamples_GBF <- 10000
nsamples_DCGBF <- 10000
ndata <- 1000
time_choice <- 0.5
prior_means <- rep(0, 14)
prior_variances <- rep(1, 14)
C <- 8
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
CESS_j_threshold <- 0.2
diffusion_estimator <- 'NB'

##### Loading in Data #####

load_smart_grid_data <- function(file, seed = NULL) {
  original_data <- read.csv(file)
  original_data$stabf <- as.numeric(original_data$stabf == 'stable')
  original_data <- original_data[complete.cases(original_data),]
  if (!is.null(seed)) {
    set.seed(seed)
    original_data <- original_data[sample(1:nrow(original_data)),]
  }
  X <- subset(original_data, select = -c(stabf))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = original_data,
              'y' = original_data$stabf,
              'X' = design_mat))
}

smart_grid <- load_smart_grid_data('scripts/logistic_regression/smart_grid/smart_grid_stability_augmented.csv')

##### Sampling from full posterior #####

full_posterior <- hmc_sample_BLR_2(y = smart_grid$y,
                                   X = smart_grid$X,
                                   C = 1,
                                   prior_means = prior_means,
                                   prior_variances = prior_variances,
                                   iterations = nsamples_MCF + 10000,
                                   warmup = 10000,
                                   chains = 1,
                                   seed = seed,
                                   output = T)

##### Sampling from sub-posterior C=8 #####

data_split_8 <- split_data(smart_grid$data, y_col_index = 14, X_col_index = 1:13, C = C, as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_BLR_2(nsamples = nsamples_MCF,
                                           data_split = data_split_8,
                                           C = C, 
                                           prior_means = prior_means,
                                           prior_variances = prior_variances,
                                           warmup = 10000,
                                           seed = seed,
                                           output = T)

##### Applying other methodologies #####

# print('Applying other methodologies')
consensus_mat_8 <- consensus_scott(S = 8, samples_to_combine = sub_posteriors_8, indep = F)
consensus_sca_8 <- consensus_scott(S = 8, samples_to_combine = sub_posteriors_8, indep = T)
neiswanger_true_8 <- neiswanger(S = 8,
                                samples_to_combine = sub_posteriors_8,
                                anneal = TRUE)
neiswanger_false_8 <- neiswanger(S = 8,
                                 samples_to_combine = sub_posteriors_8,
                                 anneal = FALSE)
weierstrass_importance_8 <- weierstrass(Samples = sub_posteriors_8,
                                        method = 'importance')
weierstrass_rejection_8 <- weierstrass(Samples = sub_posteriors_8,
                                       method = 'reject')

integrated_abs_distance(full_posterior, consensus_mat_8$samples)
integrated_abs_distance(full_posterior, consensus_sca_8$samples)
integrated_abs_distance(full_posterior, neiswanger_true_8$samples)
integrated_abs_distance(full_posterior, neiswanger_false_8$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_8$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_8$samples)

# ##### NB (Hypercube Centre) #####
# print('NB Fusion (hypercube centre)')
# NB_hc_8 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples_MCF, 2),
#                                      m_schedule = rep(2, 2),
#                                      time_schedule = rep(time_choice, 2),
#                                      base_samples = sub_posteriors_8,
#                                      L = 3,
#                                      dim = 14,
#                                      data_split = data_split_8,
#                                      prior_means = prior_means,
#                                      prior_variances = prior_variances,
#                                      C = C,
#                                      precondition = TRUE,
#                                      resampling_method = 'resid',
#                                      ESS_threshold = ESS_threshold,
#                                      cv_location = 'hypercube_centre',
#                                      diffusion_estimator = 'NB',
#                                      seed = seed,
#                                      n_cores = n_cores)
# NB_hc_8$particles <- resample_particle_y_samples(particle_set = NB_hc_8$particles[[1]],
#                                                  multivariate = TRUE,
#                                                  resampling_method = 'resid',
#                                                  seed = seed)
# NB_hc_8$proposed_samples <- NB_hc_8$proposed_samples[[1]]
# print(integrated_abs_distance(full_posterior, NB_hc_8$particles$y_samples))

# ##### Generalised Bayesian Fusion #####
# 
# ##### all at once #####
# GBF_8 <- list('reg' = bal_binary_GBF_BLR(N_schedule = nsamples_GBF,
#                                          m_schedule = 8,
#                                          time_mesh = NULL,
#                                          base_samples = sub_posteriors_8,
#                                          L = 2,
#                                          dim = 14,
#                                          data_split = data_split_8,
#                                          prior_means = prior_means,
#                                          prior_variances = prior_variances,
#                                          C = C,
#                                          precondition = TRUE,
#                                          resampling_method = 'resid',
#                                          ESS_threshold = ESS_threshold,
#                                          adaptive_mesh = FALSE,
#                                          mesh_parameters = list('condition' = 'SH',
#                                                                 'CESS_0_threshold' = CESS_0_threshold,
#                                                                 'CESS_j_threshold' = CESS_j_threshold,
#                                                                 'vanilla' = FALSE),
#                                          diffusion_estimator = diffusion_estimator,
#                                          seed = seed))
# GBF_8$adaptive <-  bal_binary_GBF_BLR(N_schedule = nsamples_GBF,
#                                       m_schedule = 8,
#                                       time_mesh = NULL,
#                                       base_samples = sub_posteriors_8,
#                                       L = 2,
#                                       dim = 14,
#                                       data_split = data_split_8,
#                                       prior_means = prior_means,
#                                       prior_variances = prior_variances,
#                                       C = C,
#                                       precondition = TRUE,
#                                       resampling_method = 'resid',
#                                       ESS_threshold = ESS_threshold,
#                                       adaptive_mesh = TRUE,
#                                       mesh_parameters = list('condition' = 'SH',
#                                                              'CESS_0_threshold' = CESS_0_threshold,
#                                                              'CESS_j_threshold' = CESS_j_threshold,
#                                                              'vanilla' = FALSE),
#                                       diffusion_estimator = diffusion_estimator,
#                                       seed = seed)
# 
# # regular mesh
# GBF_8$reg$particles <- resample_particle_y_samples(particle_set = GBF_8$reg$particles[[1]],
#                                                    multivariate = TRUE,
#                                                    resampling_method = 'resid',
#                                                    seed = seed)
# print(integrated_abs_distance(full_posterior, GBF_8$reg$particles$y_samples))
# # adaptive mesh
# GBF_8$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_8$adaptive$particles[[1]],
#                                                         multivariate = TRUE,
#                                                         resampling_method = 'resid',
#                                                         seed = seed)
# print(integrated_abs_distance(full_posterior, GBF_8$adaptive$particles$y_samples))

##### bal binary combining two sub-posteriors at a time #####
balanced_C8 <- list('reg' = bal_binary_GBF_BLR(N_schedule = rep(nsamples_DCGBF, 2),
                                               m_schedule = rep(2, 2),
                                               time_mesh = NULL,
                                               base_samples = sub_posteriors_8,
                                               L = 3,
                                               dim = 14,
                                               data_split = data_split_8,
                                               prior_means = prior_means,
                                               prior_variances = prior_variances,
                                               C = C,
                                               precondition = TRUE,
                                               resampling_method = 'resid',
                                               ESS_threshold = ESS_threshold,
                                               adaptive_mesh = FALSE,
                                               mesh_parameters = list('condition' = 'SH',
                                                                      'CESS_0_threshold' = CESS_0_threshold,
                                                                      'CESS_j_threshold' = CESS_j_threshold,
                                                                      'vanilla' = FALSE),
                                               diffusion_estimator = diffusion_estimator,
                                               seed = seed,
                                               print_progress_iters = 100))
balanced_C8$adaptive <- bal_binary_GBF_BLR(N_schedule = rep(nsamples_DCGBF, 2),
                                           m_schedule = rep(2, 2),
                                           time_mesh = NULL,
                                           base_samples = sub_posteriors_8,
                                           L = 3,
                                           dim = 14,
                                           data_split = data_split_8,
                                           prior_means = prior_means,
                                           prior_variances = prior_variances,
                                           C = C,
                                           precondition = TRUE,
                                           resampling_method = 'resid',
                                           ESS_threshold = ESS_threshold,
                                           adaptive_mesh = TRUE,
                                           mesh_parameters = list('condition' = 'SH',
                                                                  'CESS_0_threshold' = CESS_0_threshold,
                                                                  'CESS_j_threshold' = CESS_j_threshold,
                                                                  'vanilla' = FALSE),
                                           diffusion_estimator = diffusion_estimator,
                                           seed = seed,
                                           print_progress_iters = 100)

# regular mesh
balanced_C8$reg$particles <- resample_particle_y_samples(particle_set = balanced_C8$reg$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
balanced_C8$reg$proposed_samples <- balanced_C8$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C8$reg$particles$y_samples))
# adaptive mesh
balanced_C8$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C8$adaptive$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
balanced_C8$adaptive$proposed_samples <- balanced_C8$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C8$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, GBF_8$reg$particles$y_samples)
integrated_abs_distance(full_posterior, GBF_8$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C8$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C8$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_8$particles$y_samples)

save.image('SG8.RData')
