library(DCFusion)
library(HMCGLMR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples_MCF <- 10000
nsamples_GBF <- 10000
warmup <- 10000
time_choice <- 1
phi_rate <- 1
prior_means <- rep(0, 10)
prior_variances <- rep(10, 10)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.1
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

load_bs_data <- function(file, seed = NULL) {
  original_data <- read.csv(file)
  bike_sharing <- data.frame(rented = original_data$cnt,
                             spring = as.numeric(original_data$season == 2),
                             summer = as.numeric(original_data$season == 3),
                             fall = as.numeric(original_data$season == 4),
                             weekend = as.numeric(original_data$weekday %in% c(6, 0)),
                             holiday = original_data$holiday,
                             rush_hour = as.numeric(!(original_data$weekday %in% c(6, 0)) 
                                                    & (original_data$hr %in% c(7, 8, 9, 16, 17, 18, 19))),
                             weather = as.numeric(original_data$weathersit == 1),
                             atemp = original_data$atemp,
                             # humid = original_data$hum,
                             wind = original_data$windspeed)
  bike_sharing <- bike_sharing[complete.cases(bike_sharing),]
  if (!is.null(seed)) {
    set.seed(seed)
    bike_sharing <- bike_sharing[sample(1:nrow(bike_sharing)),]
  }
  X <- subset(bike_sharing, select = -c(rented))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind(subset(design_mat, select = -c(intercept)),
                             'rented' = bike_sharing$rented),
              'y' = bike_sharing$rented,
              'X' = design_mat))
}

bike_sharing <- load_bs_data('scripts/count_data_regression/bike_sharing/hour.csv', seed = seed)

# MASS::glm.nb(rented~., data = as.data.frame(bike_sharing$data))

##### Sampling from full posterior #####

full_posterior <- hmc_sample_GLMR(likelihood = 'NB',
                                  y = bike_sharing$y,
                                  X = bike_sharing$X,
                                  C = 1,
                                  phi = phi_rate,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  iterations = nsamples_MCF + 10000,
                                  warmup = 10000,
                                  chains = 1,
                                  seed = seed,
                                  output = T)

apply(full_posterior, 2, mean)

##### Sampling from sub-posterior C=256 #####

data_split_256 <- split_data(bike_sharing$data,
                             y_col_index = 10,
                             X_col_index = 1:9,
                             C = 256,
                             as_dataframe = F)
sub_posteriors_256 <- hmc_base_sampler_GLMR(likelihood = 'NB',
                                            nsamples = nsamples_MCF,
                                            warmup = 10000,
                                            data_split = data_split_256,
                                            C = 256,
                                            phi = phi_rate,
                                            prior_means = prior_means,
                                            prior_variances = prior_variances,
                                            seed = seed,
                                            output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_256 <- consensus_scott(S = 256, samples_to_combine = sub_posteriors_256, indep = F)
consensus_sca_256 <- consensus_scott(S = 256, samples_to_combine = sub_posteriors_256, indep = T)
neiswanger_true_256 <- neiswanger(S = 256,
                                  samples_to_combine = sub_posteriors_256,
                                  anneal = TRUE)
neiswanger_false_256 <- neiswanger(S = 256,
                                   samples_to_combine = sub_posteriors_256,
                                   anneal = FALSE)
weierstrass_importance_256 <- weierstrass(Samples = sub_posteriors_256,
                                          method = 'importance')
weierstrass_rejection_256 <- weierstrass(Samples = sub_posteriors_256,
                                         method = 'reject')

# ##### NB (Hypercube Centre) #####
# print('NB Fusion (hypercube centre)')
# NB_hc_256 <- bal_binary_fusion_SMC_BNBR(N_schedule = rep(nsamples_MCF, 8),
#                                        m_schedule = rep(2, 8),
#                                        time_schedule = rep(time_choice, 8),
#                                        base_samples = sub_posteriors_256,
#                                        L = 9,
#                                        dim = 10,
#                                        phi_rate = phi_rate,
#                                        data_split = data_split_256,
#                                        prior_means = prior_means,
#                                        prior_variances = prior_variances,
#                                        C = 256,
#                                        precondition = TRUE,
#                                        resampling_method = 'resid',
#                                        ESS_threshold = ESS_threshold,
#                                        diffusion_estimator = diffusion_estimator,
#                                        seed = seed,
#                                        n_cores = n_cores,
#                                        print_progress_iters = 5)
# NB_hc_256$particles <- resample_particle_y_samples(particle_set = NB_hc_256$particles[[1]],
#                                                   multivariate = TRUE,
#                                                   resampling_method = 'resid',
#                                                   seed = seed)
# NB_hc_256$proposed_samples <- NB_hc_256$proposed_samples[[1]]
# print(integrated_abs_distance(full_posterior, NB_hc_256$particles$y_samples))

# ##### Generalised Bayesian Fusion #####
# 
# ##### all at once #####
# GBF_256 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = nsamples_GBF,
#                                           m_schedule = 256,
#                                           time_mesh = NULL,
#                                           base_samples = sub_posteriors_256,
#                                           L = 2,
#                                           dim = 10,
#                                           phi_rate = phi_rate,
#                                           data_split = data_split_256,
#                                           prior_means = prior_means,
#                                           prior_variances = prior_variances,
#                                           C = 256,
#                                           precondition = TRUE,
#                                           resampling_method = 'resid',
#                                           ESS_threshold = ESS_threshold,
#                                           adaptive_mesh = FALSE,
#                                           mesh_parameters = list('condition' = 'SH',
#                                                                  'CESS_0_threshold' = CESS_0_threshold,
#                                                                  'CESS_j_threshold' = CESS_j_threshold,
#                                                                  'vanilla' = FALSE),
#                                           diffusion_estimator = diffusion_estimator,
#                                           seed = seed,
#                                           print_progress_iters = 5),
#               'adaptive' = bal_binary_GBF_BNBR(N_schedule = nsamples_GBF,
#                                                m_schedule = 256,
#                                                time_mesh = NULL,
#                                                base_samples = sub_posteriors_256,
#                                                L = 2,
#                                                dim = 10,
#                                                phi_rate = phi_rate,
#                                                data_split = data_split_256,
#                                                prior_means = prior_means,
#                                                prior_variances = prior_variances,
#                                                C = 256,
#                                                precondition = TRUE,
#                                                resampling_method = 'resid',
#                                                ESS_threshold = ESS_threshold,
#                                                adaptive_mesh = TRUE,
#                                                mesh_parameters = list('condition' = 'SH',
#                                                                       'CESS_0_threshold' = CESS_0_threshold,
#                                                                       'CESS_j_threshold' = CESS_j_threshold,
#                                                                       'vanilla' = FALSE),
#                                                diffusion_estimator = diffusion_estimator,
#                                                seed = seed,
#                                                print_progress_iters = 5))
# 
# # regular mesh
# GBF_256$reg$particles <- resample_particle_y_samples(particle_set = GBF_256$reg$particles[[1]],
#                                                    multivariate = TRUE,
#                                                    resampling_method = 'resid',
#                                                    seed = seed)
# print(integrated_abs_distance(full_posterior, GBF_256$reg$particles$y_samples)) 
# # adaptive mesh
# GBF_256$adaptive$particles <- resample_particle_y_samples(particle_set = GBF_256$adaptive$particles[[1]],
#                                                         multivariate = TRUE,
#                                                         resampling_method = 'resid',
#                                                         seed = seed)
# print(integrated_abs_distance(full_posterior, GBF_256$adaptive$particles$y_samples))

##### bal binary combining two sub-posteriors at a time #####
balanced_C256 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = rep(nsamples_GBF, 8),
                                                  m_schedule = rep(2, 8),
                                                  time_mesh = NULL,
                                                  base_samples = sub_posteriors_256,
                                                  L = 9,
                                                  dim = 10,
                                                  phi_rate = phi_rate,
                                                  data_split = data_split_256,
                                                  prior_means = prior_means,
                                                  prior_variances = prior_variances,
                                                  C = 256,
                                                  precondition = TRUE,
                                                  resampling_method = 'resid',
                                                  ESS_threshold = ESS_threshold,
                                                  adaptive_mesh = FALSE,
                                                  mesh_parameters = list('condition' = 'SH',
                                                                         'CESS_0_threshold' = CESS_0_threshold,
                                                                         'CESS_j_threshold' = CESS_j_threshold,
                                                                         'vanilla' = FALSE),
                                                  record = TRUE,
                                                  diffusion_estimator = diffusion_estimator,
                                                  seed = seed,
                                                  print_progress_iters = 50))
balanced_C256$adaptive <- bal_binary_GBF_BNBR(N_schedule = rep(nsamples_GBF, 8),
                                              m_schedule = rep(2, 8),
                                              time_mesh = NULL,
                                              base_samples = sub_posteriors_256,
                                              L = 9,
                                              dim = 10,
                                              phi_rate = phi_rate,
                                              data_split = data_split_256,
                                              prior_means = prior_means,
                                              prior_variances = prior_variances,
                                              C = 256,
                                              precondition = TRUE,
                                              resampling_method = 'resid',
                                              ESS_threshold = ESS_threshold,
                                              adaptive_mesh = TRUE,
                                              mesh_parameters = list('condition' = 'SH',
                                                                     'CESS_0_threshold' = CESS_0_threshold,
                                                                     'CESS_j_threshold' = CESS_j_threshold,
                                                                     'vanilla' = FALSE),
                                              record = TRUE,
                                              diffusion_estimator = diffusion_estimator,
                                              seed = seed,
                                              print_progress_iters = 50)

# regular mesh
balanced_C256$reg$particles <- resample_particle_y_samples(particle_set = balanced_C256$reg$particles[[1]],
                                                           multivariate = TRUE,
                                                           resampling_method = 'resid',
                                                           seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C256$reg$particles$y_samples))
# adaptive mesh
balanced_C256$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C256$adaptive$particles[[1]],
                                                                multivariate = TRUE,
                                                                resampling_method = 'resid',
                                                                seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C256$adaptive$particles$y_samples))

##### IAD #####

# integrated_abs_distance(full_posterior, GBF_256$reg$particles$y_samples)
# integrated_abs_distance(full_posterior, GBF_256$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C256$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C256$adaptive$particles$y_samples)
# integrated_abs_distance(full_posterior, Poisson_hc_256$particles$y_samples)
# integrated_abs_distance(full_posterior, NB_hc_256$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_256$samples)
integrated_abs_distance(full_posterior, consensus_sca_256$samples)
integrated_abs_distance(full_posterior, neiswanger_true_256$samples)
integrated_abs_distance(full_posterior, neiswanger_false_256$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_256$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_256$samples)

save.image('BS256_NB.RData')
