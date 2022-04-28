library(DCFusion)
library(HMCBRR)
library(readxl)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples_MCF <- 100000
nsamples_DCGBF <- 10000
warmup <- 10000
time_choice <- 1
nu <- 5
sigma <- 1
prior_means <- rep(0, 5)
prior_variances <- rep(10, 5)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

# Features consist of hourly average ambient variables
# - Temperature (T) in the range 1.81°C and 37.11°C,
# - Ambient Pressure (AP) in the range 992.89-1033.30 millibar,
# - Relative Humidity (RH) in the range 25.56% to 100.16%
# - Exhaust Vacuum (V) in the range 25.36-81.56 cm Hg
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

full_posterior <-  hmc_sample_BRR(noise_error = 'student_t',
                                  y = power_plant$y,
                                  X = power_plant$X,
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

##### Sampling from sub-posterior C=64 #####

data_split_64 <- split_data(power_plant$data,
                            y_col_index = 1,
                            X_col_index = 2:5,
                            C = 64,
                            as_dataframe = F)
sub_posteriors_64 <- hmc_base_sampler_BRR(noise_error = 'student_t',
                                          nsamples = nsamples_MCF,
                                          data_split = data_split_64,
                                          C = 64,
                                          nu = nu,
                                          sigma = sigma,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          warmup = warmup,
                                          seed = seed,
                                          output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_64 <- consensus_scott(S = 64, samples_to_combine = sub_posteriors_64, indep = F)
consensus_sca_64 <- consensus_scott(S = 64, samples_to_combine = sub_posteriors_64, indep = T)
neiswanger_true_64 <- neiswanger(S = 64,
                                 samples_to_combine = sub_posteriors_64,
                                 anneal = TRUE)
neiswanger_false_64 <- neiswanger(S = 64,
                                  samples_to_combine = sub_posteriors_64,
                                  anneal = FALSE)
weierstrass_importance_64 <- weierstrass(Samples = sub_posteriors_64,
                                         method = 'importance')
weierstrass_rejection_64 <- weierstrass(Samples = sub_posteriors_64,
                                        method = 'reject')

# ##### Poisson (Hypercube Centre) #####
# print('Poisson Fusion (hypercube centre)')
# Poisson_hc_64 <- bal_binary_fusion_SMC_BRR(N_schedule = rep(nsamples_MCF, 6),
#                                            m_schedule = rep(2, 6),
#                                            time_schedule = rep(time_choice, 6),
#                                            base_samples = sub_posteriors_64,
#                                            L = 7,
#                                            dim = 5,
#                                            data_split = data_split_64,
#                                            nu = nu,
#                                            sigma = sigma,
#                                            prior_means = prior_means,
#                                            prior_variances = prior_variances,
#                                            C = 64,
#                                            precondition = TRUE,
#                                            resampling_method = 'resid',
#                                            ESS_threshold = ESS_threshold,
#                                            cv_location = 'hypercube_centre',
#                                            diffusion_estimator = 'Poisson',
#                                            seed = seed,
#                                            n_cores = n_cores)
# Poisson_hc_64$particles <- resample_particle_y_samples(particle_set = Poisson_hc_64$particles[[1]],
#                                                        multivariate = TRUE,
#                                                        resampling_method = 'resid',
#                                                        seed = seed)
# Poisson_hc_64$proposed_samples <- Poisson_hc_64$proposed_samples[[1]]
# print(integrated_abs_distance(full_posterior, Poisson_hc_64$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_64 <- bal_binary_fusion_SMC_BRR(N_schedule = rep(nsamples_MCF, 6),
                                      m_schedule = rep(2, 6),
                                      time_schedule = rep(time_choice, 6),
                                      base_samples = sub_posteriors_64,
                                      L = 7,
                                      dim = 5,
                                      data_split = data_split_64,
                                      nu = nu,
                                      sigma = sigma,
                                      prior_means = prior_means,
                                      prior_variances = prior_variances,
                                      C = 64,
                                      precondition = TRUE,
                                      resampling_method = 'resid',
                                      ESS_threshold = ESS_threshold,
                                      cv_location = 'hypercube_centre',
                                      diffusion_estimator = 'NB',
                                      seed = seed,
                                      n_cores = n_cores)
NB_hc_64$particles <- resample_particle_y_samples(particle_set = NB_hc_64$particles[[1]],
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed)
NB_hc_64$proposed_samples <- NB_hc_64$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_64$particles$y_samples))

##### Generalised Bayesian Fusion #####

##### bal binary combining two sub-posteriors at a time #####
balanced_C64 <- list('reg' = bal_binary_GBF_BRR(N_schedule = rep(nsamples_DCGBF, 6),
                                                m_schedule = rep(2, 6),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_64,
                                                L = 7,
                                                dim = 5,
                                                data_split = data_split_64,
                                                nu = nu,
                                                sigma = sigma,
                                                prior_means = prior_means,
                                                prior_variances = prior_variances,
                                                C = 64,
                                                precondition = TRUE,
                                                resampling_method = 'resid',
                                                ESS_threshold = ESS_threshold,
                                                adaptive_mesh = FALSE,
                                                mesh_parameters = list('condition' = 'SH',
                                                                       'CESS_0_threshold' = CESS_0_threshold,
                                                                       'CESS_j_threshold' = CESS_j_threshold,
                                                                       'vanilla' = FALSE),
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed))
balanced_C64$adaptive <- bal_binary_GBF_BRR(N_schedule = rep(nsamples_DCGBF, 6),
                                            m_schedule = rep(2, 6),
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_64,
                                            L = 7,
                                            dim = 5,
                                            data_split = data_split_64,
                                            nu = nu,
                                            sigma = sigma,
                                            prior_means = prior_means,
                                            prior_variances = prior_variances,
                                            C = 64,
                                            precondition = TRUE,
                                            resampling_method = 'resid',
                                            ESS_threshold = ESS_threshold,
                                            adaptive_mesh = TRUE,
                                            mesh_parameters = list('condition' = 'SH',
                                                                   'CESS_0_threshold' = CESS_0_threshold,
                                                                   'CESS_j_threshold' = CESS_j_threshold,
                                                                   'vanilla' = FALSE),
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed)

# regular mesh
balanced_C64$reg$particles <- resample_particle_y_samples(particle_set = balanced_C64$reg$particles[[1]],
                                                          multivariate = TRUE,
                                                          resampling_method = 'resid',
                                                          seed = seed)
balanced_C64$reg$proposed_samples <- balanced_C64$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C64$reg$particles$y_samples))
# adaptive mesh
balanced_C64$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C64$adaptive$particles[[1]],
                                                               multivariate = TRUE,
                                                               resampling_method = 'resid',
                                                               seed = seed)
balanced_C64$adaptive$proposed_samples <- balanced_C64$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C64$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, balanced_C64$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C64$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_64$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_64$samples)
integrated_abs_distance(full_posterior, consensus_sca_64$samples)
integrated_abs_distance(full_posterior, neiswanger_true_64$samples)
integrated_abs_distance(full_posterior, neiswanger_false_64$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_64$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_64$samples)

save.image('PP64.RData')
