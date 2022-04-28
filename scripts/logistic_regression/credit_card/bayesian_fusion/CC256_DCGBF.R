library(DCFusion)
library(HMCBLR)

seed <- 2016
nsamples <- 30000
nsamples_GBF <- 5000
time_choice <- 0.5
C <- 256
n_cores <- parallel::detectCores()
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'

##### Loading in Data #####

original_data <- read.csv('credit_cards.csv')
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

##### Sampling from full posterior #####

credit_card_data <- list()
credit_card_data$y <- credit_cards_full$y
credit_card_data$X <- as.matrix(cbind(rep(1, nrow(credit_cards_full[,2:5])), credit_cards_full[,2:5]))
colnames(credit_card_data$X)[1] <- 'intercept'
full_data_count <- unique_row_count(credit_card_data$y, credit_card_data$X)$full_data_count

full_posterior <- hmc_sample_BLR(full_data_count = full_data_count,
                                 C = 1,
                                 prior_means = rep(0, 5),
                                 prior_variances = rep(1, 5),
                                 iterations = nsamples + 10000,
                                 warmup = 10000,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=64 #####

data_split_256 <- split_data(credit_cards_full, y_col_index = 1, X_col_index = 2:5, C = 256, as_dataframe = F)
sub_posteriors_256 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                           data_split = data_split_256,
                                           C = 256, 
                                           prior_means = rep(0, 5),
                                           prior_variances = rep(1, 5),
                                           warmup = 10000,
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

integrated_abs_distance(full_posterior, consensus_mat_256$samples)
integrated_abs_distance(full_posterior, consensus_sca_256$samples)
integrated_abs_distance(full_posterior, neiswanger_true_256$samples)
integrated_abs_distance(full_posterior, neiswanger_false_256$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_256$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_256$samples)

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_256 <- bal_binary_fusion_SMC_BLR(N_schedule = rep(nsamples, 8),
                                       m_schedule = rep(2, 8),
                                       time_schedule = rep(time_choice, 8),
                                       base_samples = sub_posteriors_256,
                                       L = 9,
                                       dim = 5,
                                       data_split = data_split_256,
                                       prior_means = rep(0, 5),
                                       prior_variances = rep(1, 5),
                                       C = 256,
                                       precondition = TRUE,
                                       resampling_method = 'resid',
                                       ESS_threshold = 0.5,
                                       cv_location = 'hypercube_centre',
                                       diffusion_estimator = 'NB',
                                       seed = seed,
                                       n_cores = n_cores,
                                       print_progress_iters = 50)
NB_hc_256$particles <- resample_particle_y_samples(particle_set = NB_hc_256$particles[[1]],
                                                   multivariate = TRUE,
                                                   resampling_method = 'resid',
                                                   seed = seed)
NB_hc_256$proposed_samples <- NB_hc_256$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, NB_hc_256$particles$y_samples))

##### bal binary combining two sub-posteriors at a time #####
balanced_C256 <- list('reg' = bal_binary_GBF_BLR(N_schedule = rep(nsamples_GBF, 8),
                                                 m_schedule = rep(2, 8),
                                                 time_mesh = NULL,
                                                 base_samples = sub_posteriors_256,
                                                 L = 9,
                                                 dim = 5,
                                                 data_split = data_split_256,
                                                 prior_means = rep(0, 5),
                                                 prior_variances = rep(1, 5),
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
                                                 print_progress_iters = 50,
                                                 seed = seed))
balanced_C256$adaptive <-  bal_binary_GBF_BLR(N_schedule = rep(nsamples_GBF, 8),
                                              m_schedule = rep(2, 8),
                                              time_mesh = NULL,
                                              base_samples = sub_posteriors_256,
                                              L = 9,
                                              dim = 5,
                                              data_split = data_split_256,
                                              prior_means = rep(0, 5),
                                              prior_variances = rep(1, 5),
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
                                              print_progress_iters = 50,
                                              seed = seed)

# regular mesh
balanced_C256$reg$particles <- resample_particle_y_samples(particle_set = balanced_C256$reg$particles[[1]],
                                                           multivariate = TRUE,
                                                           resampling_method = 'resid',
                                                           seed = seed)
balanced_C256$reg$proposed_samples <- balanced_C256$reg$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C256$reg$particles$y_samples))
# adaptive mesh
balanced_C256$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C256$adaptive$particles[[1]],
                                                                multivariate = TRUE,
                                                                resampling_method = 'resid',
                                                                seed = seed)
balanced_C256$adaptive$proposed_samples <- balanced_C256$adaptive$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior, balanced_C256$adaptive$particles$y_samples))

save.image('CC256_DCGBF.RData')
