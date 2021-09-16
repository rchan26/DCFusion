library(hierarchicalFusion)
library(HMCBLR)

seed <- 2016
nsamples <- 30000
time_choice <- 0.5
n_cores <- parallel::detectCores()

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

data_split_64 <- split_data(credit_cards_full, y_col_index = 1, X_col_index = 2:5, C = 64, as_dataframe = F)
sub_posteriors_64 <- hmc_base_sampler_BLR(nsamples = nsamples,
                                          data_split = data_split_64,
                                          C = 64, 
                                          prior_means = rep(0, 5),
                                          prior_variances = rep(1, 5),
                                          warmup = 10000,
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

##### Applying Fusion #####

##### Poisson (Hypercube Centre) #####
print('Poisson Fusion (hypercube centre)')
Poisson_hc_64 <- hierarchical_fusion_SMC_BLR(N_schedule = rep(nsamples, 6),
                                             m_schedule = rep(2, 6),
                                             time_schedule = rep(time_choice, 6),
                                             base_samples = sub_posteriors_64,
                                             L = 7,
                                             dim = 5,
                                             data_split = data_split_64,
                                             prior_means = rep(0, 5),
                                             prior_variances = rep(1, 5),
                                             C = 64,
                                             precondition = TRUE,
                                             resampling_method = 'resid',
                                             ESS_threshold = 0.5,
                                             cv_location = 'hypercube_centre',
                                             diffusion_estimator = 'Poisson',
                                             seed = seed,
                                             n_cores = n_cores,
                                             print_progress_iters = 500)
Poisson_hc_64$particles <- resample_particle_y_samples(particle_set = Poisson_hc_64$particles[[1]],
                                                       multivariate = TRUE,
                                                       resampling_method = 'resid',
                                                       seed = seed)
Poisson_hc_64$proposed_samples <- Poisson_hc_64$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior,
                              Poisson_hc_64$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_64 <- hierarchical_fusion_SMC_BLR(N_schedule = rep(nsamples, 6),
                                        m_schedule = rep(2, 6),
                                        time_schedule = rep(time_choice, 6),
                                        base_samples = sub_posteriors_64,
                                        L = 7,
                                        dim = 5,
                                        data_split = data_split_64,
                                        prior_means = rep(0, 5),
                                        prior_variances = rep(1, 5),
                                        C = 64,
                                        precondition = TRUE,
                                        resampling_method = 'resid',
                                        ESS_threshold = 0.5,
                                        cv_location = 'hypercube_centre',
                                        diffusion_estimator = 'NB',
                                        seed = seed,
                                        n_cores = n_cores,
                                        print_progress_iters = 500)
NB_hc_64$particles <- resample_particle_y_samples(particle_set = NB_hc_64$particles[[1]],
                                                  multivariate = TRUE,
                                                  resampling_method = 'resid',
                                                  seed = seed)
NB_hc_64$proposed_samples <- NB_hc_64$proposed_samples[[1]]
print(integrated_abs_distance(full_posterior,
                              NB_hc_64$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, Poisson_hc_64$particles$y_samples)
integrated_abs_distance(full_posterior, NB_hc_64$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_64$samples)
integrated_abs_distance(full_posterior, consensus_sca_64$samples)
integrated_abs_distance(full_posterior, neiswanger_true_64$samples)
integrated_abs_distance(full_posterior, neiswanger_false_64$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_64$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_64$samples)

##### Save data #####

save.image('CC64.RData')
