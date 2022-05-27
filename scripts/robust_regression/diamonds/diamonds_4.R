library(DCFusion)
library(HMCBRR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples_MCF <- 30000
nsamples_GBF <- 30000
nsamples_DCGBF <- 10000
warmup <- 10000
time_choice <- 1
nu <- 5
sigma <- 1
prior_means <- rep(0, 15)
prior_variances <- rep(10, 15)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.5
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

load_diamonds_data <- function(file, seed = NULL) {
  original_data <- read.csv(file)
  diamonds <- data.frame(price = original_data$price/1000,
                         carat = original_data$carat,
                         cut_good = as.numeric(original_data$cut == "Good"),
                         cut_ideal = as.numeric(original_data$cut == "Ideal"),
                         cut_prem = as.numeric(original_data$cut == "Premium"),
                         cut_vgood = as.numeric(original_data$cut == "Very Good"),
                         # badCol = as.numeric(original_data$color %in% c("J", "I", "H")),
                         medCol = as.numeric(original_data$color %in% c("G", "F")),
                         goodCol = as.numeric(original_data$color %in% c("E", "D")),
                         # colE = as.numeric(original_data$color == "E"),
                         # colF = as.numeric(original_data$color == "F"),
                         # colG = as.numeric(original_data$color == "G"),
                         # colH = as.numeric(original_data$color == "H"),
                         # colI = as.numeric(original_data$color == "I"),
                         # colJ = as.numeric(original_data$color == "J"),
                         # badClarity = as.numeric(original_data$clarity %in% c("I1", "SI1", "SI2")),
                         medClarity = as.numeric(original_data$clarity %in% c("VS2", "VS1")),
                         goodClarity = as.numeric(original_data$clarity %in% c("VVS2", "VVS1", "IF")),
                         # clarityIF = as.numeric(original_data$clarity == "IF"),
                         # claritySI1 = as.numeric(original_data$clarity == "SI1"),
                         # claritySI2 = as.numeric(original_data$clarity == "SI2"),
                         # clarityVS1 = as.numeric(original_data$clarity == "VS1"),
                         # clarityVS2 = as.numeric(original_data$clarity == "VS2"),
                         # clarityVVS1 = as.numeric(original_data$clarity == "VVS1"),
                         # clarityVVS2 = as.numeric(original_data$clarity == "VVS2"),
                         depth = original_data$depth,
                         table = original_data$table,
                         x_size = original_data$x,
                         y_size = original_data$y,
                         z_size = original_data$z)
  diamonds <- diamonds[complete.cases(diamonds),]
  if (!is.null(seed)) {
    set.seed(seed)
    diamonds <- diamonds[sample(1:nrow(diamonds)),]
  }
  means <- list('carat' = mean(diamonds$carat),
                'depth' = mean(diamonds$depth),
                'table' = mean(diamonds$table),
                'x_size' = mean(diamonds$x_size),
                'y_size' = mean(diamonds$y_size),
                'z_size' = mean(diamonds$z_size))
  sds <- list('carat' = sd(diamonds$carat),
              'depth' = sd(diamonds$depth),
              'table' = sd(diamonds$table),
              'x_size' = sd(diamonds$x_size),
              'y_size' = sd(diamonds$y_size),
              'z_size' = sd(diamonds$z_size))
  diamonds$carat <- (diamonds$carat-means$carat) / sds$carat
  diamonds$depth <- (diamonds$depth-means$depth) / sds$depth
  diamonds$table <- (diamonds$table-means$table) / sds$table
  diamonds$x_size <- (diamonds$x_size-means$x_size) / sds$x_size
  diamonds$y_size <- (diamonds$y_size-means$y_size) / sds$y_size
  diamonds$z_size <- (diamonds$z_size-means$z_size) / sds$z_size
  X <- subset(diamonds, select = -c(price))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind(subset(design_mat, select = -c(intercept)),
                             'price' = diamonds$price),
              'y' = diamonds$price,
              'X' = design_mat,
              'variable_means' = means,
              'variable_sds' = sds))
}

diamonds <- load_diamonds_data('scripts/robust_regression/diamonds/diamonds.csv')

# original_data <- read.csv('scripts/robust_regression/diamonds/diamonds.csv')
# original_data <- original_data[,-1]
# original_data$price <- original_data$price/1000
# lm(price ~., data = original_data)

##### Sampling from full posterior #####

full_posterior <-  hmc_sample_BRR(noise_error = 'student_t',
                                  y = diamonds$y,
                                  X = diamonds$X,
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

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(diamonds$data,
                           y_col_index = 15,
                           X_col_index = 1:14,
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

# compare_samples_bivariate(sub_posteriors_4, colours = c('red', 'blue', 'green', 'orange'), c(-4,4))

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

# ##### Poisson (Hypercube Centre) #####
# print('Poisson Fusion (hypercube centre)')
# Poisson_hc_4 <- bal_binary_fusion_SMC_BRR(N_schedule = rep(nsamples_MCF, 2),
#                                           m_schedule = rep(2, 2),
#                                           time_schedule = rep(time_choice, 2),
#                                           base_samples = sub_posteriors_4,
#                                           L = 3,
#                                           dim = 15,
#                                           data_split = data_split_4,
#                                           nu = nu,
#                                           sigma = sigma,
#                                           prior_means = prior_means,
#                                           prior_variances = prior_variances,
#                                           C = 4,
#                                           precondition = TRUE,
#                                           resampling_method = 'resid',
#                                           ESS_threshold = ESS_threshold,
#                                           cv_location = 'hypercube_centre',
#                                           diffusion_estimator = 'Poisson',
#                                           seed = seed,
#                                           n_cores = n_cores)
# Poisson_hc_4$particles <- resample_particle_y_samples(particle_set = Poisson_hc_4$particles[[1]],
#                                                       multivariate = TRUE,
#                                                       resampling_method = 'resid',
#                                                       seed = seed)
# Poisson_hc_4$proposed_samples <- Poisson_hc_4$proposed_samples[[1]]
# print(integrated_abs_distance(full_posterior, Poisson_hc_4$particles$y_samples))

##### NB (Hypercube Centre) #####
print('NB Fusion (hypercube centre)')
NB_hc_4 <- bal_binary_fusion_SMC_BRR(N_schedule = rep(nsamples_MCF, 2),
                                     m_schedule = rep(2, 2),
                                     time_schedule = rep(time_choice, 2),
                                     base_samples = sub_posteriors_4,
                                     L = 3,
                                     dim = 15,
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
                                         dim = 15,
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
                                         seed = seed))
GBF_4$adaptive <- bal_binary_GBF_BRR(N_schedule = nsamples_GBF,
                                     m_schedule = 4,
                                     time_mesh = NULL,
                                     base_samples = sub_posteriors_4,
                                     L = 2,
                                     dim = 15,
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
                                     seed = seed)

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
                                               dim = 15,
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
                                               seed = seed,
                                               print_progress_iters = 100))
balanced_C4$adaptive <- bal_binary_GBF_BRR(N_schedule = rep(nsamples_DCGBF, 2),
                                           m_schedule = rep(2, 2),
                                           time_mesh = NULL,
                                           base_samples = sub_posteriors_4,
                                           L = 3,
                                           dim = 15,
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
                                           seed = seed,
                                           print_progress_iters = 100)

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

save.image('diamonds4.RData')
