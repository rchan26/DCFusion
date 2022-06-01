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
prior_means <- rep(0, 18)
prior_variances <- rep(10, 18)
C <- 16
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

load_football_data <- function(list_of_files, seed = NULL) {
  original_data <- lapply(list_of_files, function(file) {
    subset(read.csv(file,fileEncoding="latin1"), select=c('HomeTeam', "AwayTeam", "FTHG", "FTAG"))})
  full_data <- do.call(rbind, original_data)
  full_data <- full_data[complete.cases(full_data),]
  if (!is.null(seed)) {
    set.seed(seed)
    full_data <- full_data[sample(1:nrow(full_data)),]
  }
  teams <- as.character(sort(unique(full_data$HomeTeam)))
  X <- matrix(data = 0, nrow = 2*nrow(full_data), ncol = length(teams)+1)
  colnames(X) <- c(teams, "Home")
  y <- rep(NA, 2*nrow(full_data))
  for (i in 1:nrow(full_data)) {
    # home goals
    X[2*(i-1)+1, c(full_data[i,]$HomeTeam, full_data[i,]$AwayTeam, "Home")] <- 1 
    y[2*(i-1)+1] <- full_data[i,]$FTHG
    # away goals
    X[2*(i-1)+2, c(full_data[i,]$HomeTeam, full_data[i,]$AwayTeam)] <- 1
    y[2*(i-1)+2] <- full_data[i,]$FTAG
  }
  # column for first team is removed and taken as baseline
  # X <- X[,-1]
  # design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  # colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind(X, 'goals' = y),
              'y' = y,
              'X' = X))
}

list_of_files <- list('scripts/count_data_regression/football_data/SC0_1112.csv',
                      'scripts/count_data_regression/football_data/SC0_1213.csv',
                      'scripts/count_data_regression/football_data/SC0_1314.csv',
                      'scripts/count_data_regression/football_data/SC0_1415.csv',
                      'scripts/count_data_regression/football_data/SC0_1516.csv',
                      'scripts/count_data_regression/football_data/SC0_1617.csv',
                      'scripts/count_data_regression/football_data/SC0_1718.csv',
                      'scripts/count_data_regression/football_data/SC0_1819.csv',
                      'scripts/count_data_regression/football_data/SC0_1920.csv',
                      'scripts/count_data_regression/football_data/SC0_2021.csv')
scottish_football <- load_football_data(list_of_files)

glm(scottish_football$y~0+scottish_football$X, family = poisson)

##### Sampling from full posterior #####

full_posterior <- hmc_sample_GLMR(likelihood = 'NB_2',
                                  y = scottish_football$y,
                                  X = scottish_football$X,
                                  C = 1,
                                  phi = phi_rate,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  iterations = nsamples_MCF + 10000,
                                  warmup = 10000,
                                  chains = 1,
                                  seed = seed,
                                  output = T)

##### Sampling from sub-posterior C=16 #####

data_split_16 <- split_data(scottish_football$data,
                           y_col_index = 19,
                           X_col_index = 1:18,
                           C = C,
                           as_dataframe = F)
# remove intercept column
for (i in 1:16) {
  data_split_16[[i]]$X <- data_split_16[[i]]$X[,-1]
  data_split_16[[i]]$full_data_count <- subset(data_split_16[[i]]$full_data_count, select = -V2)
  data_split_16[[i]]$design_count <- subset(data_split_16[[i]]$design_count, select = -V1)
}
sub_posteriors_16 <- hmc_base_sampler_GLMR(likelihood = 'NB_2',
                                          nsamples = nsamples_MCF,
                                          warmup = 10000,
                                          data_split = data_split_16,
                                          C = C,
                                          phi = phi_rate,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          seed = seed,
                                          output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_16 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_16, indep = F)
consensus_sca_16 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_16, indep = T)
neiswanger_true_16 <- neiswanger(S = C,
                                samples_to_combine = sub_posteriors_16,
                                anneal = TRUE)
neiswanger_false_16 <- neiswanger(S = C,
                                 samples_to_combine = sub_posteriors_16,
                                 anneal = FALSE)
weierstrass_importance_16 <- weierstrass(Samples = sub_posteriors_16,
                                        method = 'importance')
weierstrass_rejection_16 <- weierstrass(Samples = sub_posteriors_16,
                                       method = 'reject')

##### bal binary combining two sub-posteriors at a time #####
balanced_C16 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = rep(nsamples_GBF, 4),
                                                m_schedule = rep(2, 4),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_16,
                                                L = 5,
                                                dim = 18,
                                                phi_rate = phi_rate,
                                                data_split = data_split_16,
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
                                                record = TRUE,
                                                diffusion_estimator = diffusion_estimator,
                                                seed = seed,
                                                print_progress_iters = 100))
balanced_C16$adaptive <- bal_binary_GBF_BNBR(N_schedule = rep(nsamples_GBF, 4),
                                            m_schedule = rep(2, 4),
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_16,
                                            L = 5,
                                            dim = 18,
                                            phi_rate = phi_rate,
                                            data_split = data_split_16,
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
                                            record = TRUE,
                                            diffusion_estimator = diffusion_estimator,
                                            seed = seed,
                                            print_progress_iters = 100)

# regular mesh
balanced_C16$reg$particles <- resample_particle_y_samples(particle_set = balanced_C16$reg$particles[[1]],
                                                         multivariate = TRUE,
                                                         resampling_method = 'resid',
                                                         seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C16$reg$particles$y_samples))
# adaptive mesh
balanced_C16$adaptive$particles <- resample_particle_y_samples(particle_set = balanced_C16$adaptive$particles[[1]],
                                                              multivariate = TRUE,
                                                              resampling_method = 'resid',
                                                              seed = seed)
print(integrated_abs_distance(full_posterior, balanced_C16$adaptive$particles$y_samples))

##### IAD #####

integrated_abs_distance(full_posterior, consensus_mat_16$samples)
integrated_abs_distance(full_posterior, consensus_sca_16$samples)
integrated_abs_distance(full_posterior, neiswanger_true_16$samples)
integrated_abs_distance(full_posterior, neiswanger_false_16$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_16$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_16$samples)

save.image('SF16.RData')
