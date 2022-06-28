library(DCFusion)
library(HMCGLMR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples <- 50000
warmup <- 10000
time_choice <- 1
phi_rate <- 1
dim <- 35
prior_means <- rep(0, dim)
prior_variances <- rep(10, dim)
C <- 4
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
  X <- matrix(data = 0, nrow = 2*nrow(full_data), ncol = 2*length(teams)+1)
  colnames(X) <- c(paste(teams, "att", sep = "_"), paste(teams, "def", sep = "_"), "Home")
  y <- rep(NA, 2*nrow(full_data))
  for (i in 1:nrow(full_data)) {
    # home goals
    X[2*(i-1)+1, c(paste(full_data[i,]$HomeTeam, "att", sep = "_"), "Home")] <- 1 
    X[2*(i-1)+1, paste(full_data[i,]$AwayTeam, "def", sep = "_")] <- -1 
    y[2*(i-1)+1] <- full_data[i,]$FTHG
    # away goals
    X[2*(i-1)+2, paste(full_data[i,]$AwayTeam, "att", sep = "_")] <- 1
    X[2*(i-1)+2, paste(full_data[i,]$HomeTeam, "def", sep = "_")] <- -1
    y[2*(i-1)+2] <- full_data[i,]$FTAG
  }
  # column for first team is removed and taken as baseline
  # X <- X[,-1]
  # design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  # colnames(design_mat)[1] <- 'intercept'
  return(list('original_data' = full_data,
              'data' = cbind(X, 'goals' = y),
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
                                  iterations = nsamples + 10000,
                                  warmup = 10000,
                                  chains = 1,
                                  seed = seed,
                                  output = T)
apply(full_posterior, 2, mean)

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(scottish_football$data,
                           y_col_index = dim+1,
                           X_col_index = 1:dim,
                           C = C,
                           as_dataframe = F)
# remove intercept column
for (i in 1:4) {
  data_split_4[[i]]$X <- data_split_4[[i]]$X[,-1]
  data_split_4[[i]]$full_data_count <- subset(data_split_4[[i]]$full_data_count, select = -V2)
  data_split_4[[i]]$design_count <- subset(data_split_4[[i]]$design_count, select = -V1)
}
sub_posteriors_4 <- hmc_base_sampler_GLMR(likelihood = 'NB_2',
                                          nsamples = nsamples,
                                          warmup = 10000,
                                          data_split = data_split_4,
                                          C = C,
                                          phi = phi_rate,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          seed = seed,
                                          output = T)
for (i in 1:4) {
  print(apply(sub_posteriors_4[[i]], 2, mean))
}

##### Applying other methodologies #####

print('Applying other methodologies')
consensus_mat_4 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_4, indep = F)
consensus_sca_4 <- consensus_scott(S = C, samples_to_combine = sub_posteriors_4, indep = T)
neiswanger_true_4 <- neiswanger(S = C,
                                samples_to_combine = sub_posteriors_4,
                                anneal = TRUE)
neiswanger_false_4 <- neiswanger(S = C,
                                 samples_to_combine = sub_posteriors_4,
                                 anneal = FALSE)
weierstrass_importance_4 <- weierstrass(Samples = sub_posteriors_4,
                                        method = 'importance')
weierstrass_rejection_4 <- weierstrass(Samples = sub_posteriors_4,
                                       method = 'reject')

integrated_abs_distance(full_posterior, consensus_mat_4$samples)
integrated_abs_distance(full_posterior, consensus_sca_4$samples)
integrated_abs_distance(full_posterior, neiswanger_true_4$samples)
integrated_abs_distance(full_posterior, neiswanger_false_4$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_4$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_4$samples)

##### bal binary combining two sub-posteriors at a time #####
balanced_C4 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = rep(nsamples, 2),
                                                m_schedule = rep(2, 2),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_4,
                                                L = 3,
                                                dim = dim,
                                                phi_rate = phi_rate,
                                                data_split = data_split_4,
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
balanced_C4$adaptive <- bal_binary_GBF_BNBR(N_schedule = rep(nsamples, 2),
                                            m_schedule = rep(2, 2),
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_4,
                                            L = 3,
                                            dim = dim,
                                            phi_rate = phi_rate,
                                            data_split = data_split_4,
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

save.image('SF4_50000.RData')
