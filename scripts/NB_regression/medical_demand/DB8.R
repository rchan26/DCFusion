library(DCFusion)
library(HMCGLMR)

##### Initialise example #####
seed <- 2022
set.seed(seed)
nsamples <- 10000
warmup <- 10000
time_choice <- 1
phi_rate <- 1
prior_means <- rep(0, 10)
prior_variances <- rep(10, 10)
ESS_threshold <- 0.5
CESS_0_threshold <- 0.2
CESS_j_threshold <- 0.05
diffusion_estimator <- 'NB'
n_cores <- parallel::detectCores()

##### Loading in Data #####

load_medical_demand_data <- function(seed = NULL) {
  load('scripts/count_data_regression/medical_demand/DebTrivedi.rda')
  med_demand <- data.frame(visits = DebTrivedi$ofp,
                           hosp_stays = DebTrivedi$hosp,
                           exc_health = as.numeric(DebTrivedi$health == "excellent"),
                           avg_health = as.numeric(DebTrivedi$health == "average"),
                           n_chron = DebTrivedi$numchron,
                           age = DebTrivedi$age,
                           black = as.numeric(DebTrivedi$black == "yes"),
                           gender = as.numeric(DebTrivedi$gender == "male"),
                           employed = as.numeric(DebTrivedi$employed == "yes"),
                           priv_ins = as.numeric(DebTrivedi$privins == "yes"))
  med_demand <- med_demand[complete.cases(med_demand),]
  if (!is.null(seed)) {
    set.seed(seed)
    med_demand <- med_demand[sample(1:nrow(med_demand)),]
  }
  X <- subset(med_demand, select = -c(visits))
  design_mat <- as.matrix(cbind(rep(1, nrow(X)), X))
  colnames(design_mat)[1] <- 'intercept'
  return(list('data' = cbind(subset(design_mat, select = -c(intercept)),
                             'visits' = med_demand$visits),
              'y' = med_demand$visits,
              'X' = design_mat))
}

med_demand <- load_medical_demand_data()

##### Sampling from full posterior #####

full_posterior <- hmc_sample_GLMR(likelihood = 'NB',
                                  y = med_demand$y,
                                  X = med_demand$X,
                                  C = 1,
                                  phi = 1,
                                  prior_means = prior_means,
                                  prior_variances = prior_variances,
                                  iterations = nsamples + 10000,
                                  warmup = 10000,
                                  chains = 1,
                                  seed = seed,
                                  output = T)

##### Sampling from sub-posterior C=8 #####

data_split_8 <- split_data(med_demand$data,
                           y_col_index = 10,
                           X_col_index = 1:9,
                           C = 8,
                           as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_GLMR(likelihood = 'NB',
                                          nsamples = nsamples,
                                          warmup = 10000,
                                          data_split = data_split_8,
                                          C = 8,
                                          phi = 1,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          seed = seed,
                                          output = T)

##### Applying other methodologies #####

print('Applying other methodologies')
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

# ##### Generalised Bayesian Fusion #####

##### bal binary combining two sub-posteriors at a time #####
balanced_C8 <- list('reg' = bal_binary_GBF_BNBR(N_schedule = rep(nsamples, 3),
                                                m_schedule = rep(2, 3),
                                                time_mesh = NULL,
                                                base_samples = sub_posteriors_8,
                                                L = 4,
                                                dim = 10,
                                                phi_rate = phi_rate,
                                                data_split = data_split_8,
                                                prior_means = prior_means,
                                                prior_variances = prior_variances,
                                                C = 8,
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
balanced_C8$adaptive <- bal_binary_GBF_BNBR(N_schedule = rep(nsamples, 3),
                                            m_schedule = rep(2, 3),
                                            time_mesh = NULL,
                                            base_samples = sub_posteriors_8,
                                            L = 4,
                                            dim = 10,
                                            phi_rate = phi_rate,
                                            data_split = data_split_8,
                                            prior_means = prior_means,
                                            prior_variances = prior_variances,
                                            C = 8,
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

integrated_abs_distance(full_posterior, balanced_C8$reg$particles$y_samples)
integrated_abs_distance(full_posterior, balanced_C8$adaptive$particles$y_samples)
integrated_abs_distance(full_posterior, consensus_mat_8$samples)
integrated_abs_distance(full_posterior, consensus_sca_8$samples)
integrated_abs_distance(full_posterior, neiswanger_true_8$samples)
integrated_abs_distance(full_posterior, neiswanger_false_8$samples)
integrated_abs_distance(full_posterior, weierstrass_importance_8$samples)
integrated_abs_distance(full_posterior, weierstrass_rejection_8$samples)

save.image('DB8.RData')
