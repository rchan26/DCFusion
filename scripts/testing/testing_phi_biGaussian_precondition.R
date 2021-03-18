library(hierarchicalFusion)

#####

seed <- 1983
set.seed(seed)

# setting parameters
mean <- c(1, 2)
sd <- c(sqrt(0.5), sqrt(2))
cov_mat <- matrix(c(sd[1]^2, 0.9, 0.9, sd[2]^2), nrow = 2, ncol = 2, byrow = T)
corr <- 0.9/(sd[1]*sd[2])
true_samples <- mvrnormArma(100000, mu = mean, Sigma = cov_mat)

ea_phi_biGaussian_DL_vec(x = true_samples[5,],
                         mean_vec = mean,
                         sd_vec = sd,
                         corr = corr,
                         beta = 1/8,
                         precondition_mat = cov(true_samples),
                         transform_mat = diag(1, 2))
terms_biGaussian_X(x = true_samples[5,],
                   mean_vec = mean,
                   sd_vec = sd,
                   corr = corr,
                   beta = 1/8,
                   precondition_mat = cov(true_samples),
                   sqrt_precondition_mat = expm::sqrtm(cov(true_samples)))
terms_biGaussian_Z(x = true_samples[1,],
                   mean_vec = mean,
                   sd_vec = sd,
                   corr = corr,
                   beta = 1/8,
                   sqrt_precondition_mat = expm::sqrtm(cov(true_samples)))

ea_phi_biGaussian_DL_vec(x = true_samples[1,],
                         mean_vec = mean,
                         sd_vec = sd,
                         corr = corr,
                         beta = 1/8,
                         precondition_mat = diag(1, 2),
                         transform_mat = diag(1, 2))
terms_biGaussian_X(x = true_samples[1,],
                   mean_vec = mean,
                   sd_vec = sd,
                   corr = corr,
                   beta = 1/8,
                   precondition_mat = diag(1, 2),
                   sqrt_precondition_mat = diag(1, 2))
terms_biGaussian_Z(x = true_samples[1,],
                   mean_vec = mean,
                   sd_vec = sd,
                   corr = corr,
                   beta = 1/8,
                   sqrt_precondition_mat = diag(1, 2))

