library(hierarchicalFusion)

seed <- 1994
set.seed(seed)

input_samples <- lapply(1:8, function(i) rnorm_tempered(N = 50000,
                                                        mean = 1.2, 
                                                        sd = sqrt(0.1), 
                                                        beta = 1/8))

###### Hierarchical [Poisson] #####

hier_standard <- hierarchical_fusion_uniGaussian(N_schedule = rep(50000, 3), 
                                                         m_schedule = rep(2, 3),
                                                         time_schedule = rep(1, 3), 
                                                         base_samples = input_samples, 
                                                         L = 4, 
                                                         mean = 1.2, 
                                                         sd = sqrt(0.1), 
                                                         start_beta = 1/8,
                                                         precondition = FALSE,
                                                         seed = seed)

hier_precondition <- hierarchical_fusion_uniGaussian(N_schedule = rep(50000, 3), 
                                                             m_schedule = rep(2, 3),
                                                             time_schedule = rep(1, 3), 
                                                             base_samples = input_samples, 
                                                             L = 4, 
                                                             mean = 1.2, 
                                                             sd = sqrt(0.1), 
                                                             start_beta = 1/8,
                                                             precondition = TRUE, 
                                                             seed = seed)

# acceptance rate comparisons
acceptance_rate_plots(hier1 = hier_precondition,
                      hier2 = hier_standard,
                      time = 1,
                      hierarchical = TRUE)

curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -1, 3)
lines(density(hier_standard$samples[[1]]), col = 'green')
lines(density(hier_precondition$samples[[1]]), col = 'blue')

ks.test(hier_standard$samples[[1]], rnorm(n = 50000, mean = 1.2, sd = sqrt(0.1)))
ks.test(hier_precondition$samples[[1]], rnorm(n = 50000, mean = 1.2, sd = sqrt(0.1)))

###### Progressive [Poisson] #####

prog_standard <- progressive_fusion_uniGaussian(N_schedule = rep(50000, 7),
                                                        time_schedule = rep(1, 7), 
                                                        base_samples = input_samples,
                                                        mean = 1.2, 
                                                        sd = sqrt(0.1), 
                                                        start_beta = 1/8,
                                                        precondition = FALSE,
                                                        seed = seed)

prog_precondition <- progressive_fusion_uniGaussian(N_schedule = rep(50000, 7),
                                                            time_schedule = rep(1, 7), 
                                                            base_samples = input_samples,
                                                            mean = 1.2, 
                                                            sd = sqrt(0.1), 
                                                            start_beta = 1/8,
                                                            precondition = TRUE,
                                                            seed = seed)

# acceptance rate comparisons
acceptance_rate_plots(hier1 = prog_precondition,
                      hier2 = prog_standard,
                      time = 1,
                      hierarchical = FALSE)

curve(dnorm(x, mean = 1.2, sd = sqrt(0.1)), -1, 3)
lines(density(prog_standard$samples[[1]]), col = 'green')
lines(density(prog_precondition$samples[[1]]), col = 'blue')

ks.test(prog_standard$samples[[1]], rnorm(n = 50000, mean = 1.2, sd = sqrt(0.1)))
ks.test(prog_precondition$samples[[1]], rnorm(n = 50000, mean = 1.2, sd = sqrt(0.1)))
