library(hierarchicalFusion)

######################################## examples ########################################

seed <- 21
set.seed(seed)

# rejection sampling target
# dominating_dnorm <- function(x) 1.35*dnorm(x, mean = 0, sd = 1)
# curve(dominating_dnorm(x), -5, 5, col = 'red')
# curve(dnorm(x, mean = 0, sd = 1), -5, 5, add = T, col = 'green')
# curve(exp_4_density(x), -5, 5, add = T, col = 'blue')
test_target_mc <- sample_exp_4(N = 100000,
                               proposal_mean = 0,
                               proposal_sd = 1,
                               dominating_M = 1.35,
                               beta = 1)

######################################## beta = 1/4

# sampling from the tempered targets
# finding the right proposal distribution and constant M in rejection sampling
dominating_dnorm <- function(x) 1.4*dnorm(x, mean = 0, sd = 1.5)
curve(dominating_dnorm(x), -5, 5, col = 'red')
curve(dnorm(x, mean = 0, sd = 1.5), -5, 5, add = T, col = 'green')
curve(exp_4_density(x, 0, 1/4), -5, 5, add = T, col = 'blue')

# using rejection sampling to obtain input samples
seed <- 21
set.seed(seed)
input_samples1 <- base_rejection_sampler_exp_4(beta = 1/4,
                                               nsamples = 100000,
                                               proposal_mean = 0,
                                               proposal_sd = 1.5,
                                               dominating_M = 1.4)
curve(exp_4_density(x, beta = 1/4), -4, 4)
# check the samples look okay
for (samples in input_samples1) {
  lines(density(samples), col = 'blue')
}

# standard hierarchical fusion
test1_hier_standard <- hierarchical_fusion_exp_4(N_schedule = rep(100000, 2),
                                                 m_schedule = rep(2, 2),
                                                 time_schedule = rep(1, 2),
                                                 base_samples = input_samples1, 
                                                 mean = 0,
                                                 start_beta = 1/4,
                                                 L = 3,
                                                 precondition = FALSE,
                                                 seed = seed)

# preconditioned hierarchical fusion
test1_hier_precondition <- hierarchical_fusion_exp_4(N_schedule = rep(100000, 2),
                                                     m_schedule = rep(2, 2),
                                                     time_schedule = rep(1, 2),
                                                     base_samples = input_samples1, 
                                                     mean = 0,
                                                     start_beta = 1/4,
                                                     L = 3,
                                                     precondition = TRUE,
                                                     seed = seed)

curve(exp_4_density(x), -4, 4, ylim = c(0,0.5), lwd = 2)
lines(density(test_target_mc), col = 'black')
lines(density(test1_hier_standard$samples[[1]]), col = 'green')
lines(density(test1_hier_precondition$samples[[1]]), col = 'blue')

# acceptance rate comparisons
acceptance_rate_plots(hier1 = test1_hier_precondition,
                      hier2 = test1_hier_standard,
                      time = 1,
                      hierarchical = TRUE)

######################################## beta = 1/8

# sampling from the tempered targets
# finding the right proposal distribution and constant M in rejection sampling
# dominating_dnorm <- function(x) 1.3*dnorm(x, mean = 0, sd = 1.5)
# curve(dominating_dnorm(x), -5, 5, col = 'red')
# curve(dnorm(x, mean = 0, sd = 1.5), -5, 5, add = T, col = 'green')
# curve(exp_4_density(x, 0, 1/8), -5, 5, add = T, col = 'blue')

# using rejection sampling to obtain input samples
seed <- 21
set.seed(seed)
input_samples2 <- base_rejection_sampler_exp_4(beta = 1/8,
                                               nsamples = 100000,
                                               proposal_mean = 0,
                                               proposal_sd = 1.5,
                                               dominating_M = 1.3)
# curve(exp_4_density(x, beta = 1/8), -4, 4)
# # check the samples look okay
# for (samples in input_samples2) {
#   lines(density(samples), col = 'blue')
# }

# standard hierarchical fusion
test2_hier_standard <- hierarchical_fusion_exp_4(N_schedule = rep(100000, 3),
                                                 m_schedule = rep(2, 3),
                                                 time_schedule = rep(1, 3),
                                                 base_samples = input_samples2,
                                                 mean = 0,
                                                 start_beta = 1/8,
                                                 L = 4,
                                                 precondition = FALSE,
                                                 seed = seed)

# preconditioned hierarchical fusion
test2_hier_precondition <- hierarchical_fusion_exp_4(N_schedule = rep(100000, 3),
                                                     m_schedule = rep(2, 3),
                                                     time_schedule = rep(1, 3),
                                                     base_samples = input_samples2,
                                                     mean = 0,
                                                     start_beta = 1/8,
                                                     L = 4,
                                                     precondition = TRUE,
                                                     seed = seed)

curve(exp_4_density(x), -4, 4, ylim = c(0,0.5))
lines(density(test_target_mc), col = 'black')
lines(density(test2_hier_standard$samples[[1]]), col = 'green')
lines(density(test2_hier_precondition$samples[[1]]), col = 'blue')

# acceptance rate comparisons
acceptance_rate_plots(hier1 = test2_hier_precondition,
                      hier2 = test2_hier_standard,
                      time = 1,
                      hierarchical = TRUE)

# #####
#
# par(mfrow=c(1,2))
#
# #####
# # hier
#
# #rho
# plot(1:length(test2_hier_standard$overall_rho), test2_hier_standard$overall_rho, ylim = c(0,1), col = 'black',
#      ylab = expression(rho), xlab = 'Level', pch = 4, cex = 0.5)
# lines(1:length(test2_hier_standard$overall_rho), test2_hier_standard$overall_rho, col = 'black')
# points(1:length(test2_hier_precondition$overall_rho), test2_hier_precondition$overall_rho, col = 'orange', pch = 4, cex = 0.5)
# lines(1:length(test2_hier_precondition$overall_rho), test2_hier_precondition$overall_rho, col = 'orange', lty = 2)
#
# # Q
# plot(1:length(test2_hier_standard$overall_Q), test2_hier_standard$overall_Q, ylim = c(0,1), col = 'black',
#      ylab = 'Q', xlab = 'Level', pch = 4, cex = 0.5)
# lines(1:length(test2_hier_standard$overall_Q), test2_hier_standard$overall_Q, col = 'black')
# points(1:length(test2_hier_precondition$overall_Q), test2_hier_precondition$overall_Q, col = 'orange', pch = 4, cex = 0.5)
# lines(1:length(test2_hier_precondition$overall_Q), test2_hier_precondition$overall_Q, col = 'orange', lty = 2)
#
# #####
#
# par(mfrow=c(1,1))
#
# # rho
# round(test2_hier_standard$overall_rho, 3)
# round(test2_hier_precondition$overall_rho, 3)
#
# # Q
# round(test2_hier_standard$overall_Q, 3)
# round(test2_hier_precondition$overall_Q, 3)

######################################## beta = 1/16

# sampling from the tempered targets
# finding the right proposal distribution and constant M in rejection sampling
# dominating_dnorm <- function(x) 1.4*dnorm(x, mean = 0, sd = 1.5)
# curve(dominating_dnorm(x), -5, 5, col = 'red')
# curve(dnorm(x, mean = 0, sd = 1.5), -5, 5, add = T, col = 'green')
# curve(exp_4_density(x, 0, 1/16), -5, 5, add = T, col = 'blue')

# using rejection sampling to obtain input samples
seed <- 21
set.seed(seed)
input_samples3 <- base_rejection_sampler_exp_4(beta = 1/16,
                                               nsamples = 100000,
                                               proposal_mean = 0,
                                               proposal_sd = 1.5,
                                               dominating_M = 1.4)
curve(exp_4_density(x, beta = 1/16), -4, 4)
# check the samples look okay
for (samples in input_samples3) {
  lines(density(samples), col = 'blue')
}
sapply(input_samples3, var)

# standard hierarchical fusion
test3_hier_standard <- hierarchical_fusion_exp_4(N_schedule = rep(100000, 4),
                                                 m_schedule = rep(2, 4),
                                                 time_schedule = rep(1, 4),
                                                 base_samples = input_samples3,
                                                 mean = 0,
                                                 start_beta = 1/16,
                                                 L = 5,
                                                 precondition = FALSE,
                                                 seed = seed)

# preconditioned hierarchical fusion
test3_hier_precondition <- hierarchical_fusion_exp_4(N_schedule = rep(100000, 4),
                                                     m_schedule = rep(2, 4),
                                                     time_schedule = rep(1, 4),
                                                     base_samples = input_samples3,
                                                     mean = 0,
                                                     start_beta = 1/16,
                                                     L = 5,
                                                     precondition = TRUE,
                                                     seed = seed)

curve(exp_4_density(x), -4, 4)
lines(density(test_target_mc), col = 'black')
lines(density(test3_hier_standard$samples[[1]]), col = 'green')
lines(density(test3_hier_precondition$samples[[1]]), col = 'blue')

# acceptance rate comparisons
acceptance_rate_plots(hier1 = test3_hier_precondition,
                      hier2 = test3_hier_standard,
                      time = 1,
                      hierarchical = TRUE)

#####

par(oma = c(2, 1, 0, 1))
par(mfrow = c(1,2))

#####
# hier

#rho
plot(1:length(test3_hier_standard$overall_rho), test3_hier_standard$overall_rho, ylim = c(0,1), col = 'black',
     ylab = expression(rho), xlab = 'Level', pch = 4, lwd = 1.5, xaxt = 'n')
axis(1, at = 1:4)
lines(1:length(test3_hier_standard$overall_rho), test3_hier_standard$overall_rho, col = 'black', lwd = 1.5)
points(1:length(test3_hier_precondition$overall_rho), test3_hier_precondition$overall_rho, col = 'magenta', pch = 2, lwd = 1.5)
lines(1:length(test3_hier_precondition$overall_rho), test3_hier_precondition$overall_rho, col = 'magenta', lty = 2, lwd = 1.5)

# Q
plot(1:length(test3_hier_standard$overall_Q), test3_hier_standard$overall_Q, ylim = c(0,1), col = 'black',
     ylab = 'Q', xlab = 'Level', pch = 4, lwd = 1.5, xaxt = 'n')
axis(1, at = 1:4)
lines(1:length(test3_hier_standard$overall_Q), test3_hier_standard$overall_Q, col = 'black', lwd = 1.5)
points(1:length(test3_hier_precondition$overall_Q), test3_hier_precondition$overall_Q, col = 'magenta', pch = 2, lwd = 1.5)
lines(1:length(test3_hier_precondition$overall_Q), test3_hier_precondition$overall_Q, col = 'magenta', lty = 2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",
       legend = c('standard', 'preconditioned'),
       lty = c(1, 2),
       xpd = TRUE,
       horiz = TRUE,
       inset = c(0, 0),
       bty = "n",
       pch = c(4, 2),
       col = c('black', 'magenta'))

par(mfrow=c(1,1))

# rho
round(test3_hier_standard$overall_rho, 3)
round(test3_hier_precondition$overall_rho, 3)

# Q
round(test3_hier_standard$overall_Q, 3)
round(test3_hier_precondition$overall_Q, 3)

save.image('hierarchical_exp_4_example.RData')
