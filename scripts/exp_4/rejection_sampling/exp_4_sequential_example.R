library(hierarchicalFusion)

######################################## examples ########################################

seed <- 21
set.seed(seed)

# rejection sampling target
# dominating_dnorm <- function(x) 1.35*dnorm(x, mean = 0, sd = 1)
# curve(dominating_dnorm(x), -5, 5, col = 'red')
# curve(dnorm(x, mean = 0, sd = 1), -5, 5, add = T, col = 'green')
# curve(target_density_exp_4(x), -5, 5, add = T, col = 'blue')
test_target_mc <- sample_exp_4(N = 100000,
                               proposal_mean = 0,
                               proposal_sd = 1,
                               dominating_M = 1.35,
                               beta = 1)

######################################## beta = 1/4

# sampling from the tempered targets
# finding the right proposal distribution and constant M in rejection sampling
# dominating_dnorm <- function(x) 1.4*dnorm(x, mean = 0, sd = 1.5)
# curve(dominating_dnorm(x), -5, 5, col = 'red')
# curve(dnorm(x, mean = 0, sd = 1.5), -5, 5, add = T, col = 'green')
# curve(tempered_target_density_exp_4(x, 0, 1/4), -5, 5, add = T, col = 'blue')

# using rejection sampling to obtain input samples
seed <- 21
set.seed(seed)
input_samples1 <- base_rejection_sampler_exp_4(beta = 1/4,
                                               nsamples = 100000,
                                               proposal_mean = 0,
                                               proposal_sd = 1.5,
                                               dominating_M = 1.4)
# curve(tempered_target_density_exp_4(x, beta = 1/4), -4, 4)
# # check the samples look okay
# for (samples in input_samples1) {
#   lines(density(samples), col = 'blue')
# }

# standard progressive fusion
test1_prog_standard <- progressive_fusion_exp_4(N_schedule = rep(100000, 3),
                                                time_schedule = rep(1, 3),
                                                base_samples = input_samples1, 
                                                mean = 0,
                                                start_beta = 1/4,
                                                precondition = FALSE,
                                                seed = seed)

# preconditioned progressive fusion
test1_prog_precondition <- progressive_fusion_exp_4(N_schedule = rep(100000, 3),
                                                    time_schedule = rep(1, 3),
                                                    base_samples = input_samples1, 
                                                    mean = 0,
                                                    start_beta = 1/4,
                                                    precondition = TRUE,
                                                    seed = seed)

curve(exp_4_density(x), -2.5, 2.5)
lines(density(test_target_mc), col = 'black')
lines(density(test1_prog_standard$samples[[1]]), col = 'green')
lines(density(test1_prog_precondition$samples[[1]]), col = 'blue')

# acceptance rate comparisons
acceptance_rate_plots(hier1 = test1_prog_precondition,
                      hier2 = test1_prog_standard,
                      time = 1,
                      hierarchical = FALSE)

######################################## beta = 1/8

# sampling from the tempered targets
# finding the right proposal distribution and constant M in rejection sampling
# dominating_dnorm <- function(x) 1.3*dnorm(x, mean = 0, sd = 1.5)
# curve(dominating_dnorm(x), -5, 5, col = 'red')
# curve(dnorm(x, mean = 0, sd = 1.5), -5, 5, add = T, col = 'green')
# curve(tempered_target_density_exp_4(x, 0, 1/8), -5, 5, add = T, col = 'blue')

# using rejection sampling to obtain input samples
seed <- 21
set.seed(seed)
input_samples2 <- base_rejection_sampler_exp_4(beta = 1/8,
                                               nsamples = 100000,
                                               proposal_mean = 0,
                                               proposal_sd = 1.5,
                                               dominating_M = 1.3)
# curve(tempered_target_density_exp_4(x, beta = 1/8), -4, 4)
# # check the samples look okay
# for (samples in input_samples2) {
#   lines(density(samples), col = 'blue')
# }

# standard progressive fusion
test2_prog_standard <- progressive_fusion_exp_4(N_schedule = rep(100000, 7),
                                                time_schedule = rep(1, 7),
                                                base_samples = input_samples2,
                                                mean = 0,
                                                start_beta = 1/8,
                                                precondition = FALSE,
                                                seed = seed)

# preconditioned progressive fusion
test2_prog_precondition <- progressive_fusion_exp_4(N_schedule = rep(100000, 7),
                                                    time_schedule = rep(1, 7),
                                                    base_samples = input_samples2,
                                                    mean = 0,
                                                    start_beta = 1/8,
                                                    precondition = TRUE,
                                                    seed = seed)

curve(exp_4_density(x), -2.5, 2.5)
lines(density(test_target_mc), col = 'black')
lines(density(test2_prog_standard$samples[[1]]), col = 'green')
lines(density(test2_prog_precondition$samples[[1]]), col = 'blue')

# acceptance rate comparisons
acceptance_rate_plots(hier1 = test2_prog_precondition,
                      hier2 = test2_prog_standard,
                      time = 1,
                      hierarchical = FALSE)

# #####
#
# par(mfrow=c(1,2))
#
# #####
# # prog
#
# #rho
# plot(1:length(test2_prog_standard$rho_acc), test2_prog_standard$rho_acc, ylim = c(0,1), col = 'black',
#      ylab = expression(rho), xlab = 'Level', pch = 4, cex = 0.5)
# lines(1:length(test2_prog_standard$rho_acc), test2_prog_standard$rho_acc, col = 'black')
# points(1:length(test2_prog_precondition$rho_acc), test2_prog_precondition$rho_acc, col = 'orange', pch = 4, cex = 0.5)
# lines(1:length(test2_prog_precondition$rho_acc), test2_prog_precondition$rho_acc, col = 'orange', lty = 2)
#
# # Q
# plot(1:length(test2_prog_standard$Q_acc), test2_prog_standard$Q_acc, ylim = c(0,1), col = 'black',
#      ylab = 'Q', xlab = 'Level', pch = 4, cex = 0.5)
# lines(1:length(test2_prog_standard$Q_acc), test2_prog_standard$Q_acc, col = 'black')
# points(1:length(test2_prog_precondition$Q_acc), test2_prog_precondition$Q_acc, col = 'orange', pch = 4, cex = 0.5)
# lines(1:length(test2_prog_precondition$Q_acc), test2_prog_precondition$Q_acc, col = 'orange', lty = 2)
#
# #####
#
# par(mfrow=c(1,1))
#
# # rho
# round(test2_prog_standard$rho_acc, 3)
# round(test2_prog_precondition$rho_acc, 3)
#
# # Q
# round(test2_prog_standard$Q_acc, 3)
# round(test2_prog_precondition$Q_acc, 3)

######################################## beta = 1/16

# sampling from the tempered targets
# finding the right proposal distribution and constant M in rejection sampling
# dominating_dnorm <- function(x) 1.4*dnorm(x, mean = 0, sd = 1.5)
# curve(dominating_dnorm(x), -5, 5, col = 'red')
# curve(dnorm(x, mean = 0, sd = 1.5), -5, 5, add = T, col = 'green')
# curve(tempered_target_density_exp_4(x, 0, 1/16), -5, 5, add = T, col = 'blue')

# using rejection sampling to obtain input samples
seed <- 21
set.seed(seed)
input_samples3 <- base_rejection_sampler_exp_4(beta = 1/16,
                                               nsamples = 100000,
                                               proposal_mean = 0,
                                               proposal_sd = 1.5,
                                               dominating_M = 1.4)
# curve(tempered_target_density_exp_4(x, beta = 1/16), -4, 4)
# # check the samples look okay
# for (samples in input_samples3) {
#   lines(density(samples), col = 'blue')
# }

# standard progressive fusion
test3_prog_standard <- progressive_fusion_exp_4(N_schedule = rep(100000, 15),
                                                time_schedule = rep(1, 15),
                                                base_samples = input_samples3,
                                                mean = 0,
                                                start_beta = 1/16,
                                                precondition = FALSE,
                                                seed = seed)

# preconditioned progressive fusion
test3_prog_precondition <- progressive_fusion_exp_4(N_schedule = rep(100000, 15),
                                                    time_schedule = rep(1, 15),
                                                    base_samples = input_samples3,
                                                    mean = 0,
                                                    start_beta = 1/16,
                                                    precondition = TRUE,
                                                    seed = seed)

curve(exp_4_density(x), -2.5, 2.5)
lines(density(test_target_mc), col = 'black')
lines(density(test3_prog_standard$samples[[1]]), col = 'green')
lines(density(test3_prog_precondition$samples[[1]]), col = 'blue')

# acceptance rate comparisons
acceptance_rate_plots(hier1 = test3_prog_precondition,
                      hier2 = test3_prog_standard,
                      time = 1,
                      hierarchical = FALSE)

# #####
#
# par(oma = c(2, 1, 0, 1))
# par(mfrow = c(1,2))
#
# #####
# # prog
#
# #rho
# plot(1:length(test3_prog_standard$rho_acc), test3_prog_standard$rho_acc, ylim = c(0,1), col = 'black',
#      ylab = expression(rho), xlab = 'Level', pch = 4, lwd = 1.5)
# lines(1:length(test3_prog_standard$rho_acc), test3_prog_standard$rho_acc, col = 'black', lwd = 1.5)
# points(1:length(test3_prog_precondition$rho_acc), test3_prog_precondition$rho_acc, col = 'magenta', pch = 2, lwd = 1.5)
# lines(1:length(test3_prog_precondition$rho_acc), test3_prog_precondition$rho_acc, col = 'magenta', lty = 2, lwd = 1.5)
#
# # Q
# plot(1:length(test3_prog_standard$Q_acc), test3_prog_standard$Q_acc, ylim = c(0,1), col = 'black',
#      ylab = 'Q', xlab = 'Level', pch = 4, lwd = 1.5)
# lines(1:length(test3_prog_standard$Q_acc), test3_prog_standard$Q_acc, col = 'black',lwd = 1.5)
# points(1:length(test3_prog_precondition$Q_acc), test3_prog_precondition$Q_acc, col = 'magenta', pch = 2, lwd = 1.5)
# lines(1:length(test3_prog_precondition$Q_acc), test3_prog_precondition$Q_acc, col = 'magenta', lty = 2, lwd = 1.5)
#
# par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
# plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# legend("bottom",
#        legend = c('standard', 'preconditioned'),
#        lty = c(1, 2),
#        xpd = TRUE,
#        horiz = TRUE,
#        inset = c(0, 0),
#        bty = "n",
#        pch = c(4, 2),
#        col = c('black', 'magenta'))
#
# par(mfrow=c(1,1))
#
# # rho
# round(test3_prog_standard$rho_acc, 3)
# round(test3_prog_precondition$rho_acc, 3)
#
# # Q
# round(test3_prog_standard$Q_acc, 3)
# round(test3_prog_precondition$Q_acc, 3)

save.image('progressive_exp_4_example.RData')
