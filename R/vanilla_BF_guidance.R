SH_T <- function() {
  # compute T in SH setting
  # compute variance of sub-posterior means (sigma_a^2)
  # compute lambda
  # compute k_1 (using desired threshold for N^{-1}CESS_{0})
  # give minimum T
}

SSH_T <- function() {
  # compute T in SSH setting
  # compute variance of sub-posterior means (sigma_a^2)
  # compute gamma = sigma_a^2
  # compute k_1 and k_2 (using desired threshold for N^{-1}CESS_{0})
  # also give the orders of k_1 and k_2
}

adaptive_mesh <- function() {
  # use full mesh guidance for step j
}

regular_mesh <- function() {
  # use regular mesh (don't compute average variation of trajectories)
}

# since the regular mesh can be computed prior to calling the algorithm,
# but adaptive mesh must be computed at each stage of the algorithm,
# then I need to have some sort of flag or logical variable which calls the adaptive
# mesh at each iteration of the algorithm. If adaptive==FALSE, we continue using
# the time mesh that is passed, otherwise, we compute one from the guidance.
# Maybe we can just set time_mesh = 'adaptive' to do this...

SH <- function() {
  
}

SSH <- function() {
  
}