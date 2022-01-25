SH_T <- function() {
  # compute T in SH setting
}

SSH_T <- function() {
  # compute T in SSH setting
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

vanilla_SH <- function() {
  
}

vanilla_SSH <- function() {
  
}