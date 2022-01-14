# # guidance for T
# # uses guidance to choose T \geq C^{3/2}/m and then we look at multiples of that
# # pass in a required ESS for initialise the algorithm
# # We keep choosing multiples of that until we reach our required ESS
# 
# # ---------- first importance sampling step
# particles <- rho_IS_multivariate(particles_to_fuse = particles_to_fuse,
#                                  dim = dim,
#                                  N = N,
#                                  m = m,
#                                  time = time_mesh[length(time_mesh)],
#                                  inv_precondition_matrices = rep(list(diag(1, dim)), m),
#                                  inverse_sum_inv_precondition_matrices = inverse_sum_matrices(rep(list(diag(1, dim)), m)),
#                                  number_of_steps = length(time_mesh),
#                                  time_mesh = time_mesh,
#                                  resampling_method = resampling_method,
#                                  seed = seed,
#                                  n_cores = n_cores,
#                                  cl = cl)
# particle_set$CESS[1]
# 
# 
# # guidance for P and n