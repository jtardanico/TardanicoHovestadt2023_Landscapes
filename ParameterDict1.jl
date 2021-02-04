Dict{String,Any}(
"bmax" => 50, # Number of timesteps in burn-in period
"tmax" => 1000, #450, # Number of timesteps in Simulation_Run
"rmax" => 10, # Number of replicates
"T_ref" => 12.5, # DEPRECATED
"grad" => 1, # Gradient steepness
"α" => 3, # Trade-off strength
"initial_pop" => 10,
"topt_µ" => 0.2171483, # mean of T_opt initialization distribution
"topt_σ" => 0.7171647, # sd of T_opt initialization distribution
"hopt_µ" => 0.09578733, # mean of H_opt initialization distribution
"hopt_σ" => 0.7809819, # sd of H_opt initialization distribution
"tsd_µ" => 0.6117247, # µ of T_sd initialization distribution
"tsd_σ" => 0.1700659, # σ of T_sd initialization distribution
"hsd_µ" => -0.1295957, # µ of H_sd initialization distribution
"hsd_σ" => 0.4754802, # σ of H_sd initialization distribution
"immi" => false, # Enables or disables immigration
"mutate" => false, # Enables or disables mutation
"mut_sd" => 0.05, # Standard deviation of mutation trait change
"mut_decay" => 0, # Exponential decay constant for mutation sd
"p_immi" => 0.005, # Probability of a patch receiving immigrants
"e_immi" => 2.5, # Mean new immigrants
#"dir" => "/home/ubuntu/outputfiles/eqtest_mut_/",
"dir" => "/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/debug_testing18/",
)
