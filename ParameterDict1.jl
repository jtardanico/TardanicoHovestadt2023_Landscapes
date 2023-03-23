# EXAMPLE PARAMETER CONFIGURATION FILE

Dict{String,Any}(
#-------------------------------------
# General Simulation parameters
"burnin" => false, # Burn-in period? Boolean
"bmax" => 50, # Number of timesteps in burn-in period
"tmax" => 10000, # Number of timesteps in main Simulation Run
"rmax" => 1, # Number of replicates. Will repeat simulation with the same landscape. Keep at 1 and do separate simulation runs.
"cmax" => 200, # Number of timesteps in the climate change period
"c_change" => 0.02, # Rate of change in T
"T_ref" => 12.5, # This is unnecessary for the actual model. Remove this.
"grad" => 0.1, # Gradient steepness
"α" => 3, # Trade-off strength
"scen" => 1, # Climate scenario. DEFUNCT PARAMETER. KEEP AT 1
"initial_pop" => 10, # Initial number of individuals in a patch
"carry_capacity" => 150, # Carring capacity of a patch
"immi" => true, # Is there immigration from outside the landscape? True or false.
"mutate" => false, # Do organisms undergo mutation. True or false
"p_mut" => 0.001, # Probability of mutation.
"mut_sd" => 0.05, # Standard deviation of trait change from mutation.
"mut_decay" => 0, # Change in mutation sd over time. Set to 0 for constant mutation sd.
"p_immi" => 0.005, # Chance of immigration to a patch.
"e_immi" => 2.5, # Expected immigrants in a patch
"seed" => 1234, # RNG seed
"output_time" => 10000, # Specify when the program should output data
#"output_interval" => 50, # Timestep interval for output of data for individual organisms. Program will only output once if this is not defined.
#-------------------------------------
# Initial trait distribution parameters
"topt_µ" => 0, # mean of T_opt initialization distribution
"topt_σ" => 1, # sd of T_opt initialization distribution
"hopt_µ" => 0, # mean of H_opt initialization distribution
"hopt_σ" => 1, # sd of H_opt initialization distribution
"tsd_µ" => 0, # µ of T_sd initialization distribution
"tsd_σ" => 1, # σ of T_sd initialization distribution
"hsd_µ" => 0, # µ of H_sd initialization distribution
"hsd_σ" => 1, # σ of H_sd initialization distribution
#--------------------------------------
# Filename and output directory parameters
"region" => "", # Name of region. Leave as "" if not using real world landscape data.
"temp_type" => "autocorrelated", # Type of landscape. Uniform, clustered, or autocorrelated
"env_type" => "autocorrelated", # Type of landscape.
"autocor_temp" => 0, # Autocorrelation degree (Hurst index) or cluster size of temperature dimension.
"autocor_env" => 0, # Autocorrelation degree (Hurst index) or cluster size of habitat dimension.
"dir" => "/home/ubuntu/output/series1/test_weak_gradient3_lowAC_trend/", # Output directory
)
