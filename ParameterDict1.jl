Dict{String,Any}(
"bmax" => 50, # Number of timesteps in burn-in period
"tmax" => 1000, #450, # Number of timesteps in Simulation_Run
"rmax" => 10, # Number of replicates
"T_ref" => 12.5, # This is unnecessary for the actual model. Remove this.
"grad" => 1, # Gradient steepness
"α" => 3, # Trade-off strength
"initial_pop" => 10,
"immi" => false,
"mutate" => false,
"p_immi" => 0.005,
"e_immi" => 2.5,
"dir" => "/home/ubuntu/outputfiles/test1/",
#"dir" => "/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/debug_testing18/",
)
