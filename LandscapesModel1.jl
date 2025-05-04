#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Landscapes Model 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---------------------
# Packages
#---------------------
#using Pkg

#Pkg.add("Plots")
#Pkg.add("Plotly")
#Pkg.add("StatsBase")
#Pkg.add("Distributions")
#Pkg.add("ArgParse")
println("Loading packages")
using Distributions
using ArgParse
using DelimitedFiles
using Plots
using Random
using StatsBase
using DataFrames
using CSV
using Printf

#--------------------------
# Input Files -- DEPRECATED, INPUT VIA SHELL INSTEAD
#--------------------------
#parasource = "/home/ubuntu/sims/model1/ParameterDict1.jl" # For use with HPC cluster
#parasource = "/home/joseph/github/LandKlifModel/ParameterDict1.jl" # For use on local machine

#---------------------
# Objects
#---------------------
# Structure of patches
# Patches are objects with positional indices (x,y), temperature, precipitation, and species populations.
# Species populations are structured as nested arrays. The x index defines the species, the y index the individual,
# and the z index the trait value for different species traits
println("Defining structs")
struct TPatch
    row::Int # Row and column indices define patch location
    col::Int
    species::Array{Array{Float32,2},1} # Species defined by array index
    temp_t::Float32
    precip_t::Float32
    habitat::Float32
end

#----------------------------------------------
# Main Program
#----------------------------------------------

# Gets and returns info about parameter dictionary for diagnostic purposes
function dict_info(par::Dict)
    println("=========================")
    println("Parameter Dict Info")
    for (key, value) in par
        println(key, ": ",value,",", typeof(value))
    end
    println("=========================")
    println("")
end
println("compiling main sim function")
# MAIN SIMULATION
# Run simulation using landscape and parameters from input files. Requires inputs for temperature, environment,
# scen,grad... are loop indices for different scenarios (e.g. climate trend, gradient strength)
function simulation_run()
    println("Starting.")
    ArgDict = read_arguments() # Read in shell input
    println(ArgDict)
    parasource = ArgDict["parasource"]
    println(parasource)
    println(typeof(parasource))
    par = include(parasource) # See parameter config file (ParaDict) for parameter overview
    dict_info(par) # Return info on parameter dict
    # READ IN PARAMETERS FROM PARAMETER FILE. SEE EXAMPLE PARAMETER CONFIG FILE FOR OVERVIEW OF CONTENTS
    if ArgDict["seed"]==nothing # Checks if a seed is provided in the shell arguements (see shell scripts for example). If none is specified, the program will chekc if one is given in the parameter file.
        if haskey(par,"seed")==true# Checks of parasource includes a key for the RNG seed. If no seed is specified, one will be chosen based on the system clock
            Random.seed!(par["seed"])
            println("Seed: ",par["seed"])
        end
    else
        Random.seed!(ArgDict["seed"])
        println("Seed: ",ArgDict["seed"])
    end
    scen = par["scen"] # Climate scenario parameter
    bmax = par["bmax"] # Length of burnin period
    tmax = par["tmax"] # length of main simulation run
    rmax = par["rmax"] # Number of replicates. Keep at 1.
    grad = par["grad"] # Strength of compositional heterogeneity
    α = par["α"] # Trade-off term
    K = par["carry_capacity"]
    T_ref = par["T_ref"]
    immi = par["immi"] # Immigration? Boolean
    mut = par["mutate"] # Mutation? Boolean
    p_mut = par["p_mut"] # Mutation probability
    mut_sd = par["mut_sd"] # standard deviation of mutation
    mut_decay = par["mut_decay"] # decay over time of standard deviation of mutation
    p_immi = par["p_immi"] # Chance of immigration in a patch
    e_immi = par["e_immi"] # Expected immigrants in a patch
    dir = par["dir"] # Check output directory. If the directory does not exist, make it
    if isdir(dir)==false
        println("Directory not found. Creating output directory.")
        mkpath(dir)
    end
    println("parasource = $parasource")
    println("Climate scenario $scen")
    println("carry capacity = $K")
    println(dir)
    filename,tempsource,envsource = set_filename(ArgDict,par,scen,grad) # Set the names for the output files
    #Random.seed!(123)
    for rep in 1:par["rmax"]
        println("Replicate $rep")
        println("Initializing world.")
        init_world(tempsource,envsource,grad,par) # Read in landscape files and intialize the landscape
        trend, mean_trend = generate_climate_trend(tmax,0,1,scen) # Generate a trend for fluctuations in T for the main simulation period
        if haskey(par,"cmax")==true # Check if there is a climate shift period. If there is, read in the length of the period and rate of change and create a second T trend for that period
            trend2, mean_trend2 = generate_climate_trend2(par["cmax"],0,1,par["c_change"])
        end
        println(length(trend))
        init_spp2() # Legacy technical debt. Do not remove or the program won't work.
        println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
        write_landscape_stats(landscape,dir,filename,rep,-1,scen,0,0,grad,par["autocor_temp"],par["autocor_env"],α,bmax)
        write_landscape_csv(landscape,dir,filename,rep,-1,scen,0,0,grad,par["autocor_temp"],par["autocor_env"],α)
        if par["burnin"] == true
            println("Starting burn-in period.") # A burn-in period with no environmental fluctuation
            for b in 1:par["bmax"]
                if b==par["bmax"]
                    println("step =",step)
                end
                #println("burn-in time step: $b")
                burnin = true
                #println("Beginning burn-in timestep $b.")
                #println("Starting dispersal routine.")
                #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
                dispersal!(landscape)
                #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
                #println("Starting reproduction routine.")
                demographics(landscape,α,0,grad,K,b)
                #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
                if mut==true
                    mutate(landscape,p_mut,mut_sd,mut_decay,b)
                end
                write_landscape_stats(landscape,dir,filename,rep,b,scen,0,0,grad,par["autocor_temp"],par["autocor_env"],α,bmax)
                if haskey(par,"output_burnin")==true
                    if par["output_burnin"]==true
                        if b==par["bmax"]
                            write_landscape_csv(landscape,dir,filename,rep,b,scen,0,0,grad,par["autocor_temp"],par["autocor_env"],α)
                        end
                    end
                end
                #println("End timestep $t.")
            end
        end
        #write_landscape_csv(landscape,dir,filename,rep,0,scen,0,grad,par["autocor_t"],par["autocor_e"],α)
        println("Starting main simulation")
        for t in 1:par["tmax"]
            if par["burnin"] == true # For correctly counting the number of timesteps
                step=t+par["bmax"]
            else
                step=t
            end
            if t==par["tmax"]
                println("step =",step)
            end
            burnin = false
            #println("Beginning timestep $t.")
            #println("Starting dispersal routine.")
            #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
            dispersal!(landscape)
            #println("Starting reproduction routine.")
            if par["immi"]==true
                demographics_immi(landscape,par,α,trend[t],K,e_immi,step)
            else
                demographics(landscape,α,0,grad,K,step)
            end
            #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
            if mut==true
                mutate(landscape,p_mut,mut_sd,mut_decay,step)
            end
            write_landscape_stats(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α,bmax)
            if haskey(par,"output_time")==true
                if haskey(par,"output_interval")==true
                    if mod((step-par["output_time"]),par["output_interval"])==0 || step==par["output_time"]
                        write_landscape_csv(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α)
                    end
                else
                    if step==par["output_time"]
                        write_landscape_csv(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α)
                    end
                end
            else
                if step==par["tmax"]
                    write_landscape_csv(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α)
                end
            end
            #if haskey(par,"output_interval")==true
            #    if haskey(par,"output_start")==true
            #        if mod(step,par["output_interval"])==0 || step>=par["output_start"] || t==tmax
            #            write_landscape_csv(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α)
            #        end
            #    end
            #elseif haskey(par,"output_interval")==false
            #    if haskey(par,"output_start")==true
            #        if step==par["output_start"]
            #            write_landscape_csv(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α)
            #        end
            #    end
            #else
            #    if t==tmax
            #        write_landscape_csv(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α)
            #    end
            #end
            #println("End timestep $t.")
        end
        if haskey(par,"cmax")==true
            println("Starting climate trend period.")
            for c in 1:par["cmax"]
                #println("Trend time step: $c")
                if par["burnin"] == true
                    step=par["tmax"]+par["bmax"]+c
                    tmax2=par["tmax"]+par["bmax"]+par["cmax"]
                    #print("step=",step)
                else
                    step=par["tmax"]+c
                    tmax2=par["tmax"]+par["cmax"]
                    #print("step=",step)
                end
                #println("Beginning burn-in timestep $b.")
                #println("Starting dispersal routine.")
                #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
                dispersal!(landscape)
                #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
                #println("Starting reproduction routine.")
                #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
                if par["immi"]==true
                    demographics_immi(landscape,par,α,trend2[c],K,e_immi,step)
                else
                    demographics(landscape,α,trend2[c],grad,K,step)
                end
                if mut==true
                    mutate(landscape,p_mut,mut_sd,mut_decay,step)
                end
                write_landscape_stats(landscape,dir,filename,rep,step,scen,trend2[c],mean_trend2[c],grad,par["autocor_temp"],par["autocor_env"],α,bmax)

                if haskey(par,"output_time")==true
                    if haskey(par,"output_interval")==true
                        if mod((step-par["output_time"]),par["output_interval"])==0 || step==par["output_time"]
                            write_landscape_csv(landscape,dir,filename,rep,step,scen,trend2[c],mean_trend2[c],grad,par["autocor_temp"],par["autocor_env"],α)
                        end
                    else
                        if step==par["output_time"]
                            write_landscape_csv(landscape,dir,filename,rep,step,scen,trend2[c],mean_trend2[c],grad,par["autocor_temp"],par["autocor_env"],α)
                        end
                    end
                else
                    if step==tmax2
                        write_landscape_csv(landscape,dir,filename,rep,step,scen,trend2[c],mean_trend2[c],grad,par["autocor_temp"],par["autocor_env"],α)
                    end
                end

                #if haskey(par,"output_interval")==true
                #    if haskey(par,"output_start")==true
                #        if mod(step,par["output_interval"])==0 || step>=par["output_start"] || c==par["cmax"]
                #            write_landscape_csv(landscape,dir,filename,rep,step,scen,trend2[c],mean_trend2[c],grad,par["autocor_temp"],par["autocor_env"],α)
                #        end
                #    end
                #elseif haskey(par,"output_interval")==false
                #    if haskey(par,"output_start")==true
                #        if step==par["output_start"]
                #            write_landscape_csv(landscape,dir,filename,rep,step,scen,trend2[c],mean_trend2[c],grad,par["autocor_temp"],par["autocor_env"],α)
                #        end
                #    end
                #else
                #    if c==par["cmax"]
                #        write_landscape_csv(landscape,dir,filename,rep,step,scen,trend2[c],mean_trend2[c],grad,par["autocor_temp"],par["autocor_env"],α)
                #    end
                #end
                #println("End timestep $t.")
            end
        end
        popcount(landscape)
        println("End time loop.")
    end
    println("Calculating patch stats.")
    patchstats(dir,filename,tempsource,envsource) # Calculate patch statistics using individual data
    println("Program complete.")
end

#----------------------------------------------
# Main Script
#----------------------------------------------
# s_clim is the climate scenario, s_grad is the gradient scenario
#simulation_run2(parasource,0)

@time begin
    # Load functions
    println("Loading Initialization.jl")
    include("Initialization.jl")

    println("Loading SimFunctions.jl")
    include("SimFunctions.jl")

    println("Loading Calculations.jl")
    include("Calculations.jl")

    println("Loading Output.jl")
    include("Output.jl")

    println("Loading PatchStatistics.jl")
    include("PatchStatistics.jl")
    
    println("Loading Misc.jl")
    include("Misc.jl")
    
    simulation_run()
end
