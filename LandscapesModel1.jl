
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
println("compiling main sim function")
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
    if haskey(par,"seed")==true# Checks of parasource includes a key for the RNG seed
        Random.seed!(par["seed"])
    end
    scen = par["scen"]
    bmax = par["bmax"]
    tmax = par["tmax"]
    rmax = par["rmax"]
    grad = par["grad"]
    α = par["α"]
    K = par["carry_capacity"]
    T_ref = par["T_ref"]
    immi = par["immi"]
    mut = par["mutate"]
    p_mut = par["p_mut"]
    mut_sd = par["mut_sd"]
    mut_decay = par["mut_decay"]
    p_immi = par["p_immi"]
    e_immi = par["e_immi"]
    dir = par["dir"]
    println("parasource = $parasource")
    println("Climate scenario $scen")
    println(dir)
    filename,tempsource,envsource = set_filename(ArgDict,par,scen,grad)
    #Random.seed!(123)
    for rep in 1:rmax
        println("Replicate $rep")
        println("Initializing world.")
        init_world(tempsource,envsource,grad,par)
        generate_climate_trend(tmax,0,1,scen)
        println(length(trend))
        init_spp2()
        println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
        write_landscape_stats(landscape,dir,filename,rep,-1,scen,0,0,grad,par["autocor_temp"],par["autocor_env"],α,bmax)
        write_landscape_csv(landscape,dir,filename,rep,-1,scen,0,0,grad,par["autocor_temp"],par["autocor_env"],α)
        if par["burnin"] == true
            println("Starting burn-in period.")
            for b in 1:bmax
                println("burn-in time step: $b")
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
                if mod(b,50)==true || b==bmax
                    write_landscape_csv(landscape,dir,filename,rep,b,scen,0,0,grad,par["autocor_temp"],par["autocor_env"],α)
                end
                #println("End timestep $t.")
            end
        end
        #write_landscape_csv(landscape,dir,filename,rep,0,scen,0,grad,par["autocor_t"],par["autocor_e"],α)
        println("Starting main simulation")
        for t in 1:tmax
            if par["burnin"] == true
                step=t+bmax
            else
                step=t
            end
            burnin = false
            println("Beginning timestep $t.")
            #println("Starting dispersal routine.")
            #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
            dispersal!(landscape)
            #println("Starting reproduction routine.")
            #println("Patch 1,1 pop: ",length(landscape[1,1].species[1][1:end,1]))
            demographics_immi(landscape,α,trend[t],grad,K,e_immi,step)
            if mut==true
                mutate(landscape,p_mut,mut_sd,mut_decay,step)
            end
            write_landscape_stats(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α,bmax)
            if mod(t,50)==0 && t>=9900 || t==tmax
                write_landscape_csv(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,par["autocor_temp"],par["autocor_env"],α)
            end
            #println("End timestep $t.")
        end
        popcount(landscape)
        println("End time loop. Analyzing landscape.")
        #diversity_analysis(landscape)
        #trait_analysis(landscape)
        #env_analysis(landscape,trend[tmax])
        #heatmaps(richness,temperatures,habitat,pops,trait_means,filename,dir)
        #write_landscape_csv(landscape,dir,filename,rep,tmax,scen,trend[tmax],grad,ArgDict["autocor_t"],ArgDict["autocor_e"],α)
    end
    patchstats(dir,filename,tempsource,envsource)
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
    #include("Misc.jl")
    simulation_run()
end
