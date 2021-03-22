
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

struct TPatch
    row::Int # Row and column indices define patch location
    col::Int
    species::Array{Array{Float32,2},1} # Species defined by array index
    temp_t::Float32
    precip_t::Float32
    habitat::Float32
end

#---------------------
# Calculation Functions
#---------------------

# Calculates an individual's average number of offspring.
function expected_fert(T_opt,T_sd,H_opt,H_sd,Fert_max, temp_t, habitat, α_t, α_h, trend)
    temp = temp_t + trend # Patch temperature
    e_fert = Fert_max * exp(-(temp-T_opt)^2/(2*T_sd^2)) * exp(-(habitat-H_opt)^2/(2*H_sd^2)) * exp(-T_sd^2/(2*α_t^2)) * exp(-H_sd^2/(2*α_h^2))
    return e_fert
end

# Calculates expected offspring. For use in calculating landscape average fitness (see write_landscape_stats)
function expected_fert2(T_opt,T_sd,H_opt,H_sd,temp_t, habitat, α_t, α_h, trend)
    #println("T_opt",size(T_opt))
    #println("T_sd",size(T_sd))
    #println("H_opt",size(H_opt))
    #println("H_sd",size(H_sd))
    #println("temp_t",size(temp_t))
    #println("habitat",size(habitat))
    #println("α",size(α_t))
    #println("trend",size(trend))
    temp = temp_t + trend # Patch temperature
    t_fert = exp(-(temp-T_opt)^2/(2*T_sd^2)) * exp(-T_sd^2/(2*α_t^2))
    h_fert = exp(-(habitat-H_opt)^2/(2*H_sd^2)) * exp(-H_sd^2/(2*α_h^2))
    e_fert = t_fert * h_fert
    return t_fert,h_fert,e_fert
end

# Function for calculating expected offspring mortality
# Terms: max fertility (Fert_max), number of offspring (offspring), carrying capacity (k).
function expected_mort(Fert_max, offspring, k)
    a = (Fert_max - 1)/(k * Fert_max)
    e_mort = 1/(1+a*offspring)
    return e_mort
end

# Function for calculating how stressful the environment of a patch is
function stress(T_opt::Float32,T_sd::Float32,H_opt::Float32,H_sd::Float32,temp_t::Float32, habitat::Float32,trend::Float32)
    temp = temp_t + trend
    S_T = exp(-(temp-T_opt)^2/(2*T_sd^2))
    S_H = exp(-(habitat-H_opt)^2/(2*H_sd^2))
    return S_T, S_H
end


# Determines carry capacity for a species; may be removed/replaced later.
function carry_capacity(k_0, S_T, S_H)
    K = round(k_0*S_T*S_H)  #K = round(k_0*S_T*S_H)
    return K
end

# Function for creating an inverse gaussian distribution for landscape cell temperature across x and y dimensions
function gaussian_landscape(row::Int,col::Int, max, mean, sd)
    patch_temp = -max * exp(-(row-mean)^2/(2*sd^2)) * exp(-(col-mean)^2/(2*sd^2)) + (max/2)
    return patch_temp
end

# Calculates weighted mean
function weightedmean(means,weights)
    wm = sum(means .* weights)/sum(weights)
    return wm
end

function shannon(i)
    p = i ./ sum(i[1:length(i)])
    h = -sum(p.*log.(p))
    return h
end

function simpson(i)
    p = i ./ sum(i[1:length(i)])
    λ =  sum(p.^2)
    return λ
end

# Exponential decay function. Rate should be a number between 0 and 1.
function decay(initial,rate,time)
    initial * exp(-(rate)*time)
end

#-------------------------
# Data Input Functions
#-------------------------

# Reads landscape temperature and environment data from files.
# Called by: init_world
function read_world_source(worldtempsource::String,worldenvsource::String)
    worldtempfile = readdlm(worldtempsource, ';')
    worldenvfile = readdlm(worldenvsource, ';')
    if length(worldtempfile[1:end,1]) == length(worldenvfile[1:end,1]) && length(worldtempfile[1,1:end]) == length(worldenvfile[1,1:end])
        landscapetemp = Float32.(worldtempfile)
        landscapeenv = Float32.(worldenvfile)

        return landscapetemp,landscapeenv
    else
        error("Dimension of landscape source files do not match.")
    end
end

# For reading in input from a shell script
function read_arguments()
    s=ArgParseSettings()
    @add_arg_table s begin
        "--parasource","-n"
            help = "parameter dictionary source file"
        "--tempsource", "-t" # Name of temperature source file in the shell script
            help = "Landscape temperature source file"
        "--envsource", "-e"
            help = "Landscape environment source file"
        "--scenario", "-s"
            help = "Climate scenario. See documentation"
            arg_type = Int
            default = Int(0)
        "--burninperiod", "-b" # Burn in period before simulation start
            help = "Burn-in period before simulation start"
            arg_type = Bool
            default = false
        "--autocorrelatedtemp", "-a" # Autocorrelated or clustered landscape? Used for output file name
            help = "Is the tempsource autocorrelated (true) or clustered (false)?"
            arg_type = Bool
            default=false
        "--autocor_t", "-c"
            help = "Degree of autocorrelation (H) or cluster size for temperature"
            arg_type = Float32
            default=Float32(NaN)
        "--autocorrelatedenv", "-i" # ^ Ditto
            help = "Is the envsource autocorrelated (true) or clustered (false)?"
            arg_type = Bool
            default=false
        "--autocor_e", "-d"
            help = "Degree of autocorrelation (H) or cluster size for environment"
            arg_type = Float32
            default=Float32(NaN)
        "--uniformt", "-u"
            help = "Is temperature uniform throughout the landscape?"
            arg_type = Bool
            default = false
        "--uniformenv", "-o"
            help = "Is the environment uniform throughout the landscape?"
            arg_type = Bool
            default = false
        end
        return parse_args(s)
end
#-------------------------
# Simulation Functions
#-------------------------

# Creates an array with dimentions n_species and n_traits which acts as a master list of all species in the simulation
function init_spp(n_species::Int)
    #println("Initializing species master list...")
    n_traits = 9
    global species_list = Array{Any,2}(undef,n_species,n_traits)
    for i in 1:n_species
        ID = i              # Trait 1: Species ID number
        T_opt = rand(Normal(mean,sd))*20+10 #8+22 # Trait 2: Temperature optimum
        T_sd = 2            # Trait 3: Temperature tolerance
        H_opt = rand()      # Trait 4: Habitat optimum
        H_sd = 0.2          # Trait 5: Habitat tolerance
        Disp_l = 0.1        # Trait 6: Local dispersal probability
        Disp_g = 0.01       # Trait 7: Global dispersal probability
        Fert_max = 10       # Trait 8: Maximum number of offspring
        dispersed = false   # Trait 9: Whether or not individual has already dispersed
        global species_list[i,1:end] = [ID, T_opt, T_sd, H_opt, H_sd, Disp_l, Disp_g, Fert_max, dispersed]
    end
end

# Creates a 1x9 array listing species traits. For use with function "simulation_run2"
function init_spp2()
    #println("Initializing species master list...")
    n_traits = 10
    global species_list = Array{Any,2}(undef,1,n_traits)
        ID = 1              # Trait 1: Species ID number
        T_opt = "10-30" #8+22 # Trait 2: Temperature optimum
        T_sd = "random"            # Trait 3: Temperature tolerance
        H_opt = "random"      # Trait 4: Habitat optimum
        H_sd = "random"          # Trait 5: Habitat tolerance
        Disp_l = "random"        # Trait 6: Local dispersal probability
        Disp_g = "random"       # Trait 7: Global dispersal probability
        Fert_max = "5-30"       # Trait 8: Maximum number of offspring
        dispersed = "N/A"   # Trait 9: Whether or not individual has already dispersed
        lineage = "N/A"    # Trait 10: Lineage identifier
        global species_list[1,1:end] = [ID, T_opt, T_sd, H_opt, H_sd, Disp_l, Disp_g, Fert_max, dispersed, lineage]
end

# Initialization of species populations
# Function creates a 1D array of 2D arrays of length n_species,
# each with dimensions n_pop by n_traits.
# Each row in a 2D array is an individual of a species; each 2D array is a population.
function init_pops(n_pop::Int)
    #println("Initializing species populations...")
    n_species = length(species_list[1:end,1])
    n_traits = length(species_list[1,1:end])
    global species_0 = Array{Array{Float32,2},1}(undef, n_species)
    for i in 1:n_species
        n_pop = 10 # Placeholder value
        population = Array{Float32,2}(undef,n_pop,n_traits)
        for j in 1:n_pop
            population[j,1:end] = species_list[i,1:end]
        end
        global species_0[i] = population
    end
end

# Initialized a population of individuals with randomized trait values
function init_popgrad(n_pop::Int,par::Dict)
    #println("init_popgrad")
    n_traits = 10
    patchpop = Array{Array{Float32,2},1}(undef,1)
    population = Array{Float32,2}(undef,n_pop,n_traits)
    for i in 1:n_pop
        ID = i              # Trait 1: Species ID number
        T_opt = rand(Normal(par["topt_µ"],par["topt_σ"])) #8+22 # Trait 2: Temperature optimum
        T_sd = rand(LogNormal(par["tsd_µ"],par["tsd_σ"]))    # Trait 3: Temperature tolerance
        H_opt = rand(Normal(par["hopt_µ"],par["hopt_σ"]))    # Trait 4: Habitat optimum
        H_sd = rand(LogNormal(par["hsd_µ"],par["hsd_σ"])) # Trait 5: Habitat tolerance
        Disp_l = rand()      # Trait 6: Dispersal probability
        Disp_g = rand()       # Trait 7: Global dispersal probability
        Fert_max = 15       # Trait 8: Maximum number of offspring
        dispersed = false   # Trait 9: Whether or not individual has already dispersed
        lineage = rand(Float32) # Trait 10: Lineage identifier
        population[i,1:end] = [ID, T_opt, T_sd, H_opt, H_sd, Disp_l, Disp_g, Fert_max, dispersed, lineage]
    end
    patchpop[1] = population
    return patchpop
end

# Generates a temperature trend. Several trend types are possible:
# const: Constant mean with random variation
# rise: Rising trend with no random variation
# rise2: Rising mean with random variation
# none: Default trend if no other is specified. Constant with no random variation
function generate_climate_trend(nsteps::Int,mean,sd,trend_type::Int)
    #println("generate_climate_trend")
    #println("nsteps=$nsteps")
    if trend_type == 1
        global trend = rand(Normal(mean,sd),nsteps) # Random variation, constant mean
        global mean_trend = zeros(nsteps) # mean trend
    elseif trend_type == 2
        t = Array{Float32,1}(undef,nsteps)
        for i in 1:nsteps
            t[i] = 0.025*i
        end
        global trend = t
        global mean_trend = t # mean trend
    elseif trend_type == 3
        t = Array{Float32,1}(undef,nsteps)
        for i in 1:nsteps
            t[i] = 0.025*i + rand(Normal(mean,sd))
            t_mean =  0.025*i # mean trend
        end
        global trend = t # Random variation & rising mean
        global mean_trend = t_mean
    elseif trend_type == 0
        global trend = zeros(nsteps)
        global mean_trend = zeros(nsteps)
    else
        println("Invalid scenario. Valid scenarios are:")
        println("0: Constant mean, no random variation")
        println("1: Constant mean temperature, constant variance")
        println("2: Increasing mean temperature, no random variation")
        println("3: Increasing mean, constant variance")
        exit(code=0)
    end
end

# Initializes a landscape from file inputs.
# 'gradient' defines landscape gradient steepness
function init_world(worldtempsource::String,worldenvsource::String,gradient,parameters::Dict)
    #println("init_world")
    temperatures, environments = read_world_source(worldtempsource,worldenvsource)
    temperatures = temperatures * gradient
    environments = environments * gradient
    nrows = length(temperatures[1:end,1])
    ncols = length(temperatures[1,1:end])
    m_t, v_t = mean_and_var(temperatures)
    m_e, v_e = mean_and_var(environments)
    global landscape = Array{TPatch, 2}(undef, nrows, ncols)
    for i in 1:nrows
        for j in 1:ncols
            row = i
            col = j
            patchpop = init_popgrad(100,parameters)
            #if (i==1 ) && (j==1)
            #    patchpop = init_popgrad(100,0.5,0.5)
            #else
            #    patchpop = Array{Array{Float32,2},1}(undef,1)
            #    patchpop[1] = Array{Float32,2}(undef,0,10)
            #end
            temp = temperatures[i,j]
            habitat = environments[i,j]
            precip = 0
            patch = TPatch(row,col,patchpop,temp,precip,habitat)
            landscape[i,j] = patch
        end
    end
end

# Generates a landscape with dimensions dim_x and dim_y and adds cell index values
# as well as temperature and precipitation values. For testing purposes.
function generate_world(dim_x::Int,dim_y::Int)
    #println("Creating landscape...")
    global landscape = Array{TPatch, 2}(undef, dim_x,dim_y)
    for i in 1:dim_x
        for j in 1:dim_y
            row = i
            col = j
            init_pops(10)          # Placeholder value; will change later
            temp = gaussian_landscape(i,j,25,25,12)  # Placeholder function; Gaussian function of x and y dims w/ mean at i,j=50,50
            precip = 500  # ^ Ditto
            habitat = rand()
            patch = TPatch(row,col,species_0,temp,precip,habitat)
            landscape[i,j] = patch
        end
    end
end

# Makes random changes to the traits of each individual in the landscape. Standard deviation of the changes is provided by mutsd.
function mutate(landscape::Array{TPatch,2},p_mut,mut_sd,mut_decay,timestep)
    function capped_poisson(mean,limit)
        p = rand(Poisson(mean))
        if p > limit
            #println("p > limit")
            return limit
        else
            #println("p = $p")
            return p
        end    
    end
    #println("Mutating")
    mut_t = decay(mut_sd,mut_decay,timestep)
    for i in 1:length(landscape[1:end,1]) # loop over landscpae length
        for j in 1:length(landscape[1,1:end]) # loop over landscape width
            if length(landscape[i,j].species[1:end])>0
                for p in 1:length(landscape[i,j].species[1:end]) # loop over species
                    if length(landscape[i,j].species[p][1:end,1])>0
                        len = length(landscape[i,j].species[p][1:end,1]) # Get length of individuals array
                        #println("len = $len")
                        #println(typeof(len))
                        #println("n = $n")
                        #println(typeof(n))
                        indices = sample(collect(1:1:len;),capped_poisson((p_mut*len),len);replace=false) # Get random sample of array indices. Sample n is a random Poisson number capped at pop size
                        #println(typeof(indices))
                        #for q in length(landscape[i,j].species[p][1:end,1]) # loop over individuals
                        for q in 1:length(indices) # Loop over indices
                            #println("Mutating at index ",indices[q])
                            #println("rand normal =", typeof(rand(Normal(0,mut_t))))
                            #println("trait =",landscape[i,j].species[p][indices[q],2])
                            landscape[i,j].species[p][indices[q],2] = landscape[i,j].species[p][indices[q],2] .+ rand(Normal(0,mut_t))# Temperature optimum
                            landscape[i,j].species[p][indices[q],3] = landscape[i,j].species[p][indices[q],3] .* (rand(LogitNormal(0,mut_t)))# Temperature tolerance
                            #println("tsd = $(landscape[i,j].species[p][indices[q],3])")
                            landscape[i,j].species[p][indices[q],4] = landscape[i,j].species[p][indices[q],4] .+ rand(Normal(0,mut_t))# Habitat optimum
                            #println("hopt = $(landscape[i,j].species[p][indices[q],4])")
                            landscape[i,j].species[p][indices[q],5] = landscape[i,j].species[p][indices[q],5] .* (rand(LogitNormal(0,mut_t)))# Habitat tolerance
                            #println("hsd = $(landscape[i,j].species[p][indices[q],5])")
                             d = landscape[i,j].species[p][indices[q],6] .+ rand(Normal(0,mut_t))# Dispersal chance
                             #println("d = ",typeof(d))
                             #println("d=$d")
                             if d .< 0
                                 d = 0
                                 #println("d set to $d")
                             elseif d .> 1
                                 d = 1
                                 #println("d set to $d")
                             end
                             landscape[i,j].species[p][indices[q],7] = copy(d)
                             #println("displ = $(landscape[i,j].species[p][indices[q],6])")
                             d = landscape[i,j].species[p][indices[q],7] .+ rand(Normal(0,mut_t)) # Global dispersal
                             if d .< 0
                                 d = 0
                                 #println("d set to $d")
                             elseif d .> 1
                                 d = 1
                                 #println("d set to $d")
                             end
                             landscape[i,j].species[p][indices[q],7] = copy(d)
                             #println("dispg = $(landscape[i,j].species[p][q,7])")
                        end # end individuals loop
                    end  # end individual pop length check
                end # end species loop
            end # end species length check
        end # end landscape width loop
    end # end landscape length loop
end # end function


# Calculates number of offspring for each individual. Offspring form the next generation of individuals.
function demographics(landscape::Array{TPatch, 2},niche_tradeoff, trend, grad, K::Int,burnin, immi,p_immi,e_immi)
    #println("demographics")
    for i in 1:length(landscape[1:end,1]) # begin landscape length loop
        for j in 1:length(landscape[1,1:end]) # Begin landscape width loop
            offspring = Array{Array{Int,1},1}(undef,length(landscape[i,j].species[1:end])) # List of offspring for all species
            total_offspring = 0 # Tally of all offpring of all species
            for p in 1:length(landscape[i,j].species[1:end]) # Loop over species
                sp_pop = length(landscape[i,j].species[p][1:end,1])
                #println("Species $p pop. = $sp_pop")
                if length(landscape[i,j].species[p][1:end,1])>0 # Check species population size
                    #println("YEET!")
                    expected = expected_fert.(landscape[i,j].species[p][1:end,2],
                               landscape[i,j].species[p][1:end,3],
                               landscape[i,j].species[p][1:end,4],
                               landscape[i,j].species[p][1:end,5],
                               landscape[i,j].species[p][1:end,8],
                               landscape[i,j].temp_t,
                               landscape[i,j].habitat, niche_tradeoff,
                               niche_tradeoff, trend) # Calculate expected offpring, see init_spp for trait key
                    #println("expected offspring = $expected")
                    global x = expected # Move to earlier position in loop structure.
                    offspring[p] = rand.(Poisson.(Float64.(expected))) # Calculate offspring produced
                    #println("Offspring = $offspring")
                    species_offspring = sum(offspring[p]) # For diagnostic purposes. Will be removed once testing is finished
                    #println("Total offspring of species $p = $species_offspring")
                    total_offspring = total_offspring + species_offspring
                end
            end # End loop over species
            for p in 1:length(landscape[i,j].species[1:end]) # Begin second loop over species
                if length(landscape[i,j].species[p][1:end,1])>0 # Check species population size
                    #=S_T, S_H = stress(landscape[i,j].species[p][1,2], # Calculating environmental stress for a species in a patch
                               landscape[i,j].species[p][1,3],
                               landscape[i,j].species[p][1,4],
                               landscape[i,j].species[p][1,5],
                               landscape[i,j].temp_t,
                               landscape[i,j].habitat, trend)=#
                    #println("S_T = $S_T, S_H = $S_H")
                    #K = carry_capacity(300,S_T,S_H) # Calculates carrying capacity of a patch
                    #println("K = $K")
                    p_surv = expected_mort.(landscape[i,j].species[p][1:end,8],total_offspring,K) # Calculate expected mortality of oppspring
                    if K > 0 # Note: May be unecessary in future model iterations. Consider removal.
                        surviving = rand.(Binomial.(offspring[p], Float64.(p_surv))) # Calculate surviving offspring
                    else
                        surviving = 0
                    end
                    #println("Species $p has $surviving surviving offspring.")
                    if sum(surviving) > 0 # Check number of survivors
                        if burnin==false && immi==true # Check burn in and immi conditions
                            immigrants = rand(Poisson(e_immi))
                            lenx = sum(surviving) + immigrants # make e_immi a dictionary parameter
                            newgen = Array{Float32,2}(undef,lenx,length(landscape[i,j].species[p][1,1:end])) # Creates an array of length sum(offspring) + new immigrants with data for species p & width of n traits
                            ind = 1 # Keeps count of individual offspring added to newgen array
                            for q in 1:length(surviving) # Goes down index of 'surviving'
                                if surviving[q] > 0
                                    for r in 1:surviving[q] # Loop from 1 to number of surviving offspring
                                        newgen[ind,1:end] = landscape[i,j].species[p][q,1:end]  #
                                        newgen[ind,9] = false
                                        #println("$(newgen[ind,9])")
                                        ind += 1
                                    end # End loop over survivng offspring
                                end # End if statement
                            end # End loop over index of 'survivng'
                            for q in 1:immigrants # Begin loop over immigrants
                                ind = q + sum(surviving) # Determines where in the array immigrants get put.
                                                         # q + sum(surviving) ensures they do not overwrite existing organisms.
                                immigrant = Array{Float32,1}(undef,length(landscape[i,j].species[p][1,1:end]))
                                ID = i              # Trait 1: Species ID number
                                T_opt = (rand()*grad*1.5)-((grad*1.5)/2)+trend #8+22 # Trait 2: Temperature optimum
                                T_sd = rand(LogNormal(0,1))    # Trait 3: Temperature tolerance
                                H_opt = (rand()*grad*1.5)-((grad*1.5)/2)      # Trait 4: Habitat optimum
                                H_sd = rand(LogNormal(0,1)) # Trait 5: Habitat tolerance
                                Disp_l = rand()      # Trait 6: Dispersal probability
                                Disp_g = rand()       # Trait 7: Global dispersal probability
                                Fert_max = 15       # Trait 8: Maximum number of offspring
                                dispersed = false   # Trait 9: Whether or not individual has already dispersed
                                lineage = rand(Float32) # Trait 10: Lineage identifier
                                immigrant = [ID, T_opt, T_sd, H_opt, H_sd, Disp_l, Disp_g, Fert_max, dispersed, lineage]
                                newgen[ind,1:end] = copy(immigrant)
                            end # End loop over immigrants
                        else
                            newgen = Array{Float32,2}(undef,sum(surviving),length(landscape[i,j].species[p][1,1:end]))
                            ind = 1 # Keeps count of individual offspring added to newgen array
                            for q in 1:length(surviving) # Goes down index of 'surviving'
                                if surviving[q] > 0
                                    for r in 1:surviving[q] # Loop from 1 to number of surviving offspring
                                        newgen[ind,1:end] = copy(landscape[i,j].species[p][q,1:end])  #
                                        newgen[ind,9] = false
                                        #println("$(newgen[ind,9])")
                                        ind += 1
                                    end # End loop over surviving offspring
                                end # End if statement
                            end # End loop over index of 'survivng'
                        end # End if-else statement
                        global landscape[i,j].species[p] = copy(newgen)
                    else
                        array = Array{Float32,2}(undef,0,length(landscape[i,j].species[p][1,1:end])) # If total offspring is 0, replaces landscape[i,j].species[p]
                        global landscape[i,j].species[p] = copy(array)                               # with a 0 by 8 array.
                        #println("Set length of landscape[$i,$j].species[$p] to zero")
                    end # End if-else statement
                else
                    #println("No individuals of species $p present")
                    array = Array{Float32,2}(undef,0,length(species_list[1,1:end]))
                    global landscape[i,j].species[p] = copy(array)                               # with a 0 by 8 array.
                    #println("Set length of landscape[$i,$j].species[$p] to zero")
                end # End if-else statement
            end # End second loop over species
        end # End landscape width loop
    end # End landscape length loop
end # End function

# FUnction for deterimining target patch for nearest neighbor dispersal. Patch is row no. ± 1 and column no. ± 1.
# Function will try to get a new target if the current target patch is the same as the natal patch.
# The function will also check if the target is outside the landscape.
function get_target_local(nrows::Int,ncols::Int,i::Int,j::Int)
    #println("get_target_local")
    #println("Getting local dispersal target patch...")
    target_row = i + rand(-1:1)
    target_col = j + rand(-1:1)
    #println("Target row = $target_row")
    #println("Target col = $target_col")
    outofbounds = false
    while target_row == i && target_col == j
        #println("Got same target as patch of origin. Rerolling...")
        target_row = i + rand(-1:1)
        target_col = j + rand(-1:1)
        if target_row > length(landscape[1:end,1]) || target_row < 1 || target_col > length(landscape[1,1:end]) || target_col < 1
            outofbounds = true
            #println("Individual went out of bounds")
        end
    end
    if target_row > length(landscape[1:end,1]) || target_row < 1 || target_col > length(landscape[1,1:end]) || target_col < 1
        outofbounds = true
        #println("Individual went out of bounds")
    end
    return target_row, target_col, outofbounds
end

# Function for determining target patch for global dispersal. Patch is determined by by randomly selecting
# row and column numbers between 1 and the max extent of the x and y dimensions of the landscape. Function will try
# to select a new target if the current target is the same as the natal patch.
# Note: It is currently not possible for global dispersal to move an individual out of the landscape. May be subject to change.
function get_target_global(nrows::Int,ncols::Int,i::Int,j::Int)
    #println("get_target_global")
    #println("Getting global dispersal target patch...")
    target_row = rand(1:length(landscape[1:end,1]))
    target_col = rand(1:length(landscape[1,1:end]))
    #println("Target patch @ row:$target_row, col:$target_col")
    outofbounds = false
    while target_row == i && target_col == j
        #println("Got same target as patch of origin. Rerolling...")
        target_row = rand(1:length(landscape[1:end,1]))
        target_col = rand(1:length(landscape[1,1:end]))
        if target_row != i || target_col != j
            #println("New target patch @ row:$target_row, col:$target_col")
        end
    end
    if target_row > length(landscape[1:end,1]) || target_row < 1 || target_col > length(landscape[1,1:end]) || target_col < 1
        outofbounds = true
        #println("Individual went out of bounds")
    end
    return target_row, target_col, outofbounds
end

function disperse_local(trow_l::Int, tcol_l::Int, outofbounds_l::Bool, i::Int, j::Int, k::Int, l::Int)
    #println("disperse_local")
    #println("outofbounds = $outofbounds_l")
    if outofbounds_l == false
        global landscape[i,j].species[k][l,9] = true # Change dispersed flag to true
        #println("migrant dispersed = $(landscape[i,j].species[k][l,9])")
        migrant = reshape(landscape[i,j].species[k][l,1:end],1,length(landscape[i,j].species[k][l,1:end])) # Reshapes vector from a 3 by 1 array to a 1 by 3 array. Needed for vcat to work.
        #println("migrant = $migrant")
        global landscape[trow_l,tcol_l].species[k] = vcat(landscape[trow_l,tcol_l].species[k],migrant) # Adds the individual to the target patch.
        #newindex = length(landscape[trow_l,tcol_l].species[k])
        #println("individual at $trow_l $tcol_l $k = $(landscape[trow_l,tcol_l].species[k][newindex,1:end])")
        #println("Added individual of species $k to patch $trow_l, $tcol_l ")
        global landscape[i,j].species[k] = landscape[i,j].species[k][setdiff(1:end,l),:] # Removes the individual from the natal patch
        #println("Removed individual of species $k from patch $i, $j")
    else
        global landscape[i,j].species[k] = landscape[i,j].species[k][setdiff(1:end,l),:]
        #println("Removed individual of species $k from patch $i, $j")
    end
end

function disperse_global(trow_g::Int,tcol_g::Int,outofbounds_g::Bool,i::Int,j::Int,k::Int,l::Int)
    #println("Initiating global dispersal")
    #println("outofbounds = $outofbounds_g")
    if outofbounds_g == false
        global landscape[i,j].species[k][l,9] = true
        migrant = reshape(landscape[i,j].species[k][l,1:end],1,length(landscape[i,j].species[k][l,1:end])) # Reshapes vector from a 3 by 1 array to a 1 by 3 array. Needed for vcat to work.
        global landscape[trow_g,tcol_g].species[k] = vcat(landscape[trow_g,tcol_g].species[k],migrant) # Adds the individual to the target patch.
        #println("Added individual of species $k to patch $trow_g, $tcol_g")
        global landscape[i,j].species[k] = landscape[i,j].species[k][setdiff(1:end,l),:] # Removes the individual from the natal patch
        #println("Removed individual of species $k from patch $i, $j")
    else
        global landscape[i,j].species[k] = landscape[i,j].species[k][setdiff(1:end,l),:]
        #println("Removed individual of species $k from patch $i, $j")
    end
end

function dispersal!(landscape::Array{TPatch,2})
    #println("dispersal!")
    nrows = length(landscape[1:end,1])
    ncols = length(landscape[1,1:end])
    for i in 1:length(landscape[1:end,1]) # Begin landscape row loop
        #println("Row $i")
        for j in 1:length(landscape[1,1:end]) # Begin landscape column loops
            #println("Col $j")
            if length(landscape[i,j].species)>0 # Check number of species
                #println("species length check")
                for k in 1:length(landscape[i,j].species) # loop across species
                    #println("species loop")
                    #sp_pop = length(landscape[i,j].species[k][1:end,1])
                    #println("Species $k pop. = $sp_pop")
                    popsize = length(landscape[i,j].species[k][1:end,1])
                    #println("pop = $popsize")
                    if length(landscape[i,j].species[k][1:end,1])>0 # Check population size
                        #println("pop size check")
                        l = 0
                        while l < length(landscape[i,j].species[k][1:end,1])
                            #println("while loop")
                            l += 1
                            #println("$(landscape[i,j].species[k][l,9])")
                            #d_global = rand()
                            #println("Individual $l: d_local = $d_local, d_global = $d_global")
                            if landscape[i,j].species[k][l,9] == false
                                #println("Individual has not dispersed. Checking roll against disp")
                                d_local = rand()
                                #println("disp = $(landscape[i,j].species[k][l,6]), roll = $d_local")
                                #println("roll <= disp = $(d_local<=landscape[i,j].species[k][l,6])")
                                if d_local<=landscape[i,j].species[k][l,6] # Check vs. dispersal probability (trait 6)
                                    #println("Individual disperses, checking mode")
                                    d_global = rand()
                                    #println("Roll = $d_global")
                                    #println("Roll <= disp_g = $(d_global <= landscape[i,j].species[k][l,7])")
                                    if d_global <= landscape[i,j].species[k][l,7] # Check vs. global dispersal probability (trait 7)
                                        #println("Getting global dispersal target")
                                        trow_g, tcol_g, outofbounds_g = get_target_global(nrows,ncols,i,j)
                                        disperse_global(trow_g,tcol_g,outofbounds_g,i,j,k,l)
                                    else
                                        #println("Getting local dispersal target")
                                        trow_l, tcol_l, outofbounds_l = get_target_local(nrows,ncols,i,j)
                                        disperse_local(trow_l,tcol_l,outofbounds_l,i,j,k,l)
                                    end # End if-else statement
                                end # End if statement
                            end # End if statement
                        end # End while loop
                    end # End if statement
                end # End loop over species
            end # End if statement
        end # End landscape column loop
    end # End landscape row loop
end # End function


#-------------------------
# Data Output Functions
#-------------------------

function set_filename(argumentsdict::Dict, par::Dict,clim,grad)
    filename = string("S", clim)
    println(filename)
    if grad == 1
        filename = string(filename, "_G0.5")
        println(filename)
    elseif grad == 2
        filename = string(filename, "_G1")
        println(filename)
    elseif grad == 3
        filename = string(filename, "_G2")
        println(filename)
    end
    if argumentsdict["burninperiod"] == true
        filename = string(filename,"_B", par["bmax"])
        println(filename)
    end
    filename = string(filename, "_T", par["tmax"])
    println(filename)
    if argumentsdict["uniformt"] == true
        filename = string(filename,"_UT")
        println(filename)
    else
        if argumentsdict["autocorrelatedtemp"] == true
            filename = string(filename,"_AT",argumentsdict["autocor_t"])
            println(filename)
        else
            filename = string(filename,"_ClT",argumentsdict["autocor_t"])
            println(filename)
        end
    end
    if argumentsdict["uniformenv"] == true
        filename = string(filename,"_UE")
        println(filename)
    else
        if argumentsdict["autocorrelatedenv"] == true
            filename = string(filename,"_Ae",argumentsdict["autocor_e"])
            println(filename)
        else
            filename = string(filename,"_Cle",argumentsdict["autocor_e"])
            println(filename)
        end
    end

    tempsource = argumentsdict["tempsource"]
    envsource = argumentsdict["envsource"]

    n = 1
    outname = string(filename,"_",n)
    searchname = string(filename,"_",n,".txt")

    files = string.(par["dir"],readdir(par["dir"]))
    while in(true,occursin.(Regex(searchname),files))==true
        n = n+1
        searchname = string(filename,"_",(n),".txt")
        outname = string(filename,"_",(n))
    end
    println(outname)
    return outname,tempsource,envsource
end

# Counts the total population in the landscape
function popcount(landscape)
    x = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                x = x + length(landscape[i,j].species[l][1:end,1])
            end
        end
    end
    #println("landscape population: $x")
    return x
end

#expected = expected_fert.(landscape[i,j].species[p][1:end,2],
#           landscape[i,j].species[p][1:end,3],
#           landscape[i,j].species[p][1:end,4],
#           landscape[i,j].species[p][1:end,5],
#           landscape[i,j].species[p][1:end,8],
#           landscape[i,j].temp_t,
#           landscape[i,j].habitat, niche_tradeoff,
#           niche_tradeoff, trend)

function fitnessmeans(landscape,trend,meantrend,α)
    obs = 0
    sumft = 0
    sumfh = 0
    sumfit = 0
    sumtdiff = 0
    sumavgtdiff = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                obs=obs+length(landscape[i,j].species[l][1:end,1])
                if size(landscape[i,j].species[l])[1]>0
                    for p in 1:length(landscape[i,j].species[l][1:end,1])
                        ft,fh,fit = expected_fert2.(landscape[i,j].species[l][p,2],
                                       landscape[i,j].species[l][p,3],
                                       landscape[i,j].species[l][p,4],
                                       landscape[i,j].species[l][p,5],
                                       landscape[i,j].temp_t,
                                       landscape[i,j].habitat,
                                       α,α,trend)
                        tdiff = landscape[i,j].species[l][p,2] .- (landscape[i,j].temp_t .+ trend)
                        avg_tdiff = landscape[i,j].species[l][p,2] .- (landscape[i,j].temp_t .+ meantrend)
                        sumtdiff = sumtdiff + sum(tdiff)
                        sumavgtdiff = sumavgtdiff + sum(avg_tdiff)
                        sumft = sumft + ft
                        sumfh = sumfh + fh
                        sumfit = sumfit + fit
                    end
                end
            end
        end
    end
    meantdiff = sumtdiff/obs
    meanavgtdiff = sumavgtdiff/obs
    meanft = sumft/obs
    meanfh = sumfh/obs
    meanfit = sumfit/obs
    ssqft = 0
    ssqfh = 0
    ssqfit = 0
    ssqtdiff = 0
    ssqavgtdiff = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                if size(landscape[i,j].species[l])[1]>0
                    for p in 1:length(landscape[i,j].species[l][1:end,1])
                        ft,fh,fit = expected_fert2.(landscape[i,j].species[l][p,2],
                                       landscape[i,j].species[l][p,3],
                                       landscape[i,j].species[l][p,4],
                                       landscape[i,j].species[l][p,5],
                                       landscape[i,j].temp_t,
                                       landscape[i,j].habitat,
                                       α,α,trend)
                        tdiff = landscape[i,j].species[l][p,2] .- (landscape[i,j].temp_t .+ trend)
                        avg_tdiff = landscape[i,j].species[l][p,2] .- (landscape[i,j].temp_t .+ meantrend)
                        ssqtdiff = ssqtdiff + sum((tdiff - meantdiff)^2)
                        ssqavgtdiff = ssqavgtdiff + sum((avg_tdiff - meanavgtdiff)^2)
                        ssqft = ssqft + sum((ft - meanft)^2)
                        ssqfh = ssqfh + sum((fh - meanfh)^2)
                        ssqfit = ssqfit + sum((fit - meanfit)^2)
                    end
                end
            end
        end
    end
    vartdiff = ssqtdiff/obs
    varavgtdiff = ssqavgtdiff/obs
    varft = ssqft/obs
    varfh = ssqfh/obs
    varfit = ssqfit/obs

    return meanft,meanfh,meanfit,meantdiff,meanavgtdiff,varft,varfh,varfit,vartdiff,varavgtdiff
end

# Calculates landscape-wide arithmetic means and variance for traits
function traitmeans(landscape)
    obs = 0
    sumtopt = 0
    sumtsd = 0
    sumhopt = 0
    sumhsd = 0
    sumdisp = 0
    sumdispg = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                obs=obs+length(landscape[i,j].species[l][1:end,1])
                sumtopt = sumtopt + sum(landscape[i,j].species[l][1:end,2])
                sumtsd = sumtsd + sum(landscape[i,j].species[l][1:end,3])
                sumhopt = sumhopt + sum(landscape[i,j].species[l][1:end,4])
                sumhsd = sumhsd + sum(landscape[i,j].species[l][1:end,5])
                sumdisp = sumdisp + sum(landscape[i,j].species[l][1:end,6])
                sumdispg = sumdispg + sum(landscape[i,j].species[l][1:end,7])
            end
        end
    end
    meantopt = sumtopt/obs
    meantsd = sumtsd/obs
    meanhopt = sumhopt/obs
    meanhsd = sumhsd/obs
    meandisp = sumdisp/obs
    meandispg = sumdispg/obs
    ssqtopt = 0
    ssqtsd = 0
    ssqhopt = 0
    ssqhsd = 0
    ssqdisp = 0
    ssqdispg = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                ssqtopt = ssqtopt + sum((landscape[i,j].species[l][1:end,2] .- meantopt).^2)
                ssqtsd = ssqtsd + sum((landscape[i,j].species[l][1:end,3] .- meantsd).^2)
                ssqhopt = ssqhopt + sum((landscape[i,j].species[l][1:end,4] .- meanhopt).^2)
                ssqhsd = ssqhsd + sum((landscape[i,j].species[l][1:end,5] .- meanhsd).^2)
                ssqdisp = ssqdisp + sum((landscape[i,j].species[l][1:end,6] .- meandisp).^2)
                ssqdispg = ssqdispg + sum((landscape[i,j].species[l][1:end,7] .- meandispg).^2)
            end
        end
    end
    vartopt = ssqtopt/obs
    vartsd = ssqtsd/obs
    varhopt = ssqhopt/obs
    varhsd = ssqhsd/obs
    vardisp = ssqdisp/obs
    vardispg = ssqdispg/obs
    return meantopt,meantsd,meanhopt,meanhsd,meandisp,meandispg,vartopt,vartsd,varhopt,varhsd,vardisp,vardispg
end

# Calculates species richness for a single patch
function landscape_div(landscape::Array{TPatch,2})
    lineages = Array{String,1}(undef,0)
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                lineages = [lineages;landscape[i,j].species[l][1:end,10]]
            end
        end
    end

    rich = length(unique(lineages))
    if rich > 0
        div = DataFrame()
        div.lineages = lineages
        div.tally = ones(length(div.lineages))
        div = combine(groupby(div,[:lineages])) do div
            DataFrame(count=sum(div.tally))
        end
        simp = simpson(div.count)
        shan = shannon(div.count)
        return rich,simp,shan
    else
        rich = 0
        simp = "NA"
        shan = "NA"
        return rich,simp,shan
    end
end

# Counts population for each species in each patch and outputs them to an array with dimensions
# rows by cols by n_species.
function patch_spp_pops(landscape::Array{TPatch,2},n_species)
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    patch_pops = Array{Array{Int,2},1}(undef, n_species)
    for p in 1:n_species
        pops = Array{Int,2}(undef,rows,cols)
        for i in 1:rows
            for j in 1:cols
                pops[i,j] = length(landscape[i,j].species[p][1:end,1])
            end
        end
        patch_pops[p] = pops
    end
    return patch_pops
end

# Calculates trait value means for each patch in the landscape and outputs a 3d array with dimensions
# rows by cols by n_traits.
# Note: This function could be more elegantly written. Consider rewriting this function.
function trait_analysis(landscape::Array{TPatch,2})
    #println("Trait Analysis")
    rows = length(landscape[1:end,1])
    #println("rows = $rows")
    cols = length(landscape[1,1:end])
    #println("cols = $cols")
    n_species = length(species_list[1:end,1])
    n_traits = length(species_list[1,2:8]) # 7 traits of interest,
    global trait_means = Array{Array{Float32,2}}(undef,n_traits)
    for k in 1:n_traits
        #println("trait loop = $k")
        means = Array{Float32,1}(undef,n_species)
        weights = Array{Int,1}(undef,n_species)
        mean_traits_values = Array{Float32,2}(undef,rows,cols)
        for i in 1:rows
            #println("Row = $i")
            for j in 1:cols
                if n_species > 1
                    #println("Col = $j")
                    for p in 1:n_species
                        #println("Species = $p")
                        weights[p] = length(landscape[i,j].species[p][1:end,1])
                        #println("Weight = $(length(landscape[i,j].species[p][1:end,1]))")
                        if length(landscape[i,j].species[p][1:end,1]) > 0
                            means[p] = mean(landscape[i,j].species[p][1:end,k+1])
                            #println("Mean = $(means[p])")
                        else
                            means[p] = 0
                            #println("Mean = $(means[p])")
                        end
                    end # End species loop
                end
                #println("Weighted mean = $(weightedmean(means,weights))")
                if n_species > 1
                    mean_traits_values[i,j] = weightedmean(means,weights)
                else
                    mean_traits_values[i,j] = mean(landscape[i,j].species[1][1:end,k+1])
                end
            end # End col loop
        end # End row loop
        global trait_means[k] = mean_traits_values
    end # End traits loop
end

function diversity_analysis(landscape::Array{TPatch,2})
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    n_species = length(species_list[1:end,1])
    global pops = Array{Array{Int,2},1}(undef, n_species)
    global richness = Array{Int,2}(undef,rows,cols)
    #global mean_traits_values = Array{Array{Float32,2},1}(undef,6)
    for i in 1:rows
        for j in 1:cols
            richness[i,j] = patch_species_richness(landscape,i,j)
        end
    end
    global pops = patch_spp_pops(landscape,n_species)
end

function env_analysis(landscape::Array{TPatch,2},trend_t)
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    temperatures = Array{Float32,2}(undef,rows,cols)
    habitat = Array{Float32,2}(undef,rows,cols)
    for i in 1:rows
        for j in 1:cols
            temperatures[i,j] = landscape[i,j].temp_t + trend_t
            habitat[i,j] = landscape[i,j].habitat
        end
    end
    global temperatures = temperatures
    global habitat = habitat
end

# Calculates unweighted mean stress for each patch
function mean_stress(landscape::Array{TPatch,2}, trend)
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    n_species = length(species_list[1:end,1])
    stress_t = Array{Float32,2}(undef,rows,cols) # Temperature stress
    stress_h = Array{Float32,2}(undef,rows,cols) # Habitat stress
    stress_o = Array{Float32,2}(undef,rows,cols) # Overall stress
    for i in 1:rows
        for j in 1:cols
            for p in 1:n_species
                S_T, S_H = stress.(landscape[i,j].species[p][1:end,2], # Calculating environmental stress for a species in a patch
                           landscape[i,j].species[p][1:end,3],
                           landscape[i,j].species[p][1:end,4],
                           landscape[i,j].species[p][1:end,5],
                           landscape[i,j].temp_t,
                           landscape[i,j].habitat, trend)
                S_O = (S_T .* S_H) # Calculate overall stress
                stress_t[i,j] = mean(S_T)
                stress_h[i,j] = mean(S_H)
                stress_o[i,j] = mean(S_O)
            end
        end
    end
    return stress_t,stress_h,stress_o
end

function write_landscape_stats(landscape,directory,filename,replicate,timestep,s_clim,trend,mean_trend,grad,H_t,H_h,α,bmax)
    outputname = string(directory,filename,"trend.txt")
    col_names1 = ["Replicate" "Timestep" "Pop" "rich" "simp" "shan" "T_opt" "T_sd" "H_opt" "H_sd" "disp" "disp_g" "VT_opt" "VT_sd" "VH_opt" "VH_sd" "Vdisp" "Vdisp_g"]
    col_names2 = ["ft" "fh" "fit" "tdiff" "avgtdiff" "varft" "varfh" "varfit" "var_tdiff" "var_avgtdiff" "clim_scen" "trend" "mean_trend" "grad" "H_t" "H_h" "alpha"]
    col_names = hcat(col_names1,col_names2)
    T_opt, T_sd, H_opt, H_sd, disp, disp_g, VT_opt, VT_sd, VH_opt, VH_sd, Vdisp, Vdisp_g = traitmeans(landscape)
    ft, fh, fit,tdiff,avgtdiff, varft, varfh, varfit, vartdiff, varavgtdiff = fitnessmeans(landscape,trend,mean_trend,α)
    rich, simp, shan = landscape_div(landscape)
    pop = popcount(landscape)
    open(outputname,"a") do IO
        if timestep==-1 && replicate==1
            writedlm(IO,col_names)
        end
        # Change to hcat
        writedlm(IO, [replicate timestep pop rich simp shan T_opt T_sd H_opt H_sd disp disp_g VT_opt VT_sd VH_opt VH_sd Vdisp Vdisp_g ft fh fit tdiff avgtdiff varft varfh varfit vartdiff varavgtdiff s_clim trend mean_trend grad H_t H_h α])
    end
end

# Writes the full landscape data to a .csv file, including every individual in every patch.
# Caution: Very large data file.
function write_landscape_csv(landscape,directory,filename,replicate,timestep,s_clim,trend,mean_trend,grad,H_t,H_h,α)
    outputname = string(directory,filename,".txt")
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    n_species = length(landscape[1,1].species[1:end,1])
    col_names = ["ID" "Replicate" "Timestep" "H_t" "H_h" "alpha" "clim_scen" "gradient" "x" "y" "T_opt" "T_sd" "H_opt" "H_sd" "disp_l" "disp_g" "fert_max" "LineageID" "temp_t" "trend" "mean_trend" "precip_t" "habitat"]
    open(outputname, "a") do IO
        if timestep==-1 && replicate==1
            writedlm(IO, col_names)
        end
        for i in 1:rows
            for j in 1:cols
                for k in 1:n_species
                    for l in 1:length(landscape[i,j].species[k][1:end,1])
                        #println("lin ID = ",string(landscape[i,j].species[k][l,10]))
                        lin = parse(Int,(@sprintf("%.0f",landscape[i,j].species[k][l,10]*10^15)))
                        lineageID = string(lin, base=62)
                        temperature = landscape[i,j].temp_t + trend
                        # Change to hcat
                        writedlm(IO, [k replicate timestep H_t H_h α s_clim grad i j landscape[i,j].species[k][l,2] landscape[i,j].species[k][l,3] landscape[i,j].species[k][l,4] landscape[i,j].species[k][l,5] landscape[i,j].species[k][l,6] landscape[i,j].species[k][l,7] landscape[i,j].species[k][l,8] lineageID temperature trend mean_trend landscape[i,j].precip_t landscape[i,j].habitat])
                    end
                end
            end
        end
        close(IO)
    end
end

# Creates heatmap graphs with patch values for environmental variables, richness, trait means,
# and individual species populations in .png format.
function heatmaps(richness,temperatures,habitats,pops,trait_means,filename,directory)
    type = ".png"
    rich = string(directory,"Rich_",filename,type)
    println(rich)
    temps = string(directory,"Temp_",filename,type)
    println(temps)
    hab = string(directory,"Habitat_",filename,type)
    println(hab)
    heatmap(richness)
    png(rich)
    heatmap(temperatures)
    png(temps)
    heatmap(habitats)
    png(hab)
    for i in 1:length(species_list[1:end,1])
        spdist = string("Spp_$i","_")
        spfilename = string(directory,spdist,filename,type)
        println(spfilename)
        heatmap(pops[i])
        png(spfilename)
    end
    for i in 1:length(trait_means)
        traits = string("Trait_$i","_")
        trfilename = string(directory,traits,filename,type)
        println(trfilename)
        heatmap(trait_means[i])
        png(trfilename)
    end
end

#------------------------------
# File manipulations
#------------------------------

# Merges data output files from different replicates into a single file containing all replicates
# Args: directory - location of files to be merged
#       filename -
function merge_output_files(directory,filename)
    for i in 1:length(infiles2)
        file = DataFrame(CSV.File(infiles2[i]))
        file.Replicate .= i
        if i==1
            CSV.write(outfile,file,append=false)
        else
            CSV.write(outfile,file,append=true)
        end
    end
end

#------------------------------
# Testing/Diagnostic Functions
#------------------------------

# Command line functions for testing/diagnostic purposes

# Prints the trait values for a species in the console
function get_sp_traits(sp_index::Int)
    if sp_index < 1 || sp_index > length(species_list[1:end,1])
        println("Please provide a valid species index value.")
    else
        print("$(species_list[sp_index,1:end])")

    end
end

# Lists out traits for all species in the console
function get_sp_list()
    n_sp = length(species_list[1:end,1])
    println("N_species = $n_sp")
    println("")
    for i in 1:length(species_list[1:end,1])
        for j in 1:length(species_list[1,1:end])
            print("$(species_list[i,j]) ")
        end
    println("")
    end
end

# Lists out patch environmental values in the console
function get_patch_values(x::Int, y::Int)
    println("Patch $x,$y:")
    println("Temperature = $(landscape[x,y].temp_t)")
    println("Precip_T = $(landscape[x,y].precip_t)")
    println("Habitat = $(landscape[x,y].habitat)")
end

function checkarray(array)
    for i in 1:length(array)
        println("$i : $(array[i])")
    end
end



#----------------------------------------------
# Main Program
#----------------------------------------------

# Run simulation using landscape and parameters from input files. Requires inputs for temperature, environment,
# scen,grad... are loop indices for different scenarios (e.g. climate trend, gradient strength)
function simulation_run()
    println("Starting.")
    ArgDict = read_arguments() # Read in shell input
    println(ArgDict)
    parasource = ArgDict["parasource"]
    println(parasource)
    println(typeof(parasource))
    par = include(parasource)
    #include("PatchStatistics.jl")
    scen = ArgDict["scenario"]
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
        write_landscape_stats(landscape,dir,filename,rep,-1,scen,0,0,grad,ArgDict["autocor_t"],ArgDict["autocor_e"],α,bmax)
        write_landscape_csv(landscape,dir,filename,rep,-1,scen,0,0,grad,ArgDict["autocor_t"],ArgDict["autocor_e"],α)
        println("Starting burn-in period.")
        if ArgDict["burninperiod"] == true
            for b in 1:bmax
                burnin = true
                #println("Beginning burn-in timestep $b.")
                #println("Starting dispersal routine.")
                dispersal!(landscape)
                #println("Starting reproduction routine.")
                demographics(landscape,α,0,grad,K,burnin,immi,p_immi,e_immi)
                if mut==true
                    mutate(landscape,p_mut,mut_sd,mut_decay,b)
                end
                write_landscape_stats(landscape,dir,filename,rep,b,scen,0,0,grad,ArgDict["autocor_t"],ArgDict["autocor_e"],α,bmax)
                if mod(b,50)==true || b==bmax
                    write_landscape_csv(landscape,dir,filename,rep,b,scen,0,0,grad,ArgDict["autocor_t"],ArgDict["autocor_e"],α)
                end
                #println("End timestep $t.")
            end
        end
        #write_landscape_csv(landscape,dir,filename,rep,0,scen,0,grad,ArgDict["autocor_t"],ArgDict["autocor_e"],α)
        println("Starting main simulation")
        for t in 1:tmax
            if ArgDict["burninperiod"] == true
                step=t+bmax
            else
                step=t
            end        
            burnin = false
            #println("Beginning timestep $t.")
            #println("Starting dispersal routine.")
            dispersal!(landscape)
            #println("Starting reproduction routine.")
            demographics(landscape,α,trend[t],grad,K,burnin,immi,p_immi,e_immi)
            if mut==true
                mutate(landscape,p_mut,mut_sd,mut_decay,step)
            end
            write_landscape_stats(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,ArgDict["autocor_t"],ArgDict["autocor_e"],α,bmax)
            if mod(t,50)==0 && t>=9900 || t==tmax
                write_landscape_csv(landscape,dir,filename,rep,step,scen,trend[t],mean_trend[t],grad,ArgDict["autocor_t"],ArgDict["autocor_e"],α)
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
    patchstats(dir,filename)
    println("Program complete.")
end

#----------------------------------------------
# Main Script
#----------------------------------------------
# s_clim is the climate scenario, s_grad is the gradient scenario
#simulation_run2(parasource,0)

@time begin
    include("PatchStatistics.jl")
    simulation_run()
end
