
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

#--------------------------
# Input Files
#--------------------------

parasource = "/home/joseph/github/LandKlifModel/ParameterDict1.jl"

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

function expected_fert(T_opt,T_sd,H_opt,H_sd,Fert_max, temp_t, habitat, α_t, α_h, t_ref, trend)
    temp = temp_t + t_ref + trend # Patch temperature
    e_fert = Fert_max * exp(-(temp-T_opt)^2/(2*T_sd^2)) * exp(-(habitat-H_opt)^2/(2*H_sd^2)) * exp(-T_sd^2/(2*α_t^2)) * exp(-H_sd^2/(2*α_h^2))
    return e_fert
end

# Function for calculating expected offspring mortality
# Terms: max fertility (Fert_max), number of offspring (offspring), carrying capacity (k).
function expected_mort(Fert_max, offspring, k)
    a = (Fert_max - 1)/(k * Fert_max)
    e_mort = 1/(1+a*offspring)
    return e_mort
end

# Function for calculating how stressful the environment of a patch is
function stress(T_opt,T_sd,H_opt,H_sd,temp_t, habitat,t_ref,trend)
    temp = temp_t + t_ref + trend
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

function read_arguments()
    s=ArgParseSettings()
    @add_arg_table s begin
        "--tempsource", "-t" # Name of temperature source file in the shell script
            help = "Landscape temperature source file"
        "--envsource", "-e"
            help = "Landscape environment source file"
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
        T_opt = rand()*20+10 #8+22 # Trait 2: Temperature optimum
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
function init_popgrad(n_pop::Int)
    n_traits = 10
    patchpop = Array{Array{Float32,2},1}(undef,1)
    population = Array{Float32,2}(undef,n_pop,n_traits)
    for i in 1:n_pop
        ID = i              # Trait 1: Species ID number
        T_opt = rand()*20+10 #8+22 # Trait 2: Temperature optimum
        T_sd = rand()*10    # Trait 3: Temperature tolerance
        H_opt = rand()      # Trait 4: Habitat optimum
        H_sd = rand() # Trait 5: Habitat tolerance
        Disp_l = rand()      # Trait 6: Dispersal probability
        Disp_g = rand()       # Trait 7: Global dispersal probability
        Fert_max = 15       # Trait 8: Maximum number of offspring
        dispersed = false   # Trait 9: Whether or not individual has already dispersed
        lineage = abs(rand(Int)) # Trait 10: Lineage identifier
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
    if trend_type == 1
        global trend = rand(Normal(mean,sd),nsteps) # Random variation, constant mean
    elseif trend_type == 2
        t = Array{Float32,1}(undef,nsteps)
        for i in 1:nsteps
            t[i] = 0.025*i
        end
        global trend = t
    elseif trend_type == 3
        t = Array{Float32,1}(undef,nsteps)
        for i in 1:nsteps
            t[i] = 0.025*i + rand(Normal(mean,sd))
        end
        global trend = t # Random variation & rising mean
    elseif trend_type == 0
        global trend = zeros(nsteps)
    else
        println("Invalid scenario. Valid scenarios are:")
        println("0: Constant mean, no random variation")
        println("1: Constant mean temperature, constant variance")
        println("2: Increasing mean temperature, no random variation")
        println("3: Increasing mean, constant variance")
        exit(code=0)
    end
end

function gradient_scenario(scenario)
    if scenario == 1
        gradient = 0.5
    elseif scenario == 2
        gradient = 1
    elseif scenario == 3
        gradient = 2
    else
        error("Invalid scenario. No such scenario implemented.")
    end
    return gradient
end
# Initializes a landscape from file inputs.
# 'gradient' defines landscape gradient steepness
function init_world(worldtempsource::String,worldenvsource::String,gradient)
    temperatures, environments = read_world_source(worldtempsource,worldenvsource)
    temperatures = temperatures * gradient
    environments = environments * gradient
    nrows = length(temperatures[1:end,1])
    ncols = length(temperatures[1,1:end])
    global landscape = Array{TPatch, 2}(undef, nrows, ncols)
    for i in 1:nrows
        for j in 1:ncols
            row = i
            col = j
            patchpop = init_popgrad(100)
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

# Calculates number of offspring for each individual. Offspring form the next generation of individuals.
function demographics(landscape::Array{TPatch, 2},niche_tradeoff, T_ref, trend, K::Int)
    #println("Starting demographics routine...")
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            offspring = Array{Array{Int,1},1}(undef,length(landscape[i,j].species[1:end]))
            total_offspring = 0
            for p in 1:length(landscape[i,j].species[1:end])
                sp_pop = length(landscape[i,j].species[p][1:end,1])
                #println("Species $p pop. = $sp_pop")
                if length(landscape[i,j].species[p][1:end,1])>0
                    expected = expected_fert.(landscape[i,j].species[p][1:end,2],
                               landscape[i,j].species[p][1:end,3],
                               landscape[i,j].species[p][1:end,4],
                               landscape[i,j].species[p][1:end,5],
                               landscape[i,j].species[p][1:end,8],
                               landscape[i,j].temp_t,
                               landscape[i,j].habitat, niche_tradeoff,
                               niche_tradeoff, T_ref, trend)
                    #println("expected offspring = $expected")
                    global x = expected
                    offspring[p] = rand.(Poisson.(Float64.(expected)))
                    #println("Offspring = $offspring")
                    species_offspring = sum(offspring[p]) # For diagnostic purposes. Will be removed once testing is finished
                    #println("Total offspring of species $p = $species_offspring")
                    total_offspring = total_offspring + species_offspring
                end
            end
            for p in 1:length(landscape[i,j].species[1:end])
                if length(landscape[i,j].species[p][1:end,1])>0
                    S_T, S_H = stress(landscape[i,j].species[p][1,2], # Calculating environmental stress for a species in a patch
                               landscape[i,j].species[p][1,3],
                               landscape[i,j].species[p][1,4],
                               landscape[i,j].species[p][1,5],
                               landscape[i,j].temp_t,
                               landscape[i,j].habitat,
                               T_ref, trend)
                    #println("S_T = $S_T, S_H = $S_H")
                    #K = carry_capacity(300,S_T,S_H) # Calculates carrying capacity of a patch
                    #println("K = $K")
                    p_surv = expected_mort.(landscape[i,j].species[p][1:end,8],total_offspring,K)
                    if K > 0 # Note: May be unecessary in future model iterations. Consider removal.
                        surviving = rand.(Binomial.(offspring[p], Float64.(p_surv)))
                    else
                        surviving = 0
                    end
                    #println("Species $p has $surviving surviving offspring.")
                    if sum(surviving) > 0
                        newgen = Array{Float32,2}(undef,sum(surviving),length(landscape[i,j].species[p][1,1:end])) # Creates an array of length sum(offspring) with data for species p
                        ind = 1 # Keeps count of individual offspring added to newgen array
                        for q in 1:length(surviving) # Goes down index of 'surviving'
                            if surviving[q] > 0
                                for r in 1:surviving[q] # Loop from 1 to number of surviving offspring
                                    newgen[ind,1:end] = landscape[i,j].species[p][q,1:end]  #
                                    ind += 1
                                end
                            end
                        end
                        global landscape[i,j].species[p] = copy(newgen)
                    else
                        array = Array{Float32,2}(undef,0,length(landscape[i,j].species[p][1,1:end])) # If total offspring is 0, replaces landscape[i,j].species[p]
                        global landscape[i,j].species[p] = array                               # with a 0 by 8 array.
                        #println("Set length of landscape[$i,$j].species[$p] to zero")
                    end
                else
                    #println("No individuals of species $p present")
                    array = Array{Float32,2}(undef,0,length(species_list[1,1:end]))
                    global landscape[i,j].species[p] = array                               # with a 0 by 8 array.
                    #println("Set length of landscape[$i,$j].species[$p] to zero")
                end
            end
        end
    end
end

# FUnction for deterimining target patch for nearest neighbor dispersal. Patch is row no. ± 1 and column no. ± 1.
# Function will try to get a new target if the current target patch is the same as the natal patch.
# The function will also check if the target is outside the landscape.
function get_target_local(nrows::Int,ncols::Int,i::Int,j::Int)
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
    #println("Initiating local dispersal")
    #println("outofbounds = $outofbounds_l")
    if outofbounds_l == false
        global landscape[i,j].species[k][l,9] = true # Change dispersed flag to true
        migrant = reshape(landscape[i,j].species[k][l,1:end],1,length(landscape[i,j].species[k][l,1:end])) # Reshapes vector from a 3 by 1 array to a 1 by 3 array. Needed for vcat to work.
        global landscape[trow_l,tcol_l].species[k] = vcat(landscape[trow_l,tcol_l].species[k],migrant) # Adds the individual to the target patch.
        #println("Added individual of species $k to patch $trow_l, $tcol_l")
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
    nrows = length(landscape[1:end,1])
    ncols = length(landscape[1,1:end])
    for i in 1:length(landscape[1:end,1])
        #println("Row $i")
        for j in 1:length(landscape[1,1:end]) # Row & column loops
            #println("Col $j")
            if length(landscape[i,j].species)>0 # Check number of species
                for k in 1:length(landscape[i,j].species) # loop across species
                    sp_pop = length(landscape[i,j].species[k][1:end,1])
                    #println("Species $k pop. = $sp_pop")
                    if length(landscape[i,j].species[k][1:end,1])>0 # Check population size
                        l = 0
                        while l < length(landscape[i,j].species[k][1:end,1])
                            l += 1
                            d_local = rand()
                            #d_global = rand()
                            #println("Individual $l: d_local = $d_local, d_global = $d_global")
                            if landscape[i,j].species[k][l,9] == false
                                if d_local<=landscape[i,j].species[k][l,6] # Check vs. dispersal probability (trait 6)
                                    d_global = rand()
                                    if d_global <= landscape[i,j].species[k][l,7] # Check vs. global dispersal probability (trait 7)
                                        trow_g, tcol_g, outofbounds_g = get_target_global(nrows,ncols,i,j)
                                        disperse_global(trow_g,tcol_g,outofbounds_g,i,j,k,l)
                                    else
                                        trow_l, tcol_l, outofbounds_l = get_target_local(nrows,ncols,i,j)
                                        disperse_local(trow_l,tcol_l,outofbounds_l,i,j,k,l)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


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
    return filename,tempsource,envsource
end

# Calculates species richness for a single patch
function patch_species_richness(landscape::Array{TPatch,2},row::Int,col::Int)
    rich = 0
    for p in 1:length(species_list[1:end,1])
        if length(landscape[row,col].species[p][1:end,1]) > 0
            rich += 1
        end
    end
    return rich
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

# Calculates unweighted mean stress
function mean_stress(landscape::Array{TPatch,2},T_ref, trend)
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
                           landscape[i,j].habitat,
                           T_ref, trend)
                S_O = (S_T .* S_H) # Calculate overall stress
                stress_t[i,j] = mean(S_T)
                stress_h[i,j] = mean(S_H)
                stress_o[i,j] = mean(S_O)
            end
        end
    end
    return stress_t,stress_h,stress_o
end

# Currently broken
function write_sp_list(species_list,directory)
    rows = length(species_list[1:end,1])
    cols = length(species_list[1,1:(end-1)])
    col_names = ["ID" "T_opt" "T_sd" "H_opt" "H_sd" "disp_l" "disp_g" "fert_max"]
    outputname = string(directory,"Sp_list",".txt")
    open(outputname, "a") do IO
        writedlm(IO, col_names)
        for i in 1:rows
            for j in 1:cols
                writedlm(IO, species_list[i,j])
            end
        end
    end
end

# Writes the full landscape data to a .csv file, including every individual in every patch.
# Caution: Very large data file.
function write_landscape_csv(landscape,directory,filename,timestep,s_clim,grad,H_t,H_h)
    outputname = string(directory,filename,".txt")
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    n_species = length(landscape[1,1].species[1:end,1])
    col_names = ["ID" "Timestep" "H_t" "H_h" "clim_scen" "gradient" "x" "y" "T_opt" "T_sd" "H_opt" "H_sd" "disp_l" "disp_g" "fert_max" "LineageID" "temp_t" "precip_t" "habitat"]
    open(outputname, "a") do IO
        writedlm(IO, col_names)
        for i in 1:rows
            for j in 1:cols
                for k in 1:n_species
                    for l in 1:length(landscape[i,j].species[k][1:end,1])
                        lin = Int(landscape[i,j].species[k][l,10])
                        lineageID = string(lin, base=16)
                        writedlm(IO, [k timestep H_t H_h s_clim grad i j landscape[i,j].species[k][l,2] landscape[i,j].species[k][l,3] landscape[i,j].species[k][l,4] landscape[i,j].species[k][l,5] landscape[i,j].species[k][l,6] landscape[i,j].species[k][l,7] landscape[i,j].species[k][l,8] lineageID landscape[i,j].temp_t landscape[i,j].precip_t landscape[i,j].habitat])
                    end
                end
            end
        end
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

#----------------------------------------------
# Main Program
#----------------------------------------------

# Used for testing the program [!! UNDER CONSTRUCTION !!]
# Other functions needed:
# - Shell independent filename function.
function simulation_test(parasource::String, s1::Int, grad::Int)
    par = include(parasource)
    bmax = par["bmax"]
    tmax = par["tmax"]
    α = par["α3"]
    T_ref = par["T_ref"]
    dir = par["dir"]
    filename,tempsource,envsource = set_filename(argumentsdict,par,s1,grad)
    Random.seed!(123)
    println("Initializing world.")
    grad = gradient_scenario(grad)
    init_world(tempsource,s1)
    generate_climate_trend(tmax,0,1,s1)
    init_spp2()
    println("Starting burn-in period.")
    for b in 1:bmax
        println("Beginning timestep $b.")
        #println("Starting dispersal routine.")
        dispersal!(landscape)
        #println("Starting reproduction routine.")
        demographics(landscape,α,T_ref,trend[1],300)
        #println("End timestep $t.")
    end
    println("Starting ")
    for t in 1:tmax
        println("Beginning timestep $t.")
        #println("Starting dispersal routine.")
        dispersal!(landscape)
        #println("Starting reproduction routine.")
        demographics(landscape,α,T_ref,trend[t],300)
        #println("End timestep $t.")
    end
    println("End time loop. Analyzing landscape.")
    diversity_analysis(landscape)
    trait_analysis(landscape)
    env_analysis(landscape,trend[tmax])
    heatmaps(richness,temperatures,habitat,pops,trait_means,filename,dir)
    write_landscape_csv(landscape,dir,filename)
    #write_sp_list(species_list,dir)
    println("Program complete.")
end

# Run simulation with automatically generated landscape. Loops over 4 climate trends.
# 10 species populations.
function simulation_run1(parasource::String)
    bmax = par["bmax"]
    tmax = par["tmax"]
    α = par["α2"]
    T_ref = par["T_ref"]
    dir = par["dir"]
    start_pop_size = par["initial_pop"]
    for s in 0:3
        Random.seed!(123)
        println("Initializing species.")
        init_spp(start_pop_size)
        println("Initializing world.")
        generate_world(50,50)
        generate_climate_trend(tmax,0,1,s)
        println("Starting burn-in period.")
        for b in 1:bmax
            println("Beginning timestep $b.")
            #println("Starting dispersal routine.")
            dispersal!(landscape)
            #println("Starting reproduction routine.")
            demographics(landscape,α,T_ref,trend[1],300)
            #println("End timestep $t.")
        end
        println("Starting time loop.")
        for t in 1:tmax
            println("Beginning timestep $t.")
            #println("Starting dispersal routine.")
            dispersal!(landscape)
            #println("Starting reproduction routine.")
            demographics(landscape,α,T_ref,trend[t],300)
            #println("End timestep $t.")
        end
        println("End time loop. Analyzing landscape.")
        diversity_analysis(landscape)
        trait_analysis(landscape)
        env_analysis(landscape,trend[tmax])
        heatmaps(richness,temperatures,habitat,pops,trait_means,tmax,s,dir)
    end
    #write_landscape_csv(landscape,5,"Output")
    #write_sp_list(species_list,dir)
    println("Program complete.")
end

# Run simulation using landscape and parameters from input files. Requires inputs for temperature, environment,
# s1,grad... are loop indices for different scenarios (e.g. climate trend, gradient strength)
function simulation_run2(parasource::String, s1::Int)
    println("Starting.")
    println("Climate scenario $s1")
    println("parasource = $parasource")
    par = include(parasource)
    ArgDict = read_arguments()
    bmax = par["bmax"]
    tmax = par["tmax"]
    grad_str = par["grad"]
    α = par["α1.5"]
    T_ref = par["T_ref"]
    dir = par["dir"]
    println(dir)
    filename,tempsource,envsource = set_filename(ArgDict,par,s1,grad_str)
    Random.seed!(123)
    println("Initializing world.")
    grad = gradient_scenario(grad_str)
    init_world(tempsource,envsource,grad)
    generate_climate_trend(tmax,0,1,s1)
    init_spp2()
    println("Starting burn-in period.")
    if ArgDict["burninperiod"] == true
        for b in 1:bmax
            println("Beginning timestep $b.")
            #println("Starting dispersal routine.")
            dispersal!(landscape)
            #println("Starting reproduction routine.")
            demographics(landscape,α,T_ref,trend[1],300)
            #println("End timestep $t.")
        end
    end
    println("Starting ")
    for t in 1:tmax
        println("Beginning timestep $t.")
        #println("Starting dispersal routine.")
        dispersal!(landscape)
        #println("Starting reproduction routine.")
        demographics(landscape,α,T_ref,trend[t],300)
        #println("End timestep $t.")
    end
    println("End time loop. Analyzing landscape.")
    diversity_analysis(landscape)
    trait_analysis(landscape)
    env_analysis(landscape,trend[tmax])
    heatmaps(richness,temperatures,habitat,pops,trait_means,filename,dir)
    write_landscape_csv(landscape,dir,filename,tmax,s1,grad,ArgDict["autocor_t"],ArgDict["autocor_e"])
    #write_sp_list(species_list,dir)
    println("Program complete.")
end

#----------------------------------------------
# Main Script
#----------------------------------------------
# s_clim is the climate scenario, s_grad is the gradient scenario
@time for s_clim in 0:3
            simulation_run2(parasource,s_clim)
        end
