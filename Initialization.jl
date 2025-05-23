
#-----------------------------------------
# Initialization and Data Input Functions
#-----------------------------------------

#-------------------------
# Data Input Functions
#-------------------------

# Reads landscape T and H data from files.
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
            help = "Landscape T source file"
        "--envsource", "-e"
            help = "Landscape H source file"
        "--scenario", "-s"
            help = "Climate scenario"
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
        "--seed", "-r"
            help = "A seed for the RNG. If specified in the shell, the program will always use this seed. Seed may alternatively be specified in parameters."
            arg_type = Int64
        end
        return parse_args(s)
end

#--------------------------------------
# Simulation Initialization Functions
#--------------------------------------

# Creates a 1x9 array listing species traits. For use with function "simulation_run2"
function init_spp2()
    #println("Initializing species master list...")
    n_traits = 13
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
        origin_patch_x = "x"
        origin_patch_y = "y"
        origin_time = "timestep"
        global species_list[1,1:end] = [ID, T_opt, T_sd, H_opt, H_sd, Disp_l, Disp_g, Fert_max, dispersed, lineage, origin_patch_x, origin_patch_y, origin_time]
end

# Initialized a population of individuals with randomized trait values. Function scales distribution of niche optima to the parameter "grad"
#
function init_popgrad(n_pop::Int,par::Dict,x::Int,y::Int)
    #println("init_popgrad")
    n_traits = 13
    patchpop = Array{Array{Float32,2},1}(undef,1)
    population = Array{Float32,2}(undef,n_pop,n_traits)
    for i in 1:n_pop
        ID = i              # Trait 1: Species ID number
        T_opt = rand(Normal(par["topt_µ"],par["topt_σ"])) * par["grad"] #8+22 # Trait 2: Temperature optimum
        T_sd = rand(LogNormal(par["tsd_µ"],par["tsd_σ"]))    # Trait 3: Temperature tolerance
        H_opt = rand(Normal(par["hopt_µ"],par["hopt_σ"])) * par["grad"]   # Trait 4: Habitat optimum
        H_sd = rand(LogNormal(par["hsd_µ"],par["hsd_σ"])) # Trait 5: Habitat tolerance
        Disp_l = rand()      # Trait 6: Dispersal probability
        Disp_g = rand()       # Trait 7: Global dispersal probability
        Fert_max = 15       # Trait 8: Maximum number of offspring
        dispersed = false   # Trait 9: Whether or not individual has already dispersed
        lineage = rand(Float32) # Trait 10: Lineage identifier
        origin_patch_x = x # Trait 11: X index of patch of first appearance
        origin_patch_y = y # Trait 12: Y index of patch of first appearance
        origin_time = -1 # Trait 13: Timestep of first appearance
        population[i,1:end] = [ID, T_opt, T_sd, H_opt, H_sd, Disp_l, Disp_g, Fert_max, dispersed, lineage, origin_patch_x ,origin_patch_y, origin_time]
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
        trend = rand(Normal(mean,sd),nsteps) # Random variation, constant mean
        mean_trend = zeros(nsteps) # mean trend
    elseif trend_type == 2
        t = Array{Float32,1}(undef,nsteps)
        for i in 1:nsteps
            t[i] = 0.025*i
        end
        trend = t
        mean_trend = t # mean trend
    elseif trend_type == 3
        t = Array{Float32,1}(undef,nsteps)
        for i in 1:nsteps
            t[i] = 0.025*i + rand(Normal(mean,sd))
            t_mean =  0.025*i # mean trend
        end
        trend = t # Random variation & rising mean
        mean_trend = t_mean
    elseif trend_type == 0
        trend = zeros(nsteps)
        mean_trend = zeros(nsteps)
    else
        println("Invalid scenario. Valid scenarios are:")
        println("0: Constant mean, no random variation")
        println("1: Constant mean temperature, constant variance")
        println("2: Increasing mean temperature, no random variation")
        println("3: Increasing mean, constant variance")
        exit(code=0)
    end
    return trend, mean_trend
end

# Generates a second climate trend
function generate_climate_trend2(nsteps::Int,mean,sd,change)
    t = Array{Float32,1}(undef,nsteps)
    t_mean = Array{Float32,1}(undef,nsteps)
    for i in 1:nsteps
        t[i] = change*i + rand(Normal(mean,sd))
        t_mean[i] =  change*i # mean trend
    end
    trend2 = t
    mean_trend2 = t_mean
    return trend2, mean_trend2
end

# Initializes a landscape from file inputs.
# 'gradient' defines landscape gradient steepness
function init_world(worldtempsource::String,worldenvsource::String,gradient,parameters::Dict)
    #println("init_world")
    if haskey(parameters,"initial_pop")==true
        init_pop = parameters["initial_pop"]
    end
    temperatures, environments = read_world_source(worldtempsource,worldenvsource)
    typeof(temperatures) # For debugging diagnostics
    typeof(environments)
    temperatures = temperatures * gradient
    environments = environments * gradient
    nrows = length(temperatures[1:end,1])
    println("nrows = $nrows")
    ncols = length(temperatures[1,1:end])
    println("ncols = $ncols")
    m_t, v_t = mean_and_var(temperatures)
    m_e, v_e = mean_and_var(environments)
    global landscape = Array{TPatch, 2}(undef, nrows, ncols)
    if haskey(parameters,"initial_pop")==true
        for i in 1:nrows
            for j in 1:ncols
                row = i
                col = j
                patchpop = init_popgrad(init_pop,parameters,i,j)

                #if (i==1 ) && (j==1) # Use this to initialize a population in only one patch
                #    patchpop = init_popgrad(init_pop,parameters,i,j)
                #    println(typeof(patchpop))
                #    len = length(patchpop[1][1:end,1])
                #    wid = length(patchpop[1][1,1:end])
                #    println("Initial patchpop dims: $len, $wid")
                #else
                #    patchpop = Array{Array{Float32,2},1}(undef,1)
                #    patchpop[1] = Array{Float32,2}(undef,0,13)
                #end
                temp = temperatures[i,j]
                habitat = environments[i,j]
                precip = 0
                patch = TPatch(row,col,patchpop,temp,precip,habitat)
                #println("patch: $i, $j, ", length(patch.species[1][1:end,1]))
                landscape[i,j] = patch
            end
        end
    else
        for i in 1:nrows
            for j in 1:ncols
                row = i
                col = j
                patchpop = Array{Array{Float32,2},1}(undef,1)
                patchpop[1] = Array{Float32,2}(undef,0,13)
                temp = temperatures[i,j]
                habitat = environments[i,j]
                precip = 0
                patch = TPatch(row,col,patchpop,temp,precip,habitat)
                #println("patch: $i, $j, ", length(patch.species[1][1:end,1]))
                landscape[i,j] = patch
            end
        end
    end
end