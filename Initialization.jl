#-----------------------------------------
# Initialization and Data Input Functions
#-----------------------------------------

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

#--------------------------------------
# Simulation Initialization Functions
#--------------------------------------

# Creates an array with dimentions n_species and n_traits which acts as a master list of all species in the simulation
function init_spp(n_species::Int)
    #println("Initializing species master list...")
    n_traits = 13
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
        origin_patch_x = "x"
        origin_patch_y = "y"
        origin_time = "timestep"
        global species_list[i,1:end] = [ID, T_opt, T_sd, H_opt, H_sd, Disp_l, Disp_g, Fert_max, dispersed, origin_patch_x, origin_patch_y, origin_time]
    end
end

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
function init_popgrad(n_pop::Int,par::Dict,x::Int,y::Int)
    #println("init_popgrad")
    n_traits = 13
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
        origin_patch_x = x
        origin_patch_y = y
        origin_time = -1
        population[i,1:end] = [ID, T_opt, T_sd, H_opt, H_sd, Disp_l, Disp_g, Fert_max, dispersed, lineage, origin_patch_x ,origin_patch_y, origin_time,]
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
    init_pop = parameters["initial_pop"]
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
    for i in 1:nrows
        for j in 1:ncols
            row = i
            col = j
            #patchpop = init_popgrad(init_pop,parameters,i,j)
            if (i==1 ) && (j==1) # Use this to initialize a population in only one patch
                patchpop = init_popgrad(init_pop,parameters,i,j)
                println(typeof(patchpop))
                len = length(patchpop[1][1:end,1])
                wid = length(patchpop[1][1,1:end])
                println("Initial patchpop dims: $len, $wid")
            else
                patchpop = Array{Array{Float32,2},1}(undef,1)
                patchpop[1] = Array{Float32,2}(undef,0,10)
            end
            temp = temperatures[i,j]
            habitat = environments[i,j]
            precip = 0
            patch = TPatch(row,col,patchpop,temp,precip,habitat)
            println("patch: $i, $j, ", length(patch.species[1][1:end,1]))
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
