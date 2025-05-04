#-------------------------------
# Misc. and Deprecated Functions
#-------------------------------

#------------------------------
# Testing/Diagnostic Functions


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

function get_patchpop_dims(landscape)
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for k in 1:length(landscape[i,j].species)
                display(landscape[i,j].species[k])
                println(" ")
            end
        end
    end
end

#------------------------------------------------------------------
# Unused, deprecated or defuct functions


# CALCULATION FUNCTIONS


# DEFUNCT: Determines carry capacity for a species; may be removed/replaced later.
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


# INITIALIZATION FUNCTIONS


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

# Creates an array with dimentions n_species and n_traits which acts as a master list of all species in the simulation
# This function can be removed
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