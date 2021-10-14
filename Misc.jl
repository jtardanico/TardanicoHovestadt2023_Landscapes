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
