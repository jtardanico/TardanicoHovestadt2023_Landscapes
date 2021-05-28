#-------------------------
# Simulation Functions
#-------------------------

#-----------------------------
# Demograpic processes

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
function demographics(landscape::Array{TPatch, 2},niche_tradeoff, trend, grad, K::Int,burnin, immi,p_immi,e_immi,timestep)
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
                                origin_patch_x = i
                                origin_patch_y = j
                                origin_time = timestep
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

#----------------------------------------
# Dispersal processes

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
