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
