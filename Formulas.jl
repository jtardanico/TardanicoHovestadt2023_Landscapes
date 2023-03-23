
# This module contains mathematical functions for use with the LandscapesModel simulation program.

module LandKlifFormulas

export expected_fert, expected_mort, stress, carry_capacity, gaussian_landscape, weightedmean

# Function for calculating the expected fertility (number of offspring) of an organism given its current fitness in its patch
function expected_fert(T_opt,T_sd,H_opt,H_sd,Fert_max, temp_t, habitat, α_t, α_h, t_ref, trend)
    temp = temp_t + t_ref + trend # Patch temperature
    e_fert = Fert_max * exp(-(temp-T_opt)^2/(2*T_sd^2)) * exp(-(habitat-H_opt)^2/(2*H_sd^2)) * exp(-T_sd^2/2*α_t^2) * exp(-H_sd^2/2*α_h^2)
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

end
