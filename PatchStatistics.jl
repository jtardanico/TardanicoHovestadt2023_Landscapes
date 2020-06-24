# Landscape statistics calculation

#----------------------------------------------------------------------
# This is a Julia implementation of the R script LandscapePatchStat2.R
#----------------------------------------------------------------------
using DataFrames
using CSV

# Calculation functions
function tfit(T_opt,T_sd,H_opt,H_sd, temp_t, habitat, α)
    fitness = exp(-(temp_t-t_opt)^2/(2*T_sd^2)) * exp(-t_sd^2/(2*α^2))
    return fitness
end

function hfit(T_opt,T_sd,H_opt,H_sd, temp_t, habitat, α)
    fitness = exp(-(habitat-h_opt)^2/(2*H_sd^2)) * exp(-h_sd^2/(2*α^2))
    return fitness
end

function fitness(T_opt,T_sd,H_opt,H_sd, temp_t, habitat, α)
    fitness =  exp(-(temp_t-t_opt)^2/(2*T_sd^2)) * exp(-(habitat-h_opt)^2/(2*H_sd^2)) * exp(-t_sd^2/(2*α^2)) * exp(-h_sd^2/(2*α^2))
    return fitness
end

data = DataFrame(CSV.File("/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/test5/S0_G0.5_B50_T50_UT_Ae0.0.txt"))

data.H_t = string.(data.H_t)
data.H_h = string.(data.H_h)

data.H_t = replace!(data.H_t,"NaN" => "Uniform")
data.H_h = replace!(data.H_h,"NaN" => "Uniform")

# Patch population sizes

populations = by(data, [:Timestep,:x,:y]) do data
    DataFrame(pop=length(data.ID))
end

rich = by(data, [:Timestep,:x,:y]) do data
    DataFrame(richness=length(unique(data.LineageID)))
end

ft = by(data, [:Timestep,:x,:y]) do data
    DataFrame(ft=mean(data.))


#div = select(data,[:Timestep,:x,:y,:LineageID])
#div.tally = 1
#div = by(div,[:Timestep,:x,:y,:LineageID]) do div
#    DataFrame(count=sum(div.tally))
#end
#rich = by(div,[:Timestep,:x,:y]) do div
#    DataFrame(rich=sum())
