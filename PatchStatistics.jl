# Landscape statistics calculation

#----------------------------------------------------------------------
# This is a Julia implementation of the R script LandscapePatchStat2.R
#----------------------------------------------------------------------
using DataFrames
using CSV
using Plots

# Calculation functions
function tfit(t_opt,t_sd,temp_t, α)
    fitness = exp(-(temp_t-t_opt)^2/(2*t_sd^2)) * exp(-t_sd^2/(2*α^2))
    return fitness
end

function hfit(h_opt,h_sd,habitat,α)
    fitness = exp(-(habitat-h_opt)^2/(2*h_sd^2)) * exp(-h_sd^2/(2*α^2))
    return fitness
end

function fitness(t_opt,t_sd,h_opt,h_sd, temp_t, habitat, α)
    fitness =  exp(-(temp_t-t_opt)^2/(2*t_sd^2)) * exp(-(habitat-h_opt)^2/(2*h_sd^2)) * exp(-t_sd^2/(2*α^2)) * exp(-h_sd^2/(2*α^2))
    return fitness
end

function shannon(i)
    p = i ./ sum(i[1:length(i)])
    h = -sum(p.*log.(p))
    return h
end

function shannon(i)
    p = i ./ sum(i[1:length(i)])
    λ =  sum(p.^2)
    return λ
end

#--------------------------------------------------------
# Main script

data = DataFrame(CSV.File("/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/test5/S0_G0.5_B50_T50_UT_Ae0.0.txt"))

data.H_t = string.(data.H_t)
data.H_h = string.(data.H_h)

data.H_t = replace!(data.H_t,"NaN" => "Uniform")
data.H_h = replace!(data.H_h,"NaN" => "Uniform")

# Patch population sizes

populations = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(pop=length(data.ID))
end

rich = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(richness=length(unique(data.LineageID)))
end

data.temp_t = data.temp_t.+12.5

data.tfit = tfit.(data.T_opt,data.T_sd,data.temp_t,data.alpha)
data.hfit = hfit.(data.H_opt,data.H_sd,data.habitat,data.alpha)
data.fit = fitness.(data.T_opt,data.H_opt,data.T_sd,data.H_sd,data.temp_t,data.habitat,data.alpha)

ft = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(ft=mean(data.tfit))
end

fh = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(ft=mean(data.hfit))
end

f = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(ft=mean(data.fit))
end

topt = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(T_opt=mean(data.T_opt))
end

tsd = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(T_sd=mean(data.T_sd))
end

hopt = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(H_opt=mean(data.H_opt))
end

hsd = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(H_sd=mean(data.H_sd))
end

displ = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(disp_l=mean(data.disp_l))
end

dispg = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(disp_g=mean(data.disp_g))
end

fert = by(data, [:Timestep,:x,:y,:alpha]) do data
    DataFrame(fert=mean(data.fert_max))
end

div = select(data,[:Timestep,:x,:y,:alpha,:LineageID])
div.tally = 1
div = by(div,[:Timestep,:x,:y,:alpha,:LineageID]) do div
    DataFrame(count=sum(div.tally))
end

shan = by(div,[:Timestep,:x,:y,:alpha]) do div
    DataFrame(shannon=shannon(div.count))
end
