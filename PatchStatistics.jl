# Landscape statistics calculation

#----------------------------------------------------------------------
# This is a Julia implementation of the R script LandscapePatchStat2.R
#----------------------------------------------------------------------
using DataFrames
using CSV
using StatsBase



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

function simpson(i)
    p = i ./ sum(i[1:length(i)])
    λ =  sum(p.^2)
    return λ
end

function calculations(data)
    data = DataFrame(CSV.File(data))

    data.H_t = string.(data.H_t)
    data.H_h = string.(data.H_h)

    data.H_t = replace!(data.H_t,"NaN" => "Uniform")
    data.H_h = replace!(data.H_h,"NaN" => "Uniform")

    # Patch population sizes

    populations = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(pop=length(data.ID))
    end

    rich = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(richness=length(unique(data.LineageID)))
    end

    data.temp_t = data.temp_t.+12.5

    data.tfit = tfit.(data.T_opt,data.T_sd,data.temp_t,data.alpha)
    data.hfit = hfit.(data.H_opt,data.H_sd,data.habitat,data.alpha)
    data.fit = fitness.(data.T_opt,data.H_opt,data.T_sd,data.H_sd,data.temp_t,data.habitat,data.alpha)

    ft = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(ft=mean(data.tfit))
    end

    fh = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(fh=mean(data.hfit))
    end

    f = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(f=mean(data.fit))
    end

    topt = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(T_opt=mean(data.T_opt))
    end

    tsd = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(T_sd=mean(data.T_sd))
    end

    hopt = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(H_opt=mean(data.H_opt))
    end

    hsd = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(H_sd=mean(data.H_sd))
    end

    displ = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(disp_l=mean(data.disp_l))
    end

    dispg = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(disp_g=mean(data.disp_g))
    end

    fert = by(data, [:x,:y,:Timestep,:H_t,:H_h,:alpha]) do data
        DataFrame(fert=mean(data.fert_max))
    end

    div = select(data,[:x,:y,:Timestep,:H_t,:H_h,:alpha,:LineageID])
    div.tally = 1
    div = by(div,[:x,:y,:Timestep,:H_t,:H_h,:alpha,:LineageID]) do div
        DataFrame(count=sum(div.tally))
    end

    shan = by(div,[:x,:y,:Timestep,:H_t,:H_h,:alpha]) do div
        DataFrame(shannon=shannon(div.count))
    end

    simp = by(div,[:x,:y,:Timestep,:H_t,:H_h,:alpha]) do div
        DataFrame(simpson=simpson(div.count))
    end

    out = innerjoin(populations,rich, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,shan, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,simp, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,ft, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,fh, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,f, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,topt, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,tsd, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,hopt, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,hsd, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,displ, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,dispg, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    out = innerjoin(out,fert, on = [:x,:y,:Timestep,:H_t,:H_h,:alpha])
    return out
end

#--------------------------------------------------------
# Main script
@time begin
    inputdir = "/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/test12/"
    outputdir = "/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/G1Patches/"
    outfilename = "patchstats_h1_test12.txt"
    infiles = string.(inputdir,readdir(inputdir))
    outfile = string(outputdir,outfilename)
    for i in 1:length(infiles)
        patchstats = calculations(infiles[i])
        if i==1
            CSV.write(outfile,patchstats,append=true)
        else
            CSV.write(outfile,patchstats,append=false)
        end
    end
end
