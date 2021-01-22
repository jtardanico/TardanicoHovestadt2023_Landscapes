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

    println(data)
    data = DataFrame(CSV.File(data))


    data.H_t = string.(data.H_t)
    data.H_h = string.(data.H_h)

    data.H_t = replace!(data.H_t,"NaN" => "Uniform")
    data.H_h = replace!(data.H_h,"NaN" => "Uniform")

        # Patch population sizes
    populations = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
            DataFrame(pop=length(data.ID))
    end

    data = innerjoin(populations,data, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])

    populations = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
            DataFrame(pop=mean(data.pop))
    end

    rich = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(richness=length(unique(data.LineageID)))
    end

    data.tdiff = data.T_opt .- data.temp_t
    data.mt_tdiff = data.T_opt .- data.mean_trend

    data.time_tfit = tfit.(data.T_opt,data.T_sd,data.temp_t,data.alpha)

    data.mean_tfit = tfit.(data.T_opt,data.T_sd,data.mean_trend,data.alpha)

    data.hfit = hfit.(data.H_opt,data.H_sd,data.habitat,data.alpha)

    data.time_fit = fitness.(data.T_opt,data.H_opt,data.T_sd,data.H_sd,data.temp_t,data.habitat,data.alpha)

    data.trend_fit = fitness.(data.T_opt,data.H_opt,data.T_sd,data.H_sd,data.trend,data.habitat,data.alpha)

    data.mean_fit = fitness.(data.T_opt,data.H_opt,data.T_sd,data.H_sd,data.mean_trend,data.habitat,data.alpha)

    tdiff = by(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(t_ft=mean(data.tdiff))
    end

    mt_tdiff = by(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(t_ft=mean(data.mt_tdiff))
    end

    # time_tfit
    t_ft = by(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(t_ft=mean(data.time_tfit))
    end

    # mean_tfit
    m_ft = by(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(m_ft=mean(data.mean_tfit))
    end

    # hfit
    t_fh = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(t_fh=mean(data.hfit))
    end

    # time_fit
    #t_f = by(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
    #    DataFrame(t_f=mean(data.time_fit))
    #end

    # trend_fit
    #tr_f = by(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
    #    DataFrame(tr_f=mean(data.trend_fit))
    #end
    # mean_fit
    #m_f = by(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
    #    DataFrame(m_f=mean(data.mean_fit))
    #end

    #f = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
    #    DataFrame(f=mean(data.fit))
    #end

    topt = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(T_opt=mean(data.T_opt))
    end

    tsd = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(T_sd=mean(data.T_sd))
    end

    hopt = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(H_opt=mean(data.H_opt))
    end

    hsd = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(H_sd=mean(data.H_sd))
    end

    displ = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(disp_l=mean(data.disp_l))
    end

    dispg = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
        DataFrame(disp_g=mean(data.disp_g))
    end

    fert = by(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient]) do data
            DataFrame(fert=mean(data.fert_max))
    end

    div = select(data,[:x,:y,:Replicate,:Timestep,:LineageID])
    div.tally = 1
    div = by(div,[:x,:y,:Replicate,:Timestep,:LineageID]) do div
        DataFrame(count=sum(div.tally))
    end

    shan = by(div,[:x,:y,:Replicate,:Timestep]) do div
        DataFrame(shannon=shannon(div.count))
    end

    simp = by(div,[:x,:y,:Replicate,:Timestep]) do div
        DataFrame(simpson=simpson(div.count))
    end

    ht = data[1,:].H_t
    hh = data[1,:].H_h
    alpha = data[1,:].alpha
    cs = data[1,:].clim_scen
    grad = data[1,:].gradient

    out = innerjoin(populations,rich, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,shan, on = [:x,:y,:Replicate,:Timestep])
    out = innerjoin(out,simp, on = [:x,:y,:Replicate,:Timestep])
    out = innerjoin(out,tdiff, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,mt_tdiff, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,t_ft, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,m_ft, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,t_fh, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out.time_fit = out.t_ft*out.t_fh
    out.mean_fit = out.m_ft*out.t_fh
    #out = innerjoin(out,f, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,topt, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,tsd, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,hopt, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,hsd, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,displ, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,dispg, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,fert, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out.H_t = ht
    out.H_h = hh
    out.alpha = alpha
    out.clim_scen = cs
    out.gradient = grad
    return out
end

function patchstats(dir::String,filename::String)
    outfile = string(inputdir,filename,"patches.txt")
    infiles = string.(inputdir,readdir(inputdir))
    filter!(x -> occursin(filename,x)==true, infiles)
    filter!(x -> occursin("trend",x)==false, infiles)
    for i in 1:length(infiles)
        patchstats = calculations(infiles[i])
        if i==1
            CSV.write(outfile,patchstats,append=false)
        else
            CSV.write(outfile,patchstats,append=true)
        end
    end
end

#--------------------------------------------------------
# Main script
#@time begin
#    inputdir = "/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/HPCtest1/"
#    outputdir = "/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/debug_testing_stats/"
#    outfilename = "patchstats_HPCtest1.txt"
#    infiles = string.(inputdir,readdir(inputdir))
#    filter!(x -> occursin("trend",x)==false, infiles)
#    outfile = string(outputdir,outfilename)
#    for i in 1:length(infiles)
#        patchstats = calculations(infiles[i])
#        if i==1
#            CSV.write(outfile,patchstats,append=false)
#        else
#            CSV.write(outfile,patchstats,append=true)
#        end
#    end
#end
