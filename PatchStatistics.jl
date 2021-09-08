# Landscape statistics calculation

#----------------------------------------------------------------------
# This is a Julia implementation of the R script LandscapePatchStat2.R
#----------------------------------------------------------------------
using DataFrames
using CSV
using StatsBase
using DelimitedFiles



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

# Creates a data frame of length x*y*timesteps*replicates containing simulation parameters and baseline patch attributes for all patches.
function patch_attributes(tempsource,envsource,data)
    temp_landscape = readdlm(tempsource,';')
    env_landscape = readdlm(envsource,';')

    times = unique(data.Timestep)
    reps = unique(data.Replicate)
    ht = unique(data.H_t)
    hh = unique(data.H_h)
    α = unique(data.alpha)
    scen = unique(data.clim_scen)
    grad = unique(data.gradient)

    trends = Array{Float64,1}(undef,(length(times)*length(reps)))
    counter = 1
    for i in 1:length(reps)
        for j in 1:length(times)
             t = unique(data.trend[.&(data[!,:Replicate].==reps[i], data[!,:Timestep] .== times[j]),:])
             println(t)
             trends[counter] = t[1]
            counter = counter + 1
        end
    end

    n_patches = length(temp_landscape[1:end,1]) * length(temp_landscape[1,1:end]) * length(times)
    v_x = Array{Int,1}(undef,n_patches)
    v_y = Array{Int,1}(undef,n_patches)
    v_time = Array{Int,1}(undef,n_patches)
    v_reps = Array{Float64,1}(undef,n_patches)
    v_ht = Array{Any,1}(undef,n_patches)
    v_hh = Array{Any,1}(undef,n_patches)
    v_α = Array{Float64,1}(undef,n_patches)
    v_scen = Array{Float64,1}(undef,n_patches)
    v_grad = Array{Float64,1}(undef,n_patches)
    v_temp = Array{Float64,1}(undef,n_patches)
    v_env = Array{Float64,1}(undef,n_patches)
    v_trend = Array{Float64,1}(undef,n_patches)
    counter = 1
    t_counter = 1
    for g in 1:length(reps)
        for h in 1:length(times)

            for i in 1:length(temp_landscape[1:end,1])
                for j in 1:length(temp_landscape[1,1:end])
                    v_x[counter] = i
                    v_y[counter] = j
                    v_time[counter] = times[h]
                    v_reps[counter] = reps[g]
                    v_temp[counter] = temp_landscape[i,j]
                    v_env[counter] = env_landscape[i,j]
                    v_trend[counter] = trends[t_counter]
                    counter = counter + 1
                end
            end
            t_counter = t_counter + 1
        end
    end

    v_ht .= ht[1]
    v_hh .= hh[1]
    v_α .= α[1]
    v_scen .= scen[1]
    v_grad .= grad[1]

    attributes = DataFrame(x = v_x,
                           y = v_y,
                           Timestep = v_time,
                           Replicate = v_reps,
                           H_t = v_ht,
                           H_h = v_hh,
                           alpha = v_α,
                           clim_scen = v_scen,
                           gradient = v_grad,
                           temperature = v_temp,
                           habitat = v_env,
                           trend = v_trend)
    return attributes
end

function calculations(data,tempsource,envsource)

    println(data)
    data = DataFrame(CSV.File(data))

    atts = patch_attributes(tempsource,envsource,data)

    data.H_t = string.(data.H_t)
    data.H_h = string.(data.H_h)

    atts.H_t = string.(atts.H_t)
    atts.H_h = string.(atts.H_h)

    data.H_t = replace!(data.H_t,"NaN" => "Uniform")
    data.H_h = replace!(data.H_h,"NaN" => "Uniform")

    atts.H_t = replace!(atts.H_t,"NaN" => "Uniform")
    atts.H_h = replace!(atts.H_h,"NaN" => "Uniform")

    # Patch population sizes
    populations = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
            DataFrame(pop=length(data.ID))
    end

    data = innerjoin(populations,data, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])

    rich = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(richness=length(unique(data.LineageID)))
    end

    #temperature = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
    #    DataFrame(temperature=unique(data.temp_t))
    #end

    #habitat = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
    #    DataFrame(habitat=unique(data.habitat))
    #end

    data.trend_t = data.temp_t .+ data.trend
    data.mean_trend_t = data.temp_t .+ data.mean_trend

    data.tdiff = data.T_opt .- data.trend_t
    data.mt_tdiff = data.T_opt .- data.mean_trend_t

    data.time_tfit = tfit.(data.T_opt,data.T_sd,data.trend_t,data.alpha)

    data.mean_tfit = tfit.(data.T_opt,data.T_sd,data.mean_trend_t,data.alpha)

    data.hfit = hfit.(data.H_opt,data.H_sd,data.habitat,data.alpha)

    data.time_fit = fitness.(data.T_opt,data.H_opt,data.T_sd,data.H_sd,data.trend_t,data.habitat,data.alpha)

    data.trend_fit = fitness.(data.T_opt,data.H_opt,data.T_sd,data.H_sd,data.trend_t,data.habitat,data.alpha)

    data.mean_fit = fitness.(data.T_opt,data.H_opt,data.T_sd,data.H_sd,data.mean_trend_t,data.habitat,data.alpha)

    tdiff = combine(groupby(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(tdiff=mean(data.tdiff))
    end

    mt_tdiff = combine(groupby(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(mt_tdiff=mean(data.mt_tdiff))
    end

    # time_tfit
    t_ft = combine(groupby(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(t_ft=mean(data.time_tfit))
    end

    # mean_tfit
    m_ft = combine(groupby(data, [:x,:y,:Replicate, :Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(m_ft=mean(data.mean_tfit))
    end

    # hfit
    t_fh = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
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

    topt = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(T_opt=mean(data.T_opt))
    end

    tsd = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(T_sd=mean(data.T_sd))
    end

    hopt = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(H_opt=mean(data.H_opt))
    end

    hsd = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(H_sd=mean(data.H_sd))
    end

    displ = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(disp_l=mean(data.disp_l))
    end

    dispg = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
        DataFrame(disp_g=mean(data.disp_g))
    end

    fert = combine(groupby(data, [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])) do data
            DataFrame(fert=mean(data.fert_max))
    end

    div = select(data,[:x,:y,:Replicate,:Timestep,:LineageID])
    div.tally = ones(length(div.x))
    div = combine(groupby(div,[:x,:y,:Replicate,:Timestep,:LineageID])) do div
        DataFrame(count=sum(div.tally))
    end

    shan = combine(groupby(div,[:x,:y,:Replicate,:Timestep])) do div
        DataFrame(shannon=shannon(div.count))
    end

    simp = combine(groupby(div,[:x,:y,:Replicate,:Timestep])) do div
        DataFrame(simpson=simpson(div.count))
    end

    ht = data[1,:].H_t
    hh = data[1,:].H_h
    alpha = data[1,:].alpha
    cs = data[1,:].clim_scen
    grad = data[1,:].gradient

    out = innerjoin(populations,rich, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    #out = innerjoin(out,temperature, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    #out = innerjoin(out,habitat, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,shan, on = [:x,:y,:Replicate,:Timestep])
    out = innerjoin(out,simp, on = [:x,:y,:Replicate,:Timestep])
    out = innerjoin(out,tdiff, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,mt_tdiff, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,t_ft, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,m_ft, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,t_fh, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out.time_fit = out.t_ft.*out.t_fh
    out.mean_fit = out.m_ft.*out.t_fh
    #out = innerjoin(out,f, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,topt, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,tsd, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,hopt, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,hsd, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,displ, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,dispg, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out = innerjoin(out,fert, on = [:x,:y,:Replicate,:Timestep,:H_t,:H_h,:alpha,:clim_scen,:gradient])
    out.H_t .= ht
    out.H_h .= hh
    out.alpha .= alpha
    out.clim_scen .= cs
    out.gradient .= grad

    out = outerjoin(out,atts, on = [:x,:y,:Timestep,:Replicate,:H_t,:H_h,:alpha,:clim_scen,:gradient]) # Outer join with 'atts' ensures that empty patches have data entries.

    return out
end

function patchstats(dir::String,filename::String,tempsource::String,envsource::String)
    println("Calculating patch stats...")
    outfile = string(dir,filename,"patches.txt")
    println("Output file: $outfile")
    infiles = string.(dir,readdir(dir))
    n_files = length(infiles)
    println("No. files in directory = $n_files")
    println("Files in directory:")
    for i in 1:length(infiles)
        println(infiles[i])
    end
    println("")
    filter!(x -> occursin(filename,x)==true, infiles)
    n_files = length(infiles)
    println("No. files = $n_files") # <----------------
    println("Files containing $filename:")
    for i in 1:length(infiles)
        println(infiles[i])
    end
    println("")
    filter!(x -> occursin("trend",x)==false, infiles)
    n_files = length(infiles)
    println("No. files = $n_files")
    println("Files containing $filename -'trend': ")
    for i in 1:length(infiles)
        println(infiles[i])
    end
    println("")
    for i in 1:length(infiles)
        patchstats = calculations(infiles[i],tempsource,envsource)
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
