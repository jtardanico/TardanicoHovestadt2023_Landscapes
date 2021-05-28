#-------------------------
# Data Output Functions
#-------------------------

function set_filename(argumentsdict::Dict, par::Dict,clim,grad)
    filename = string("S", par["scen"])
    println(filename)
    filename = string(filename, "_G",par["grad"])

    if argumentsdict["burninperiod"] == true
        filename = string(filename,"_B", par["bmax"])
        println(filename)
    end
    filename = string(filename, "_T", par["tmax"])
    println(filename)

    if par["temp_type"] == "uniform"
        filename = string(filename,"_UT")
        println(filename)
    elseif par["temp_type"] == "clustered"
        filename = string(filename,"_ClT",par["autocor_temp"])
        println(filename)
    elseif par["temp_type"] == "autocorrelated"
        filename = string(filename,"_AT",par["autocor_temp"])
        println(filename)
    else
        filename = string(filename,"_",par["temp_type"])
    end

    if par["env_type"] == "uniform"
        filename = string(filename,"_UT")
        println(filename)
    elseif par["env_type"] == "clustered"
        filename = string(filename,"_ClT",par["autocor_env"])
        println(filename)
    elseif par["env_type"] == "autocorrelated"
        filename = string(filename,"_AT",par["autocor_env"])
        println(filename)
    else
        if par["env_type"] != par["temp_type"]
            filename = string(filename,"_",par["env_type"])
            println(filename)
        else
            println(filename)
        end
    end

    tempsource = argumentsdict["tempsource"]
    envsource = argumentsdict["envsource"]

    n = 1
    outname = string(filename,"_",n)
    searchname = string(filename,"_",n,".txt")

    files = string.(par["dir"],readdir(par["dir"]))
    while in(true,occursin.(Regex(searchname),files))==true
        n = n+1
        searchname = string(filename,"_",(n),".txt")
        outname = string(filename,"_",(n))
    end
    println(outname)
    return outname,tempsource,envsource
end

# Counts the total population in the landscape
function popcount(landscape)
    x = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                x = x + length(landscape[i,j].species[l][1:end,1])
            end
        end
    end
    #println("landscape population: $x")
    return x
end

function fitnessmeans(landscape,trend,meantrend,α)
    obs = 0
    sumft = 0
    sumfh = 0
    sumfit = 0
    sumtdiff = 0
    sumavgtdiff = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                obs=obs+length(landscape[i,j].species[l][1:end,1])
                if size(landscape[i,j].species[l])[1]>0
                    for p in 1:length(landscape[i,j].species[l][1:end,1])
                        ft,fh,fit = expected_fert2.(landscape[i,j].species[l][p,2],
                                       landscape[i,j].species[l][p,3],
                                       landscape[i,j].species[l][p,4],
                                       landscape[i,j].species[l][p,5],
                                       landscape[i,j].temp_t,
                                       landscape[i,j].habitat,
                                       α,α,trend)
                        tdiff = landscape[i,j].species[l][p,2] .- (landscape[i,j].temp_t .+ trend)
                        avg_tdiff = landscape[i,j].species[l][p,2] .- (landscape[i,j].temp_t .+ meantrend)
                        sumtdiff = sumtdiff + sum(tdiff)
                        sumavgtdiff = sumavgtdiff + sum(avg_tdiff)
                        sumft = sumft + ft
                        sumfh = sumfh + fh
                        sumfit = sumfit + fit
                    end
                end
            end
        end
    end
    meantdiff = sumtdiff/obs
    meanavgtdiff = sumavgtdiff/obs
    meanft = sumft/obs
    meanfh = sumfh/obs
    meanfit = sumfit/obs
    ssqft = 0
    ssqfh = 0
    ssqfit = 0
    ssqtdiff = 0
    ssqavgtdiff = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                if size(landscape[i,j].species[l])[1]>0
                    for p in 1:length(landscape[i,j].species[l][1:end,1])
                        ft,fh,fit = expected_fert2.(landscape[i,j].species[l][p,2],
                                       landscape[i,j].species[l][p,3],
                                       landscape[i,j].species[l][p,4],
                                       landscape[i,j].species[l][p,5],
                                       landscape[i,j].temp_t,
                                       landscape[i,j].habitat,
                                       α,α,trend)
                        tdiff = landscape[i,j].species[l][p,2] .- (landscape[i,j].temp_t .+ trend)
                        avg_tdiff = landscape[i,j].species[l][p,2] .- (landscape[i,j].temp_t .+ meantrend)
                        ssqtdiff = ssqtdiff + sum((tdiff - meantdiff)^2)
                        ssqavgtdiff = ssqavgtdiff + sum((avg_tdiff - meanavgtdiff)^2)
                        ssqft = ssqft + sum((ft - meanft)^2)
                        ssqfh = ssqfh + sum((fh - meanfh)^2)
                        ssqfit = ssqfit + sum((fit - meanfit)^2)
                    end
                end
            end
        end
    end
    vartdiff = ssqtdiff/obs
    varavgtdiff = ssqavgtdiff/obs
    varft = ssqft/obs
    varfh = ssqfh/obs
    varfit = ssqfit/obs

    return meanft,meanfh,meanfit,meantdiff,meanavgtdiff,varft,varfh,varfit,vartdiff,varavgtdiff
end

# Calculates landscape-wide arithmetic means and variance for traits
function traitmeans(landscape)
    obs = 0
    sumtopt = 0
    sumtsd = 0
    sumhopt = 0
    sumhsd = 0
    sumdisp = 0
    sumdispg = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                obs=obs+length(landscape[i,j].species[l][1:end,1])
                sumtopt = sumtopt + sum(landscape[i,j].species[l][1:end,2])
                sumtsd = sumtsd + sum(landscape[i,j].species[l][1:end,3])
                sumhopt = sumhopt + sum(landscape[i,j].species[l][1:end,4])
                sumhsd = sumhsd + sum(landscape[i,j].species[l][1:end,5])
                sumdisp = sumdisp + sum(landscape[i,j].species[l][1:end,6])
                sumdispg = sumdispg + sum(landscape[i,j].species[l][1:end,7])
            end
        end
    end
    meantopt = sumtopt/obs
    meantsd = sumtsd/obs
    meanhopt = sumhopt/obs
    meanhsd = sumhsd/obs
    meandisp = sumdisp/obs
    meandispg = sumdispg/obs
    ssqtopt = 0
    ssqtsd = 0
    ssqhopt = 0
    ssqhsd = 0
    ssqdisp = 0
    ssqdispg = 0
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                ssqtopt = ssqtopt + sum((landscape[i,j].species[l][1:end,2] .- meantopt).^2)
                ssqtsd = ssqtsd + sum((landscape[i,j].species[l][1:end,3] .- meantsd).^2)
                ssqhopt = ssqhopt + sum((landscape[i,j].species[l][1:end,4] .- meanhopt).^2)
                ssqhsd = ssqhsd + sum((landscape[i,j].species[l][1:end,5] .- meanhsd).^2)
                ssqdisp = ssqdisp + sum((landscape[i,j].species[l][1:end,6] .- meandisp).^2)
                ssqdispg = ssqdispg + sum((landscape[i,j].species[l][1:end,7] .- meandispg).^2)
            end
        end
    end
    vartopt = ssqtopt/obs
    vartsd = ssqtsd/obs
    varhopt = ssqhopt/obs
    varhsd = ssqhsd/obs
    vardisp = ssqdisp/obs
    vardispg = ssqdispg/obs
    return meantopt,meantsd,meanhopt,meanhsd,meandisp,meandispg,vartopt,vartsd,varhopt,varhsd,vardisp,vardispg
end

# Calculates species richness for a single patch
function landscape_div(landscape::Array{TPatch,2})
    lineages = Array{String,1}(undef,0)
    for i in 1:length(landscape[1:end,1])
        for j in 1:length(landscape[1,1:end])
            for l in length(landscape[i,j].species[1:end])
                lineages = [lineages;landscape[i,j].species[l][1:end,10]]
            end
        end
    end

    rich = length(unique(lineages))
    if rich > 0
        div = DataFrame()
        div.lineages = lineages
        div.tally = ones(length(div.lineages))
        div = combine(groupby(div,[:lineages])) do div
            DataFrame(count=sum(div.tally))
        end
        simp = simpson(div.count)
        shan = shannon(div.count)
        return rich,simp,shan
    else
        rich = 0
        simp = "NA"
        shan = "NA"
        return rich,simp,shan
    end
end

# Counts population for each species in each patch and outputs them to an array with dimensions
# rows by cols by n_species.
function patch_spp_pops(landscape::Array{TPatch,2},n_species)
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    patch_pops = Array{Array{Int,2},1}(undef, n_species)
    for p in 1:n_species
        pops = Array{Int,2}(undef,rows,cols)
        for i in 1:rows
            for j in 1:cols
                pops[i,j] = length(landscape[i,j].species[p][1:end,1])
            end
        end
        patch_pops[p] = pops
    end
    return patch_pops
end

# Calculates trait value means for each patch in the landscape and outputs a 3d array with dimensions
# rows by cols by n_traits.
# Note: This function could be more elegantly written. Consider rewriting this function.
function trait_analysis(landscape::Array{TPatch,2})
    #println("Trait Analysis")
    rows = length(landscape[1:end,1])
    #println("rows = $rows")
    cols = length(landscape[1,1:end])
    #println("cols = $cols")
    n_species = length(species_list[1:end,1])
    n_traits = length(species_list[1,2:8]) # 7 traits of interest,
    global trait_means = Array{Array{Float32,2}}(undef,n_traits)
    for k in 1:n_traits
        #println("trait loop = $k")
        means = Array{Float32,1}(undef,n_species)
        weights = Array{Int,1}(undef,n_species)
        mean_traits_values = Array{Float32,2}(undef,rows,cols)
        for i in 1:rows
            #println("Row = $i")
            for j in 1:cols
                if n_species > 1
                    #println("Col = $j")
                    for p in 1:n_species
                        #println("Species = $p")
                        weights[p] = length(landscape[i,j].species[p][1:end,1])
                        #println("Weight = $(length(landscape[i,j].species[p][1:end,1]))")
                        if length(landscape[i,j].species[p][1:end,1]) > 0
                            means[p] = mean(landscape[i,j].species[p][1:end,k+1])
                            #println("Mean = $(means[p])")
                        else
                            means[p] = 0
                            #println("Mean = $(means[p])")
                        end
                    end # End species loop
                end
                #println("Weighted mean = $(weightedmean(means,weights))")
                if n_species > 1
                    mean_traits_values[i,j] = weightedmean(means,weights)
                else
                    mean_traits_values[i,j] = mean(landscape[i,j].species[1][1:end,k+1])
                end
            end # End col loop
        end # End row loop
        global trait_means[k] = mean_traits_values
    end # End traits loop
end

function diversity_analysis(landscape::Array{TPatch,2})
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    n_species = length(species_list[1:end,1])
    global pops = Array{Array{Int,2},1}(undef, n_species)
    global richness = Array{Int,2}(undef,rows,cols)
    #global mean_traits_values = Array{Array{Float32,2},1}(undef,6)
    for i in 1:rows
        for j in 1:cols
            richness[i,j] = patch_species_richness(landscape,i,j)
        end
    end
    global pops = patch_spp_pops(landscape,n_species)
end

function env_analysis(landscape::Array{TPatch,2},trend_t)
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    temperatures = Array{Float32,2}(undef,rows,cols)
    habitat = Array{Float32,2}(undef,rows,cols)
    for i in 1:rows
        for j in 1:cols
            temperatures[i,j] = landscape[i,j].temp_t + trend_t
            habitat[i,j] = landscape[i,j].habitat
        end
    end
    global temperatures = temperatures
    global habitat = habitat
end

# Calculates unweighted mean stress for each patch
function mean_stress(landscape::Array{TPatch,2}, trend)
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    n_species = length(species_list[1:end,1])
    stress_t = Array{Float32,2}(undef,rows,cols) # Temperature stress
    stress_h = Array{Float32,2}(undef,rows,cols) # Habitat stress
    stress_o = Array{Float32,2}(undef,rows,cols) # Overall stress
    for i in 1:rows
        for j in 1:cols
            for p in 1:n_species
                S_T, S_H = stress.(landscape[i,j].species[p][1:end,2], # Calculating environmental stress for a species in a patch
                           landscape[i,j].species[p][1:end,3],
                           landscape[i,j].species[p][1:end,4],
                           landscape[i,j].species[p][1:end,5],
                           landscape[i,j].temp_t,
                           landscape[i,j].habitat, trend)
                S_O = (S_T .* S_H) # Calculate overall stress
                stress_t[i,j] = mean(S_T)
                stress_h[i,j] = mean(S_H)
                stress_o[i,j] = mean(S_O)
            end
        end
    end
    return stress_t,stress_h,stress_o
end

function write_landscape_stats(landscape,directory,filename,replicate,timestep,s_clim,trend,mean_trend,grad,H_t,H_h,α,bmax)
    outputname = string(directory,filename,"trend.txt")
    col_names1 = ["Replicate" "Timestep" "Pop" "rich" "simp" "shan" "T_opt" "T_sd" "H_opt" "H_sd" "disp" "disp_g" "VT_opt" "VT_sd" "VH_opt" "VH_sd" "Vdisp" "Vdisp_g"]
    col_names2 = ["ft" "fh" "fit" "tdiff" "avgtdiff" "varft" "varfh" "varfit" "var_tdiff" "var_avgtdiff" "clim_scen" "trend" "mean_trend" "grad" "H_t" "H_h" "alpha"]
    col_names = hcat(col_names1,col_names2)
    T_opt, T_sd, H_opt, H_sd, disp, disp_g, VT_opt, VT_sd, VH_opt, VH_sd, Vdisp, Vdisp_g = traitmeans(landscape)
    ft, fh, fit,tdiff,avgtdiff, varft, varfh, varfit, vartdiff, varavgtdiff = fitnessmeans(landscape,trend,mean_trend,α)
    rich, simp, shan = landscape_div(landscape)
    pop = popcount(landscape)
    open(outputname,"a") do IO
        if timestep==-1 && replicate==1
            writedlm(IO,col_names)
        end
        # Change to hcat
        writedlm(IO, [replicate timestep pop rich simp shan T_opt T_sd H_opt H_sd disp disp_g VT_opt VT_sd VH_opt VH_sd Vdisp Vdisp_g ft fh fit tdiff avgtdiff varft varfh varfit vartdiff varavgtdiff s_clim trend mean_trend grad H_t H_h α])
    end
end

# Writes the full landscape data to a .csv file, including every individual in every patch.
# Caution: Very large data file.

# OPTIMIZATION: WRITE DATA TO ARRAY, THEN WRITE ARRAY TO FILE <---- DO THIS!!!

function write_landscape_csv(landscape,directory,filename,replicate,timestep,s_clim,trend,mean_trend,grad,H_t,H_h,α)
    outputname = string(directory,filename,".txt")
    rows = length(landscape[1:end,1])
    cols = length(landscape[1,1:end])
    n_species = length(landscape[1,1].species[1:end,1])
    col_names = ["ID" "Replicate" "Timestep" "H_t" "H_h" "alpha" "clim_scen" "gradient" "x" "y" "T_opt" "T_sd" "H_opt" "H_sd" "disp_l" "disp_g" "fert_max" "LineageID" "temp_t" "trend" "mean_trend" "precip_t" "habitat"]
    open(outputname, "a") do IO
        if timestep==-1 && replicate==1
            writedlm(IO, col_names)
        end
        for i in 1:rows
            for j in 1:cols
                for k in 1:n_species
                    for l in 1:length(landscape[i,j].species[k][1:end,1])
                        #println("lin ID = ",string(landscape[i,j].species[k][l,10]))
                        lin = parse(Int,(@sprintf("%.0f",landscape[i,j].species[k][l,10]*10^15)))
                        lineageID = string(lin, base=62)
                        # Change to hcat
                        writedlm(IO, [k replicate timestep H_t H_h α s_clim grad i j landscape[i,j].species[k][l,2] landscape[i,j].species[k][l,3] landscape[i,j].species[k][l,4] landscape[i,j].species[k][l,5] landscape[i,j].species[k][l,6] landscape[i,j].species[k][l,7] landscape[i,j].species[k][l,8] lineageID landscape[i,j].temp_t trend mean_trend landscape[i,j].precip_t landscape[i,j].habitat])
                    end
                end
            end
        end
        close(IO)
    end
end

# Creates heatmap graphs with patch values for environmental variables, richness, trait means,
# and individual species populations in .png format.
function heatmaps(richness,temperatures,habitats,pops,trait_means,filename,directory)
    type = ".png"
    rich = string(directory,"Rich_",filename,type)
    println(rich)
    temps = string(directory,"Temp_",filename,type)
    println(temps)
    hab = string(directory,"Habitat_",filename,type)
    println(hab)
    heatmap(richness)
    png(rich)
    heatmap(temperatures)
    png(temps)
    heatmap(habitats)
    png(hab)
    for i in 1:length(species_list[1:end,1])
        spdist = string("Spp_$i","_")
        spfilename = string(directory,spdist,filename,type)
        println(spfilename)
        heatmap(pops[i])
        png(spfilename)
    end
    for i in 1:length(trait_means)
        traits = string("Trait_$i","_")
        trfilename = string(directory,traits,filename,type)
        println(trfilename)
        heatmap(trait_means[i])
        png(trfilename)
    end
end
