# Landscape Statistics

using DataFrames
using CSV
using StatsBase

data = DataFrame(CSV.File("/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/a/S0_G0.5_B50_T50_AT0.5_UE.txt"))

inputdir = "/home/joseph/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/a/"
infiles = string.(inputdir,readdir(inputdir))
c = Array{Any}(undef,length(infiles))
stats = DataFrame()
stats. = c
