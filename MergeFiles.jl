# Data file merge script

# For use merging simulation replicate output files.
# Group files by simulation parameters.

println("Loading packages")

using DataFrames
using CSV
using ArgParse

#------------------------------------
# Function definitions
#------------------------------------

function read_args()
    s=ArgParseSettings()
    @add_arg_table s begin
        "--parasource","-n"
            help = "parameter dictionary source file"
    end
    return parse_args(s)
end

#----------------------------------
# Main Script
#----------------------------------

println("Starting main script")
args = read_args()
println("Parasource = ",args["parasource"])
par = include(args["parasource"])
dir = par["dir"]
println("File directory = ",dir)

#---------------------------------------
# Use if running directly from IDE
#---------------------------------------

#using DataFrames
#using CSV
#using ArgParse

#dir = ""


#---------------------------------------

infiles = string.(dir,readdir(dir))

filekeys = copy(infiles)

s = r"_\d+\.txt"
s2 = r"_\d+trend\.txt"
s3 = r"_\d+patches\.txt"

println("Getting file name keywords")
for i in 1:length(filekeys)
    filekeys[i] = replace(filekeys[i],dir=>"")
    filekeys[i] = replace(filekeys[i],s=>"")
    filekeys[i] = replace(filekeys[i],s2=>"")
    filekeys[i] = replace(filekeys[i],s3=>"")
end

unique!(filekeys)
println(filekeys)

println("Merging data files")
for i in 1:length(filekeys)
    # Merge whole landscape data files
    println(filekeys[i])
    infiles2 = infiles[occursin.(Regex("trend"),infiles).==false]
    infiles2 = infiles2[occursin.(Regex("patches"),infiles2).==false]
    infiles2 = infiles2[occursin.(Regex(filekeys[i]),infiles2).==true]
    outfile = string(dir,"merged",filekeys[i],".txt")
    println(length(infiles2))
    #println(outfile)
    for j in 1:length(infiles2)
        println("file $j")
        file = DataFrame(CSV.File(infiles2[j]))
        file.Replicate .= j
        if j==1

            CSV.write(outfile,file,append=false)
        else
            CSV.write(outfile,file,append=true)
        end
    end
    # Merge trend data files
    #println(filekeys[i]," trend")
    infiles2 = infiles[occursin.(Regex("trend"),infiles).==true]
    infiles2 = infiles2[occursin.(Regex("patches"),infiles2).==false]
    infiles2 = infiles2[occursin.(Regex(filekeys[i]),infiles2).==true]
    outfile = string(dir,"merged",filekeys[i],"_trend.txt")
    println(length(infiles2))
    #println(outfile)
    for j in 1:length(infiles2)
        file = DataFrame(CSV.File(infiles2[j]))
        file.Replicate .= j
        if j==1
            CSV.write(outfile,file,append=false)
        else
            CSV.write(outfile,file,append=true)
        end
    end
    #println(filekeys[i]," patches")
    infiles2 = infiles[occursin.(Regex("trend"),infiles).==false]
    infiles2 = infiles2[occursin.(Regex("patches"),infiles2).==true]
    infiles2 = infiles2[occursin.(Regex(filekeys[i]),infiles2).==true]
    outfile = string(dir,"merged",filekeys[i],"_patches.txt")
    println(length(infiles2))
    for j in 1:length(infiles2)
        file = DataFrame(CSV.File(infiles2[j]))
        file.Replicate .= j
        if j==1
            CSV.write(outfile,file,append=false)
        else
            CSV.write(outfile,file,append=true)
        end
    end
    println("end script")
end
